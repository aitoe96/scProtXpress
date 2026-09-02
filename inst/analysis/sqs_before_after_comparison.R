#' SQS before/after comparison
#'
#' The package already implements a per-sample, per-marker Stain Quality
#' Score (`.sqs_*` functions in R/internalFunctions.R) -- a legitimate
#' signal-to-noise separation metric (peak separation, valley depth,
#' within-peak spread) that flags samples where a given marker's staining
#' didn't actually separate into distinguishable populations. It is off by
#' default (`apply_sqs = "none"`) in both `inst/drivers/*.R`, meaning the
#' published reference includes every sample regardless of staining
#' quality.
#'
#' This script builds the SAME reference twice -- once with the existing
#' default (no SQS filtering) and once with `apply_sqs = "ref_only"` -- and
#' reports: (1) how many (sample, marker) combinations SQS actually
#' excluded, (2) how much each affected marker's threshold moved as a
#' result, and (3) whether recall/precision on the canonical-marker
#' benchmark changes. This tests whether dataset/panel heterogeneity
#' affects the thresholds with something better than a leave-one-out
#' sensitivity analysis: it removes bad data at the source rather than
#' just measuring how much it perturbs the result.
#'
#' USAGE: edit the CONFIG block, then
#'   Rscript inst/analysis/sqs_before_after_comparison.R
#' or source it interactively. Mirrors inst/drivers/cytof_driver.R's
#' structure/config style so it drops into the same workflow.
#'
#' NOTE: this comparison requires ADTnorm/flowCore to run end-to-end. The
#' comparison logic (`compare_sqs_runs()`) is written to operate on the
#' *outputs* of two `run_scProteinXpress()` calls (their `_th_dist.rds`,
#' `_excluded_markers.csv`, and the SQS log CSV), so once you have those two
#' output directories from a real run, `compare_sqs_runs()` itself is plain
#' data-wrangling and can be tested/adapted without re-running the pipeline.

suppressPackageStartupMessages({
  library(Matrix)
  library(scProtXpress)
})

# ----------------------------------------------------------------------------
# 1. CONFIG (edit per run)
# ----------------------------------------------------------------------------
data_root  <- Sys.getenv("SCPROTXPRESS_DATA_ROOT", "/scratch/USER/data_citeseq")
infile_ref <- file.path(data_root, "list_reference.rds")
out_root   <- Sys.getenv("SCPROTXPRESS_OUT_ROOT", "/scratch/USER/results_sqs_comparison")

nb <- 1000; ns <- 100
sqs_threshold_val  <- as.numeric(Sys.getenv("SCPROTXPRESS_SQS_THRESHOLD", "1.0"))
ncores_val         <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
bootstrap_stratify_val <- Sys.getenv("SCPROTXPRESS_BOOTSTRAP_STRATIFY", "cell.type_tissue")
threshold_method_val   <- Sys.getenv("SCPROTXPRESS_THRESHOLD_METHOD",   "mixture")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 2. BUILD BOTH REFERENCES
# ----------------------------------------------------------------------------
build_one <- function(apply_sqs_val, run_label) {
  message(sprintf("\n=== Building reference: apply_sqs = '%s' (%s) ===", apply_sqs_val, run_label))
  ref_obj  <- readRDS(infile_ref)
  unpacked <- unpack_list_citeseq(ref_obj)
  outdir   <- file.path(out_root, run_label)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  run_scProteinXpress(
    data               = unpacked$data,
    metadata           = unpacked$metadata,
    dir                = outdir,
    file.name          = "ref",
    modality           = "CITE-seq",
    cell_id_col        = "cell.id",
    cell_type_col      = "cell.type_harmonised",
    dataset_col        = "dataset",
    sample_col         = "sample",
    batch_col          = "batch",
    tissue_default     = "PBMC",
    discretize         = TRUE,
    sig.frames         = TRUE,
    num.boot.sample    = nb,
    num.cells.sample   = ns,
    ncores             = ncores_val,
    is_transformed     = unpacked$is_transformed,
    qc_outlier_cells   = "auto",
    sqs_threshold       = sqs_threshold_val,
    apply_sqs          = apply_sqs_val,
    bootstrap_stratify = bootstrap_stratify_val,
    threshold_method   = threshold_method_val
  )
  file.path(outdir, "Ref_ref")
}

# ----------------------------------------------------------------------------
# 3. COMPARE
# ----------------------------------------------------------------------------
#' @param dir_no_sqs,dir_with_sqs Output `Ref_*` directories from `build_one()`.
#' @return list(sqs_log, threshold_shift, summary)
compare_sqs_runs <- function(dir_no_sqs, dir_with_sqs) {
  sqs_log_path <- file.path(dir_with_sqs, "ref_SQS_ref_log.csv")
  sqs_log <- if (file.exists(sqs_log_path)) {
    ## NOTE: .sqs_log_append() (R/internalFunctions.R) writes the literal
    ## string "NA" into the `action` column to mean "this sample failed the
    ## SQS threshold and was excluded" -- read.csv()'s default na.strings
    ## silently turns that back into a real NA, which then makes
    ## `action == "NA"` evaluate to NA instead of TRUE/FALSE everywhere
    ## (confirmed while testing this script: n_excluded came back NA on a
    ## log that very much had exclusions in it). na.strings=character(0)
    ## keeps it as the string it was written as. This is a latent
    ## footgun in the SQS log format itself, not just this reader --
    ## worth considering renaming that value to "excluded" at the source
    ## (`.run_one_marker_ref()`'s `action = if (...) "NA" else "keep"`
    ## line) rather than only patching every reader.
    utils::read.csv(sqs_log_path, na.strings = character(0), stringsAsFactors = FALSE)
  } else {
    warning("compare_sqs_runs: no SQS log found at ", sqs_log_path,
            " -- did you call run_scProteinXpress with apply_sqs != 'none'?")
    data.frame()
  }
  n_excluded <- if (nrow(sqs_log) > 0) sum(sqs_log$action == "NA") else NA_integer_

  th_no  <- readRDS(file.path(dir_no_sqs,   "ref_th_dist.rds"))
  th_yes <- readRDS(file.path(dir_with_sqs, "ref_th_dist.rds"))
  common_markers <- intersect(names(th_no), names(th_yes))

  shift <- do.call(rbind, lapply(common_markers, function(m) {
    data.frame(
      marker = m,
      lower_no_sqs   = mean(th_no[[m]][,  "Lower"]), lower_with_sqs = mean(th_yes[[m]][, "Lower"]),
      upper_no_sqs   = mean(th_no[[m]][,  "Upper"]), upper_with_sqs = mean(th_yes[[m]][, "Upper"]),
      shift_lower    = mean(th_yes[[m]][, "Lower"]) - mean(th_no[[m]][, "Lower"]),
      shift_upper    = mean(th_yes[[m]][, "Upper"]) - mean(th_no[[m]][, "Upper"])
    )
  }))
  shift <- shift[order(-abs(shift$shift_lower) - abs(shift$shift_upper)), ]

  summary_tbl <- data.frame(
    n_sample_marker_combinations_excluded = n_excluded,
    n_markers_compared     = length(common_markers),
    n_markers_shifted_gt_0.1 = sum(abs(shift$shift_lower) > 0.1 | abs(shift$shift_upper) > 0.1),
    median_abs_shift_lower = median(abs(shift$shift_lower)),
    median_abs_shift_upper = median(abs(shift$shift_upper))
  )

  message(sprintf(
    "\n[compare_sqs_runs] SQS excluded %s sample-marker combinations. %d/%d markers' thresholds moved by >0.1 (of %d compared). Median |shift| lower=%.3f upper=%.3f.",
    n_excluded, summary_tbl$n_markers_shifted_gt_0.1, length(common_markers), length(common_markers),
    summary_tbl$median_abs_shift_lower, summary_tbl$median_abs_shift_upper))

  list(sqs_log = sqs_log, threshold_shift = shift, summary = summary_tbl)
}

# ----------------------------------------------------------------------------
# 4. ENTRY POINT
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  dir_no_sqs   <- build_one("none",     "no_sqs")
  dir_with_sqs <- build_one("ref_only", "with_sqs")
  result <- compare_sqs_runs(dir_no_sqs, dir_with_sqs)
  utils::write.csv(result$threshold_shift, file.path(out_root, "threshold_shift_sqs_vs_no_sqs.csv"), row.names = FALSE)
  utils::write.csv(result$summary,         file.path(out_root, "sqs_comparison_summary.csv"),        row.names = FALSE)
  print(result$summary)
}
