#' CyTOF driver template
#'
#' Single, modular, sequential driver for the scProtXpress pipeline against
#' a CyTOF SingleCellExperiment. Mirrors the production
#' `run_scProtDefaultParams.R` driver used in Hyperion.
#'
#' Phases:
#'   1. Build reference background from a CyTOF SCE (arcsinh `exprs` assay).
#'   2. Query external datasets against the cached Ref.
#'   3. Held-out query: cells dropped by hierarchical subsampling during
#'      Ref construction.
#'
#' Usage (interactive):
#'   library(scProtXpress)
#'   source(system.file("drivers", "cytof_driver.R", package = "scProtXpress"))
#'
#' SLURM template lives at `inst/slurm/scProtXpress.sh`.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(scProtXpress)
})

# ----------------------------------------------------------------------------
# 1. CONFIG (edit per run)
# ----------------------------------------------------------------------------
data_root   <- Sys.getenv("SCPROTXPRESS_DATA_ROOT", "/scratch/USER/data")
infile_ref  <- file.path(data_root, "sce_reference.rds")
ds_query    <- file.path(data_root, "sce_query.rds")
out_root    <- Sys.getenv("SCPROTXPRESS_OUT_ROOT",  "/scratch/USER/results")
run_name    <- Sys.getenv("SCPROTXPRESS_RUN_NAME", "cytof_run")
outdir      <- file.path(out_root, run_name)
file_name_ref <- "CyTOF_ADTNorm"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

p_assay           <- "exprs"          # "exprs" => arcsinh-transformed
nb                <- 1000             # bootstrap iterations
ns                <- 1000             # cells per bootstrap sample (CyTOF default)
sqs_threshold_val <- 1.0
enable_sqs_report <- FALSE

## --- THESE THREE USED TO BE IMPLICIT (not passed to run_scProteinXpress at
## all, silently resolving to the package defaults "none"/"none"/
## "rectangle"). Made explicit here after finding that the manuscript's
## Methods describes a cell-type-by-tissue-stratified bootstrap and an
## SQS-filtered reference, neither of which this driver was actually
## requesting -- see the same note in `citeseq_driver.R`. If you change
## these, the `_provenance.json` written next to every reference records
## exactly what ran; check it before assuming the manuscript text and a
## given `Ref_*` output agree. ---------------------------------------------
bootstrap_stratify_val <- Sys.getenv("SCPROTXPRESS_BOOTSTRAP_STRATIFY", "cell.type_tissue")
apply_sqs_val          <- Sys.getenv("SCPROTXPRESS_APPLY_SQS",          "ref_only")
threshold_method_val   <- Sys.getenv("SCPROTXPRESS_THRESHOLD_METHOD",   "mixture")

ncores_val        <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

message(sprintf(
  "[CONFIG] bootstrap_stratify=%s  apply_sqs=%s (threshold=%.2f)  threshold_method=%s",
  bootstrap_stratify_val, apply_sqs_val, sqs_threshold_val, threshold_method_val))

# ----------------------------------------------------------------------------
# 2. PHASE 1 -- REFERENCE
# ----------------------------------------------------------------------------
phase_reference <- function() {
  message("[PHASE 1] Building CyTOF reference background")
  sce <- readRDS(infile_ref)
  unpacked <- unpack_sce_cytof(sce, assay_name = p_assay,
                               batch_default_col = "dataset_id")

  meta <- unpacked$metadata
  keep <- meta$cell_type != "Unclear/Other"
  meta <- meta[keep, , drop = FALSE]
  data_mat <- unpacked$data[, keep, drop = FALSE]

  run_scProteinXpress(
    data               = data_mat,
    metadata           = meta,
    dir                = outdir,
    file.name          = file_name_ref,
    modality           = "CyTOF",
    cell_id_col        = "cell_id",
    cell_type_col      = "cell_type",
    dataset_col        = "dataset_id",
    sample_col         = "sample_id",
    batch_col          = "batch",
    tissue_default     = "PBMC",
    discretize         = TRUE,
    sig.frames         = TRUE,
    num.boot.sample    = nb,
    num.cells.sample   = ns,
    ncores             = ncores_val,
    is_transformed     = unpacked$is_transformed,
    sqs_threshold      = sqs_threshold_val,
    apply_sqs          = apply_sqs_val,
    bootstrap_stratify = bootstrap_stratify_val,
    threshold_method   = threshold_method_val
  )
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# 3. PHASE 2 -- EXTERNAL QUERY
# ----------------------------------------------------------------------------
run_query <- function(infile, file.name, precBack_path) {
  message("[PHASE 2] Query ", file.name)
  sce <- readRDS(infile)
  unpacked <- unpack_sce_cytof(sce, assay_name = p_assay,
                               batch_default_col = "dataset_id")
  meta <- unpacked$metadata
  keep <- meta$cell_type != "Unclear/Other"
  meta <- meta[keep, , drop = FALSE]
  data_mat <- unpacked$data[, keep, drop = FALSE]

  out_q <- file.path(outdir, paste0("Query_", file.name))
  dir.create(out_q, recursive = TRUE, showWarnings = FALSE)

  run_scProteinXpress(
    data            = data_mat,
    metadata        = meta,
    dir             = out_q,
    file.name       = file.name,
    precBack        = precBack_path,
    modality        = "CyTOF",
    cell_id_col     = "cell_id",
    cell_type_col   = "cell_type",
    dataset_col     = "dataset_id",
    sample_col      = "sample_id",
    batch_col       = "batch",
    tissue_default  = "PBMC",
    ncores          = ncores_val,
    is_transformed  = unpacked$is_transformed,
    apply_sqs       = apply_sqs_val
  )
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# 4. ENTRY POINT
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  phase_reference()
  precBack <- file.path(outdir, paste0("Ref_", file_name_ref))
  if (file.exists(ds_query)) {
    run_query(ds_query, "external", precBack)
  }
}
