#' CITE-seq driver template (Hito 0 Q1 resolution)
#'
#' Mirror of `cytof_driver.R` for the CITE-seq path. Reference is built
#' from a list `{$counts, $cell_x_feature}` (the conventional CITE-seq
#' query object). External query also uses the same shape against a
#' pre-built Ref.
#'
#' Phases:
#'   1. `phase_reference()`: build CITE-seq reference background.
#'   2. `run_query()`: query an external CITE-seq dataset against the Ref.

suppressPackageStartupMessages({
  library(Matrix)
  library(scProtXpress)
})

# ----------------------------------------------------------------------------
# 1. CONFIG (edit per run)
# ----------------------------------------------------------------------------
data_root  <- Sys.getenv("SCPROTXPRESS_DATA_ROOT", "/scratch/USER/data_citeseq")
infile_ref <- file.path(data_root, "list_reference.rds")
ds_query   <- file.path(data_root, "list_query.rds")
out_root   <- Sys.getenv("SCPROTXPRESS_OUT_ROOT", "/scratch/USER/results_citeseq")
run_name   <- Sys.getenv("SCPROTXPRESS_RUN_NAME", "citeseq_run")
outdir     <- file.path(out_root, run_name)
file_name_ref <- "CITEseq_ADTNorm"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

nb                <- 1000      # bootstrap iterations
ns                <- 100       # CITE-seq default (smaller per-pop)
sqs_threshold_val <- 1.0
apply_sqs_val     <- "none"
ncores_val        <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

# ----------------------------------------------------------------------------
# 2. PHASE 1 -- REFERENCE
# ----------------------------------------------------------------------------
phase_reference <- function() {
  message("[PHASE 1] Building CITE-seq reference background")
  ref_obj  <- readRDS(infile_ref)
  unpacked <- unpack_list_citeseq(ref_obj)

  meta <- unpacked$metadata
  # CITE-seq queries may use `cell.type_harmonised`; cell_type_col fallback
  # in run_scProteinXpress handles this automatically.

  run_scProteinXpress(
    data             = unpacked$data,
    metadata         = meta,
    dir              = outdir,
    file.name        = file_name_ref,
    modality         = "CITE-seq",
    cell_id_col      = "cell.id",
    cell_type_col    = "cell.type_harmonised",
    dataset_col      = "dataset",
    sample_col       = "sample",
    batch_col        = "batch",
    tissue_default   = "PBMC",
    discretize       = TRUE,
    sig.frames       = TRUE,
    num.boot.sample  = nb,
    num.cells.sample = ns,
    ncores           = ncores_val,
    is_transformed   = unpacked$is_transformed,
    qc_outlier_cells = "auto",      # CITE-seq => scuttle filter active.
    sqs_threshold    = sqs_threshold_val,
    apply_sqs        = apply_sqs_val
  )
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# 3. PHASE 2 -- EXTERNAL QUERY
# ----------------------------------------------------------------------------
run_query <- function(infile, file.name, precBack_path) {
  message("[PHASE 2] CITE-seq query ", file.name)
  qobj     <- readRDS(infile)
  unpacked <- unpack_list_citeseq(qobj)

  out_q <- file.path(outdir, paste0("Query_", file.name))
  dir.create(out_q, recursive = TRUE, showWarnings = FALSE)

  run_scProteinXpress(
    data             = unpacked$data,
    metadata         = unpacked$metadata,
    dir              = out_q,
    file.name        = file.name,
    precBack         = precBack_path,
    modality         = "CITE-seq",
    cell_id_col      = "cell.id",
    cell_type_col    = "cell.type_harmonised",
    dataset_col      = "dataset",
    sample_col       = "sample",
    batch_col        = "batch",
    tissue_default   = "PBMC",
    ncores           = ncores_val,
    is_transformed   = unpacked$is_transformed,
    qc_outlier_cells = "auto",
    apply_sqs        = apply_sqs_val
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
