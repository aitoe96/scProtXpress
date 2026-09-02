#' From an already-ADTnorm-normalized list(counts, metadata) to the
#' background/threshold calculation
#'
#' Use this when you've already run ADTnorm normalization yourself (or your
#' `counts` matrix is otherwise already on the aligned/arcsinh scale
#' scProtXpress expects) and just want the bootstrap threshold step --
#' Fig. 1a's "background building" -- without re-running
#' `run_scProteinXpress()`'s own normalization (`ADTprot_norm()`), which
#' would redo peak/valley alignment you don't need.
#'
#' If your data is NOT yet normalized (raw ADT/CyTOF counts), use
#' `run_scProteinXpress()` instead (see inst/drivers/citeseq_driver.R /
#' cytof_driver.R) -- it calls ADTnorm normalization AND this same
#' `GetBackground()` step for you, plus QC and SQS filtering.
#'
#' Tested end-to-end against synthetic data shaped like this (600 cells,
#' 15 markers, 3 cell types, 2 tissues, 6 donors) -- this script runs with
#' no ADTnorm/flowCore dependency at all, since normalization is exactly
#' the part being skipped.

suppressPackageStartupMessages({
  library(Matrix)
  library(scProtXpress)
})

# ----------------------------------------------------------------------------
# 1. LOAD THE LIST
# ----------------------------------------------------------------------------
# Expected shape: list(counts = <markers x cells matrix, ADTnorm-normalized>,
#                       cell_x_feature = <one row per cell, cell-id-keyed>)
# This is the same shape unpack_list_citeseq() / citeseq_driver.R expect for
# RAW data; here `counts` is just already normalized, so we skip the
# unpack_list_citeseq() + ADTprot_norm() steps entirely.
adtnorm_data <- readRDS("path/to/your_adtnorm_normalized_list.rds")

counts         <- adtnorm_data$counts            # markers x cells
cell_x_feature <- adtnorm_data$cell_x_feature     # metadata, one row per cell

# ----------------------------------------------------------------------------
# 2. RENAME METADATA COLUMNS TO scProtXpress'S CANONICAL NAMES
# ----------------------------------------------------------------------------
# GetBackground() looks specifically for `cell.type` and `tissue` columns
# when stratifying the bootstrap (see below) -- edit the left-hand side of
# each line to match whatever your own metadata actually calls these.
meta <- cell_x_feature
colnames(meta)[colnames(meta) == "cell_type"] <- "cell.type"   # <- EDIT to your column name
colnames(meta)[colnames(meta) == "tissue_id"] <- "tissue"      # <- EDIT to your column name
colnames(meta)[colnames(meta) == "sample_id"] <- "sample"      # <- EDIT to your column name

# If every cell is from one tissue (e.g. PBMC-only data), there may be no
# tissue column at all -- add a constant one so cell.type_tissue
# stratification still works (it becomes equivalent to cell.type alone,
# but keeps the same code path / same downstream file format):
if (!"tissue" %in% colnames(meta)) meta$tissue <- "PBMC"        # <- EDIT if not PBMC
if (!"batch"  %in% colnames(meta)) meta$batch  <- meta$sample

# ----------------------------------------------------------------------------
# 3. ALIGN METADATA ROWS TO COUNTS COLUMNS
# ----------------------------------------------------------------------------
# GetBackground()/RunBootstrap() use metadata POSITIONALLY (row i of
# metadata must describe column i of counts) -- they do not re-match by
# cell.id internally, so get this right here or every downstream grouping
# is silently wrong. Edit "cell.id" below if your metadata's ID column has
# a different name.
stopifnot(all(colnames(counts) %in% meta$cell.id))
meta <- meta[match(colnames(counts), meta$cell.id), , drop = FALSE]
stopifnot(identical(meta$cell.id, colnames(counts)))   # hard stop if alignment failed

# ----------------------------------------------------------------------------
# 4. BACKGROUND CALCULATION
# ----------------------------------------------------------------------------
out_dir <- "path/to/your_output_dir"   # <- EDIT
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

th_dist <- GetBackground(
  data_norm = counts,
  dir_name  = out_dir,
  file.id   = "my_reference",          # <- EDIT: used as the output filename prefix
  metadata  = meta,

  # Cells resampled per bootstrap iteration. 100-200 is typical for
  # CITE-seq-sized populations; CyTOF drivers in this repo default to 1000.
  sample_n = 200L,

  # Bootstrap iterations. 1000 for a real reference; kept low here as a
  # quick-start default -- raise this for anything you intend to actually use.
  num.boot.sample = 200L,

  parallelize = TRUE,
  ncores      = 4L,       # <- EDIT to your machine/cluster
  fixseed     = 1234L,

  # Stratify the bootstrap by cell type AND tissue, not a flat pool across
  # everything. This is the setting the accompanying manuscript's Methods
  # describes; "none" (the package default) pools every cell for a marker
  # regardless of cell type/tissue composition, which is fine for a quick
  # look but means the resulting thresholds are more sensitive to whatever
  # happens to be over-represented in your input than you probably want
  # for a reference you'll reuse. Use "cell.type" if you don't have
  # meaningful tissue labels, or "none" to reproduce the older/simpler
  # behaviour.
  bootstrap_stratify = "cell.type_tissue",

  # Bayes-decision-boundary threshold from a fitted 1-3 component mixture
  # per marker, instead of the original single-curve rectangle heuristic.
  # Falls back to "rectangle" automatically per-marker if mclust isn't
  # installed or the fit is unstable/unimodal. Use "rectangle" to
  # reproduce the original behaviour exactly.
  threshold_method = "mixture"
)

# ----------------------------------------------------------------------------
# 5. WHAT YOU GET
# ----------------------------------------------------------------------------
# th_dist is a named list, one entry per marker, each an n_boot x 2
# (Lower, Upper) matrix of bootstrap threshold draws -- this IS "the
# background": the reference distributions in Fig. 1a, steps 3-5.
cat("Markers with thresholds:", length(th_dist), "\n")
print(head(th_dist[[1]]))

# Also written to disk under out_dir:
#   my_reference_th_dist.rds      -- th_dist itself (what Discretise_data() reads)
#   my_reference_sf.txt           -- per-feature scale factors (all 1 for this path)
#   my_reference_provenance.json  -- records file.id, package version, git commit
#                                     (if available), and every parameter above --
#                                     read this back later with
#                                     jsonlite::fromJSON() to confirm what a given
#                                     reference was actually built with, rather
#                                     than trusting a script's comments.

# ----------------------------------------------------------------------------
# NEXT STEPS (not run here)
# ----------------------------------------------------------------------------
# To classify query cells against this background:
#   Discretise_data(query_ScaledData, Ref_data = out_dir, file.id = "my_query")
# `query_ScaledData` needs to be on the SAME normalized scale as `counts`
# above -- if your query data went through the same external ADTnorm run,
# it already is; if not, align it to this reference's landmarks first
# (see ADTprot_norm(data.type = "query", bg.prot.pars = ...) /
# PrepQueryData_ADT(), which is what run_scProteinXpress() calls
# internally for query datasets).
