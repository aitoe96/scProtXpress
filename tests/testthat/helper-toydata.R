## Synthetic toy data factory.
##
## Builds a deterministic proteins x cells matrix big enough for the
## reference threshold (>= 50 cells per cell.type) and a metadata frame
## with the canonical CyTOF column names (cell_id, cell_type, dataset_id,
## sample_id, batch). Pass `style = "citeseq"` for the codebase-style
## names (cell.id, cell.type_harmonised, dataset, sample, batch).

.make_toy_counts <- function(n_markers = 20L, n_cells_per_ct = 60L,
                             n_cell_types = 3L, n_samples = 3L,
                             seed = 42L, transformed = FALSE) {
  set.seed(seed)
  n_cells <- n_cells_per_ct * n_cell_types
  cell_ids   <- sprintf("c%05d", seq_len(n_cells))
  cell_types <- rep(paste0("ct", seq_len(n_cell_types)),
                    each = n_cells_per_ct)
  samples    <- rep(paste0("s", seq_len(n_samples)),
                    length.out = n_cells)

  if (transformed) {
    counts <- matrix(stats::rnorm(n_markers * n_cells, mean = 3, sd = 1.5),
                     nrow = n_markers, ncol = n_cells)
  } else {
    counts <- matrix(stats::rpois(n_markers * n_cells, lambda = 5),
                     nrow = n_markers, ncol = n_cells)
  }
  rownames(counts) <- sprintf("CD%d", seq_len(n_markers))
  colnames(counts) <- cell_ids
  list(
    counts     = counts,
    cell_ids   = cell_ids,
    cell_types = cell_types,
    samples    = samples
  )
}

make_toy_cytof <- function(...) {
  raw <- .make_toy_counts(transformed = TRUE, ...)
  meta <- data.frame(
    cell_id    = raw$cell_ids,
    cell_type  = raw$cell_types,
    dataset_id = "DSx",
    sample_id  = raw$samples,
    batch      = raw$samples,
    tissue     = "PBMC",
    stringsAsFactors = FALSE
  )
  list(data = raw$counts, metadata = meta, is_transformed = TRUE)
}

make_toy_citeseq <- function(...) {
  raw <- .make_toy_counts(transformed = FALSE, ...)
  meta <- data.frame(
    cell.id              = raw$cell_ids,
    cell.type_harmonised = raw$cell_types,
    dataset              = "DSy",
    sample               = raw$samples,
    batch                = raw$samples,
    tissue               = "PBMC",
    stringsAsFactors = FALSE
  )
  list(
    counts         = methods::as(raw$counts, "dgCMatrix"),
    cell_x_feature = meta
  )
}
