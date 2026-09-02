## bootstrap_stratify dispatcher tests.
##
## Exercises GetBackground() in all three modes plus the argument
## validation surfaced by run_scProteinXpress(). The full run_scProteinXpress
## end-to-end smoke test is covered by test-smoke-endtoend.R and skips
## when ADTnorm is unavailable; these tests stay local because they call
## the bootstrap helpers directly.

testthat::test_that("bootstrap_stratify = 'none' produces a per-marker thresholds object", {
  skip_if_not_installed("foreach")

  toy <- .make_toy_counts(n_markers = 5L, n_cells_per_ct = 60L,
                          n_cell_types = 2L, n_samples = 2L,
                          transformed = TRUE, seed = 11L)
  meta <- data.frame(
    cell.id   = toy$cell_ids,
    cell.type = toy$cell_types,
    tissue    = "PBMC",
    stringsAsFactors = FALSE
  )

  outdir <- tempfile("scProtXpress_test_none_")
  thrs <- GetBackground(
    data_norm          = toy$counts,
    dir_name           = outdir,
    file.id            = "toy_none",
    metadata           = meta,
    sample_n           = 50L,
    parallelize        = FALSE,
    ncores             = 1L,
    num.boot.sample    = 4L,
    fixseed            = 1L,
    bootstrap_stratify = "none"
  )

  testthat::expect_type(thrs, "list")
  testthat::expect_true(all(rownames(toy$counts) %in% names(thrs)))
  testthat::expect_true(all(vapply(thrs, function(m) all(colnames(m) == c("Lower","Upper")), logical(1))))
  testthat::expect_true(file.exists(file.path(outdir, "toy_none_th_dist.rds")))
})

testthat::test_that("bootstrap_stratify = 'cell.type' still runs (stratified path)", {
  skip_if_not_installed("foreach")

  toy <- .make_toy_counts(n_markers = 5L, n_cells_per_ct = 60L,
                          n_cell_types = 2L, n_samples = 2L,
                          transformed = TRUE, seed = 12L)
  meta <- data.frame(
    cell.id   = toy$cell_ids,
    cell.type = toy$cell_types,
    stringsAsFactors = FALSE
  )

  outdir <- tempfile("scProtXpress_test_ct_")
  thrs <- GetBackground(
    data_norm          = toy$counts,
    dir_name           = outdir,
    file.id            = "toy_ct",
    metadata           = meta,
    sample_n           = 50L,
    parallelize        = FALSE,
    ncores             = 1L,
    num.boot.sample    = 4L,
    fixseed            = 2L,
    bootstrap_stratify = "cell.type"
  )

  testthat::expect_type(thrs, "list")
  testthat::expect_true(length(thrs) >= 1L)
})

testthat::test_that("bootstrap_stratify = 'cell.type_tissue' still runs (two-level stratified path)", {
  skip_if_not_installed("foreach")

  ## Sized so that each (cell.type, tissue) stratum has 30 cells and the
  ## stratified path's `length(a) >= n` filter (n = sample_n / n_tissues)
  ## clears with sample_n = 50.
  toy <- .make_toy_counts(n_markers = 5L, n_cells_per_ct = 60L,
                          n_cell_types = 2L, n_samples = 2L,
                          transformed = TRUE, seed = 13L)
  meta <- data.frame(
    cell.id   = toy$cell_ids,
    cell.type = toy$cell_types,
    tissue    = rep(c("PBMC", "BM"), length.out = length(toy$cell_ids)),
    stringsAsFactors = FALSE
  )

  outdir <- tempfile("scProtXpress_test_ct_tissue_")
  thrs <- GetBackground(
    data_norm          = toy$counts,
    dir_name           = outdir,
    file.id            = "toy_ct_tissue",
    metadata           = meta,
    sample_n           = 50L,
    parallelize        = FALSE,
    ncores             = 1L,
    num.boot.sample    = 4L,
    fixseed            = 3L,
    bootstrap_stratify = "cell.type_tissue"
  )

  testthat::expect_type(thrs, "list")
  testthat::expect_true(length(thrs) >= 1L)
})

testthat::test_that("bootstrap_stratify = 'cell.type_tissue' without tissue column warns and falls back", {
  skip_if_not_installed("foreach")

  toy <- .make_toy_counts(n_markers = 4L, n_cells_per_ct = 60L,
                          n_cell_types = 2L, n_samples = 2L,
                          transformed = TRUE, seed = 14L)
  meta <- data.frame(
    cell.id   = toy$cell_ids,
    cell.type = toy$cell_types,
    stringsAsFactors = FALSE
  )

  outdir <- tempfile("scProtXpress_test_ct_tissue_fallback_")
  testthat::expect_warning(
    thrs <- GetBackground(
      data_norm          = toy$counts,
      dir_name           = outdir,
      file.id            = "toy_fb",
      metadata           = meta,
      sample_n           = 50L,
      parallelize        = FALSE,
      ncores             = 1L,
      num.boot.sample    = 4L,
      fixseed            = 4L,
      bootstrap_stratify = "cell.type_tissue"
    ),
    regexp = "falling back to cell.type"
  )
  testthat::expect_type(thrs, "list")
})

testthat::test_that("invalid bootstrap_stratify is rejected by match.arg", {
  toy <- .make_toy_counts(n_markers = 3L, n_cells_per_ct = 60L,
                          n_cell_types = 2L, n_samples = 2L,
                          transformed = TRUE, seed = 15L)
  meta <- data.frame(
    cell.id   = toy$cell_ids,
    cell.type = toy$cell_types,
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    GetBackground(
      data_norm          = toy$counts,
      dir_name           = tempfile("scProtXpress_test_bad_"),
      file.id            = "bad",
      metadata           = meta,
      sample_n           = 50L,
      parallelize        = FALSE,
      ncores             = 1L,
      num.boot.sample    = 2L,
      fixseed            = 5L,
      bootstrap_stratify = "everything"
    )
  )
})
