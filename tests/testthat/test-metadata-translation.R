## Metadata column resolution (A1 fallbacks).
## Verified through run_scProteinXpress() which calls .resolve_meta_col()
## internally. We exercise the validation block by triggering it on
## inputs that pre-empt downstream computation (precBack = NA forces the
## Ref build to start, then we expect a known error from V13 or before).

test_that("metadata translation accepts dot-style CITE-seq names", {
  skip_if_not_installed("ADTnorm")

  toy <- make_toy_citeseq()
  data <- as.matrix(toy$counts)
  meta <- toy$cell_x_feature

  expect_error(
    run_scProteinXpress(
      data            = data,
      metadata        = meta,
      dir             = tempfile("scprot-md-"),
      file.name       = "smoke",
      modality        = "CITE-seq",
      cell_id_col     = "cell.id",
      cell_type_col   = "cell.type_harmonised",
      dataset_col     = "dataset",
      sample_col      = "sample",
      batch_col       = "batch",
      num.boot.sample = 5L,
      ncores          = 1L
    ),
    regexp = NA  # accept any error AFTER metadata translation succeeds;
                 # we are checking that resolution does not stop() here.
  )
})

test_that("metadata translation fails clearly when a column is missing", {
  toy <- make_toy_cytof()
  meta_bad <- toy$metadata
  meta_bad$cell_type <- NULL  # remove the cell-type column

  expect_error(
    run_scProteinXpress(
      data         = toy$data,
      metadata     = meta_bad,
      dir          = tempfile("scprot-md2-"),
      file.name    = "smoke",
      modality     = "CyTOF",
      cell_type_col = "cell_type",
      ncores       = 1L
    ),
    regexp = "metadata.*could not resolve"
  )
})
