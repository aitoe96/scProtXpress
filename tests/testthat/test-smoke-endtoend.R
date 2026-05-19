## End-to-end smoke test against synthetic data. ADTnorm + flowCore must
## be installed (or the test is skipped). Runs with a tiny bootstrap
## count so the whole thing fits in seconds.

test_that("CyTOF smoke run produces the expected Ref outputs", {
  skip_if_not_installed("ADTnorm")
  skip_if_not_installed("flowCore")
  skip_on_cran()

  toy   <- make_toy_cytof(n_cells_per_ct = 60L, n_markers = 6L)
  outdir <- tempfile("scprot-smoke-")
  dir.create(outdir, recursive = TRUE)

  expect_invisible(suppressMessages(
    run_scProteinXpress(
      data              = toy$data,
      metadata          = toy$metadata,
      dir               = outdir,
      file.name         = "smoke",
      modality          = "CyTOF",
      cell_id_col       = "cell_id",
      cell_type_col     = "cell_type",
      dataset_col       = "dataset_id",
      sample_col        = "sample_id",
      batch_col         = "batch",
      tissue_default    = "PBMC",
      is_transformed    = TRUE,
      num.boot.sample   = 5L,
      num.cells.sample  = 20L,
      min.cells.feature = 1L,
      min.coverage.fraction = 0.0,
      ncores            = 1L,
      do_subsample      = FALSE
    )
  ))

  ref_dir <- file.path(outdir, "Ref_smoke")
  expect_true(file.exists(file.path(ref_dir, "smoke_NormData.rds")))
  expect_true(file.exists(file.path(ref_dir, "smoke_bg_prot_pars.rds")))
  expect_true(file.exists(file.path(ref_dir, "smoke_th_dist.rds")))
})
