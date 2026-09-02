## Tests for downsample_reference_stability().

test_that("n_cells_used exactly matches the requested fractions", {
  set.seed(10)
  n <- 1000
  data_norm <- matrix(rnorm(n), nrow = 1, dimnames = list("CD1", NULL))
  metadata <- data.frame(cell.id = paste0("c", 1:n), dataset = sample(paste0("ds", 1:4), n, replace = TRUE),
                        cell.type = "x", stringsAsFactors = FALSE)

  res <- downsample_reference_stability(data_norm, metadata,
                                        fractions = c(0.2, 0.5, 1.0), subsample_by = "cells",
                                        n_replicates = 2, sample_n = 100, num.boot.sample = 5, fixseed = 1)

  expect_setequal(unique(res$n_cells_used[res$fraction == 0.2]), round(n * 0.2))
  expect_setequal(unique(res$n_cells_used[res$fraction == 0.5]), round(n * 0.5))
  expect_setequal(unique(res$n_cells_used[res$fraction == 1.0]), n)
})

test_that("fraction = 1.0 always has zero shift by construction, and only one replicate", {
  set.seed(11)
  n <- 800
  data_norm <- matrix(rnorm(n), nrow = 1, dimnames = list("CD1", NULL))
  metadata <- data.frame(cell.id = paste0("c", 1:n), dataset = sample(paste0("ds", 1:3), n, replace = TRUE),
                        cell.type = "x", stringsAsFactors = FALSE)

  res <- downsample_reference_stability(data_norm, metadata,
                                        fractions = c(0.3, 1.0), n_replicates = 3,
                                        sample_n = 100, num.boot.sample = 5, fixseed = 2)

  full <- res[res$fraction == 1.0, ]
  expect_equal(nrow(full), 1L)
  expect_equal(full$abs_shift_lower, 0)
  expect_equal(full$abs_shift_upper, 0)
})

test_that("a small fraction shows a larger average shift from the full reference than a large one", {
  set.seed(12)
  n <- 2500
  mk <- c(rnorm(n * 0.6, 1, 0.6), rnorm(n * 0.4, 5, 0.7))
  data_norm <- matrix(mk, nrow = 1, dimnames = list("CD1", NULL))
  metadata <- data.frame(cell.id = paste0("c", 1:n), dataset = sample(paste0("ds", 1:5), n, replace = TRUE),
                        cell.type = "x", stringsAsFactors = FALSE)

  res <- downsample_reference_stability(data_norm, metadata,
                                        fractions = c(0.05, 0.75, 1.0), subsample_by = "cells",
                                        n_replicates = 6, sample_n = 100, num.boot.sample = 20, fixseed = 3)

  mean_shift_small <- mean(res$abs_shift_lower[res$fraction == 0.05] + res$abs_shift_upper[res$fraction == 0.05])
  mean_shift_large  <- mean(res$abs_shift_lower[res$fraction == 0.75] + res$abs_shift_upper[res$fraction == 0.75])
  expect_gt(mean_shift_small, mean_shift_large)
})

test_that("replicates within the same fraction are genuinely independent, not duplicates of each other", {
  ## This is the specific gap that let a real bug ship undetected: earlier
  ## versions of this function called GetBackground() (which internally
  ## does set.seed(fixseed) on every call) repeatedly within the replicate
  ## loop, which silently reset the RNG to the same state after every call
  ## and collapsed replicates 2+ into exact duplicates of each other. The
  ## other tests in this file check aggregate trends across replicates,
  ## which can look fine even when the "replicates" aren't actually
  ## independent -- this test checks independence directly.
  set.seed(20)
  n <- 2000
  data_norm <- matrix(rnorm(n, 2, 1), nrow = 1, dimnames = list("CD1", NULL))
  metadata <- data.frame(cell.id = paste0("c", 1:n), dataset = sample(paste0("ds", 1:6), n, replace = TRUE),
                        cell.type = "x", stringsAsFactors = FALSE)

  res <- downsample_reference_stability(data_norm, metadata,
                                        fractions = c(0.05, 1.0), subsample_by = "cells",
                                        n_replicates = 8, sample_n = 100, num.boot.sample = 30, fixseed = 1)
  sub <- res[res$fraction == 0.05, ]
  expect_equal(nrow(sub), 8L)
  ## No two replicates' Upper threshold should be exactly identical at full
  ## numeric precision (genuine independent random subsamples essentially
  ## never coincide exactly; duplicates indicate the RNG-collapse bug).
  expect_equal(sum(duplicated(sub$Upper)), 0L)
  expect_equal(sum(duplicated(sub$Lower)), 0L)
})

test_that("downsample_reference_stability is still fully reproducible given the same fixseed", {
  set.seed(21)
  n <- 1000
  data_norm <- matrix(rnorm(n, 2, 1), nrow = 1, dimnames = list("CD1", NULL))
  metadata <- data.frame(cell.id = paste0("c", 1:n), dataset = sample(paste0("ds", 1:4), n, replace = TRUE),
                        cell.type = "x", stringsAsFactors = FALSE)
  args <- list(data_norm = data_norm, metadata = metadata, fractions = c(0.2, 1.0),
              subsample_by = "cells", n_replicates = 3, sample_n = 100, num.boot.sample = 10, fixseed = 7)
  res1 <- do.call(downsample_reference_stability, args)
  res2 <- do.call(downsample_reference_stability, args)
  expect_identical(res1$Upper, res2$Upper)
  expect_identical(res1$Lower, res2$Lower)
})

test_that("subsample_by = 'datasets' requires a dataset column and behaves correctly when present", {
  set.seed(13)
  n <- 600
  data_norm <- matrix(rnorm(n), nrow = 1, dimnames = list("CD1", NULL))
  metadata_no_ds <- data.frame(cell.id = paste0("c", 1:n), cell.type = "x", stringsAsFactors = FALSE)
  expect_error(
    downsample_reference_stability(data_norm, metadata_no_ds, subsample_by = "datasets"),
    "dataset"
  )

  metadata <- data.frame(cell.id = paste0("c", 1:n), dataset = sample(paste0("ds", 1:6), n, replace = TRUE),
                        cell.type = "x", stringsAsFactors = FALSE)
  res <- downsample_reference_stability(data_norm, metadata, fractions = c(0.5, 1.0),
                                        subsample_by = "datasets", n_replicates = 2,
                                        sample_n = 50, num.boot.sample = 5, fixseed = 4)
  expect_true(all(res$n_cells_used[res$fraction == 1.0] == n))
  ## With 6 datasets and fraction 0.5, exactly 3 datasets should be kept each
  ## replicate -- can't check dataset membership directly from the returned
  ## summary, but n_cells_used should be less than n and reproducible given
  ## the fixed seed.
  expect_true(all(res$n_cells_used[res$fraction == 0.5] < n))
})
