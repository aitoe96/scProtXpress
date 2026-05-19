## D2 follow-up (Task #14): verify qc_outlier_cells switch behaviour.

test_that("cleanData_ADT skips the outlier filter when qc_outlier_cells = 'no'", {
  toy <- make_toy_cytof()
  out_no <- cleanData_ADT(
    data_counts           = toy$data,
    nmads                 = 3,
    min.cells.feature     = 1L,
    min.coverage.fraction = 0.0,
    qc_outlier_cells      = "no"
  )
  expect_equal(ncol(out_no$Counts), ncol(toy$data))
})

test_that("cleanData_ADT applies scuttle::isOutlier when qc_outlier_cells = 'yes'", {
  skip_if_not_installed("scuttle")
  set.seed(123)
  toy <- make_toy_cytof()

  # inject a cluster of very low-signal cells (~10) that scuttle should flag.
  bad_idx <- sample(seq_len(ncol(toy$data)), 12)
  toy$data[, bad_idx] <- 0
  toy$data[1, bad_idx] <- 1  # keep them non-empty (escape the colSums==0 prune)

  out_yes <- cleanData_ADT(
    data_counts           = toy$data,
    nmads                 = 1,
    min.cells.feature     = 1L,
    min.coverage.fraction = 0.0,
    qc_outlier_cells      = "yes"
  )
  expect_lt(ncol(out_yes$Counts), ncol(toy$data))
})
