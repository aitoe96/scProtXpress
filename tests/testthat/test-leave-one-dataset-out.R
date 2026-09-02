## Tests for leave_one_dataset_out_stability().

.make_loo_test_data <- function(seed = 1) {
  set.seed(seed)
  datasets <- paste0("ds", 1:5)
  n_per_ds <- 150
  dataset_vec <- rep(datasets, each = n_per_ds)
  n_total <- length(dataset_vec)
  cell_ids <- paste0("cell_", seq_len(n_total))

  ## Stable marker: identical bimodal distribution in every dataset.
  CD_stable <- unlist(lapply(datasets, function(d) {
    c(rnorm(n_per_ds * 0.6, 1, 0.5), rnorm(n_per_ds * 0.4, 5, 0.6))[seq_len(n_per_ds)]
  }))
  ## Unstable marker: one dataset (ds5) has a real, different distribution.
  CD_unstable <- unlist(lapply(datasets, function(d) {
    if (d == "ds5") c(rnorm(n_per_ds * 0.6, 1, 0.5), rnorm(n_per_ds * 0.4, 9, 0.6))[seq_len(n_per_ds)]
    else            c(rnorm(n_per_ds * 0.6, 1, 0.5), rnorm(n_per_ds * 0.4, 5, 0.6))[seq_len(n_per_ds)]
  }))
  ## Rare marker: only measured in 2 datasets -- must be excluded at the
  ## default min_datasets_per_marker = 3.
  CD_rare <- rep(NA_real_, n_total)
  in_rare <- dataset_vec %in% c("ds1", "ds2")
  CD_rare[in_rare] <- rnorm(sum(in_rare), 2, 0.5)

  data_norm <- rbind(CD_stable = CD_stable, CD_unstable = CD_unstable, CD_rare = CD_rare)
  colnames(data_norm) <- cell_ids
  metadata <- data.frame(cell.id = cell_ids, dataset = dataset_vec, cell.type = "x",
                         stringsAsFactors = FALSE)
  list(data_norm = data_norm, metadata = metadata)
}

test_that("markers below min_datasets_per_marker are excluded from the analysis", {
  sim <- .make_loo_test_data(1)
  res <- leave_one_dataset_out_stability(sim$data_norm, sim$metadata,
                                         min_datasets_per_marker = 3,
                                         sample_n = 100, num.boot.sample = 15, fixseed = 1)
  expect_setequal(res$stability_summary$marker, c("CD_stable", "CD_unstable"))
  expect_false("CD_rare" %in% res$stability_summary$marker)
})

test_that("a marker with a real dataset-specific shift is flagged as less stable than one without", {
  sim <- .make_loo_test_data(2)
  res <- leave_one_dataset_out_stability(sim$data_norm, sim$metadata,
                                         min_datasets_per_marker = 3,
                                         sample_n = 100, num.boot.sample = 20, fixseed = 2)
  s <- res$stability_summary
  stable_shift   <- s$max_abs_shift_upper[s$marker == "CD_stable"]
  unstable_shift <- s$max_abs_shift_upper[s$marker == "CD_unstable"]
  expect_gt(unstable_shift, stable_shift)
})

test_that("errors are clear and immediate for malformed inputs", {
  sim <- .make_loo_test_data(3)
  expect_error(
    leave_one_dataset_out_stability(sim$data_norm, sim$metadata[, "cell.id", drop = FALSE]),
    "dataset"
  )
  expect_error(
    leave_one_dataset_out_stability(sim$data_norm, sim$metadata[1:10, ]),
    "ncol.*nrow|positionally"
  )
  only_two_ds <- sim$metadata
  only_two_ds$dataset[only_two_ds$dataset == "ds5"] <- "ds4"
  only_two_ds$dataset[!(only_two_ds$dataset %in% c("ds1", "ds2"))] <- "ds1"
  expect_error(
    leave_one_dataset_out_stability(sim$data_norm, only_two_ds),
    "at least 3"
  )
})
