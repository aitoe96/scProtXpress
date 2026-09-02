## Tests for .compute_sig_table_for_population() -- the dispatch logic that
## decides whether a population's sig-table uses the original pooled-cell
## test or a donor-aware one, inside run_scProteinXpress()'s sig.frames
## step; see test-calibration.R for the underlying statistical validation
## of discrete.expr.levels_donor() itself.

if (requireNamespace("scProtXpress", quietly = TRUE) &&
    !exists(".compute_sig_table_for_population", mode = "function")) {
  .compute_sig_table_for_population <- getFromNamespace(".compute_sig_table_for_population", "scProtXpress")
}

.make_donor_clustered_mat <- function(seed, D = 10, n_per_donor = 120, n_genes = 100, rho = 0.10, prob = 0.3) {
  set.seed(seed)
  donor_ids <- rep(paste0("donor", seq_len(D)), each = n_per_donor)
  cell_ids  <- paste0("cell_", seq_along(donor_ids))
  kappa <- (1 - rho) / rho
  mat <- do.call(rbind, lapply(seq_len(n_genes), function(i) {
    donor_rate <- rbeta(D, prob * kappa, (1 - prob) * kappa)
    unlist(lapply(seq_len(D), function(d) rbinom(n_per_donor, 1, donor_rate[d])))
  }))
  rownames(mat) <- paste0("CD", seq_len(n_genes)); colnames(mat) <- cell_ids
  list(mat = mat, pop_meta = data.frame(cell.id = cell_ids, sample = donor_ids, stringsAsFactors = FALSE))
}

test_that("default sig_test_method='pooled' reproduces the original Create.sig.frame()+discrete.expr.levels() call exactly", {
  sim <- .make_donor_clustered_mat(seed = 1)
  via_dispatch <- .compute_sig_table_for_population(sim$mat, sim$pop_meta, 0, tempfile(), sig_test_method = "pooled")

  sig.res <- Create.sig.frame(sim$mat, min.coverage.fraction.sigframe = 0, audit_csv_path = tempfile())
  via_direct <- discrete.expr.levels(sig.res[[1]], sig.res[[2]], prob_low = sig.res[[3]], prob_high = sig.res[[4]])

  expect_identical(via_dispatch$expr.level, via_direct$expr.level)
})

test_that("sig_test_method='design_effect' meaningfully reduces the false-positive rate under realistic donor clustering", {
  sim <- .make_donor_clustered_mat(seed = 2, rho = 0.10)
  tbl_pooled <- .compute_sig_table_for_population(sim$mat, sim$pop_meta, 0, tempfile(), sig_test_method = "pooled")
  tbl_deff   <- .compute_sig_table_for_population(sim$mat, sim$pop_meta, 0, tempfile(), sig_test_method = "design_effect")

  fpr_pooled <- mean(tbl_pooled$expr.level %in% c("medium", "high"))
  fpr_deff   <- mean(tbl_deff$expr.level %in% c("medium", "high"))

  ## Every gene here is a true null (no real signal) -- pooled should be
  ## well above nominal 0.05 under this clustering, design_effect should
  ## be close to it.
  expect_gt(fpr_pooled, 0.10)
  expect_lt(fpr_deff, 0.12)
})

test_that("a donor-aware method falls back to pooled with a warning when no donor metadata is available", {
  sim <- .make_donor_clustered_mat(seed = 3)
  tbl_pooled <- .compute_sig_table_for_population(sim$mat, sim$pop_meta, 0, tempfile(), sig_test_method = "pooled")

  expect_warning(
    tbl_fallback <- .compute_sig_table_for_population(sim$mat, NULL, 0, tempfile(), sig_test_method = "donor_t"),
    "no per-cell"
  )
  expect_identical(tbl_fallback$expr.level, tbl_pooled$expr.level)
})

test_that("a donor-aware method falls back to pooled with a warning when metadata has no 'sample' column", {
  sim <- .make_donor_clustered_mat(seed = 4)
  bad_meta <- sim$pop_meta
  colnames(bad_meta)[colnames(bad_meta) == "sample"] <- "donor_but_wrong_name"

  expect_warning(
    .compute_sig_table_for_population(sim$mat, bad_meta, 0, tempfile(), sig_test_method = "donor_wilcox"),
    "no per-cell"
  )
})

test_that("a donor-aware method falls back to pooled with a warning on unmatched cell IDs", {
  sim <- .make_donor_clustered_mat(seed = 5)
  scrambled_meta <- sim$pop_meta
  scrambled_meta$cell.id[1:5] <- paste0("not_a_real_cell_", 1:5)

  expect_warning(
    .compute_sig_table_for_population(sim$mat, scrambled_meta, 0, tempfile(), sig_test_method = "design_effect"),
    "could not align"
  )
})
