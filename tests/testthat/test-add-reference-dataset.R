## Unit tests for the bookkeeping helpers behind add_reference_dataset().
## These do NOT require ADTnorm/flowCore -- they test the merge/audit logic
## against mocked threshold-list objects, which is the part of
## add_reference_dataset() that provides the actual "incremental" guarantee
## (every marker not touched by the new dataset is left byte-for-byte
## unchanged). The ADTnorm-dependent alignment calls inside
## add_reference_dataset() itself are NOT covered here; see the file header
## of R/reference_update.R for what has and hasn't been verified.

## Internal-function access shim: see the identical block in
## test-calibration.R for why this is needed (library(scProtXpress) alone
## does not expose non-exported functions; devtools::load_all() does).
if (requireNamespace("scProtXpress", quietly = TRUE) &&
    !exists(".merge_threshold_lists", mode = "function")) {
  .merge_threshold_lists <- getFromNamespace(".merge_threshold_lists", "scProtXpress")
  .build_update_audit    <- getFromNamespace(".build_update_audit", "scProtXpress")
}

.mk_th <- function(lower_mean, upper_mean, seed) {
  set.seed(seed)
  cbind(Lower = rnorm(3, lower_mean, 0.01), Upper = rnorm(3, upper_mean, 0.01))
}

test_that("untouched markers are byte-identical after a merge", {
  old_th <- list(
    CD3  = .mk_th(1.0, 3.0, 1),
    CD4  = .mk_th(1.2, 3.2, 2),
    CD19 = .mk_th(1.1, 3.1, 3)
  )
  th_touched <- list(CD4 = .mk_th(1.8, 4.0, 4))
  updated <- .merge_threshold_lists(old_th, th_touched, touched_markers = "CD4")

  expect_identical(updated$CD3,  old_th$CD3)
  expect_identical(updated$CD19, old_th$CD19)
  expect_identical(updated$CD4,  th_touched$CD4)
  expect_false(identical(updated$CD4, old_th$CD4))
})

test_that("a brand-new marker is added without disturbing existing ones", {
  old_th <- list(CD3 = .mk_th(1.0, 3.0, 1), CD4 = .mk_th(1.2, 3.2, 2))
  th_touched <- list(CD16 = .mk_th(2.0, 5.0, 5))
  updated <- .merge_threshold_lists(old_th, th_touched, touched_markers = "CD16")

  expect_true("CD16" %in% names(updated))
  expect_equal(length(updated), 3L)
  expect_identical(updated$CD3, old_th$CD3)
  expect_identical(updated$CD4, old_th$CD4)
})

test_that("merge errors loudly if th_touched is missing a declared touched marker", {
  old_th <- list(CD3 = .mk_th(1.0, 3.0, 1))
  expect_error(
    .merge_threshold_lists(old_th, list(CD4 = .mk_th(1, 1, 2)), touched_markers = c("CD4", "CD99")),
    "missing results for"
  )
})

test_that("audit table correctly separates new markers from updated ones and computes shift", {
  old_th <- list(CD4 = .mk_th(1.0, 3.0, 1), CD8 = .mk_th(0.9, 2.8, 2))
  updated_th <- list(CD4 = .mk_th(1.8, 4.0, 3), CD8 = .mk_th(0.9, 2.8, 2), CD16 = .mk_th(2.0, 5.0, 4))
  audit <- .build_update_audit(old_th, updated_th,
                               touched_markers = c("CD4", "CD8", "CD16"),
                               markers_brand_new = "CD16")

  expect_equal(nrow(audit), 3L)
  expect_equal(audit$status[audit$marker == "CD16"], "new_marker")
  expect_true(is.na(audit$old_mean_lower[audit$marker == "CD16"]))
  expect_equal(audit$status[audit$marker == "CD4"], "existing_marker_updated")
  ## CD4 was deliberately shifted a lot; CD8 was regenerated with the same
  ## seed/params and should show ~0 shift.
  expect_gt(abs(audit$shift_lower[audit$marker == "CD4"]), 0.3)
  expect_lt(abs(audit$shift_lower[audit$marker == "CD8"]), 0.01)
})
