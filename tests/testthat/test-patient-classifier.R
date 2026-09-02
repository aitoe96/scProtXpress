## Tests for fit_patient_classifier().
##
## The headline finding this function's design responds to: LOO/k-fold
## cross-validated AUC on small samples with weak/no signal is NOT centered
## on 0.5 (see the function's own docstring). That is why n_permutations
## exists and why these tests check the PERMUTATION P-VALUE's calibration,
## not just the raw AUC's proximity to 0.5.

## Internal-function access shim: see the identical pattern in
## test-calibration.R for why this is needed.
if (requireNamespace("scProtXpress", quietly = TRUE) &&
    !exists(".compute_auc", mode = "function")) {
  .compute_auc <- getFromNamespace(".compute_auc", "scProtXpress")
}

.make_classifier_test_data <- function(seed, n = 40, signal = c("strong", "none")) {
  signal <- match.arg(signal)
  set.seed(seed)
  labels <- rep(c(0, 1), each = n / 2)
  feat <- if (identical(signal, "strong")) {
    data.frame(biomarker_prop = c(rnorm(n / 2, 0.2, 0.08), rnorm(n / 2, 0.6, 0.08)))
  } else {
    data.frame(biomarker_prop = rnorm(n, 0.4, 0.1))
  }
  list(features = feat, labels = labels)
}

test_that("a strong true signal gives high AUC and a small permutation p-value", {
  sim <- .make_classifier_test_data(1, signal = "strong")
  res <- fit_patient_classifier(sim$features, sim$labels, cv_method = "loo",
                                n_permutations = 100, fixseed = 5)
  expect_gt(res$auc, 0.9)
  expect_lt(res$perm_p_value, 0.05)
  expect_equal(length(res$null_aucs), 100L)
})

test_that("a null (no-signal) scenario gives a non-significant permutation p-value regardless of the raw AUC's distance from 0.5", {
  sim <- .make_classifier_test_data(2, signal = "none")
  res <- fit_patient_classifier(sim$features, sim$labels, cv_method = "loo",
                                n_permutations = 100, fixseed = 6)
  ## This is the crux of the whole fix: do NOT assert res$auc is close to
  ## 0.5 (confirmed it often isn't, for reasons documented in the
  ## function). Assert instead that the permutation p-value -- which
  ## compares the observed AUC to the empirical null for this exact
  ## procedure -- correctly calls it non-significant.
  expect_gt(res$perm_p_value, 0.10)
})

test_that("permutation results are reproducible given the same fixseed", {
  sim <- .make_classifier_test_data(3, signal = "strong")
  res1 <- fit_patient_classifier(sim$features, sim$labels, n_permutations = 30, fixseed = 42)
  res2 <- fit_patient_classifier(sim$features, sim$labels, n_permutations = 30, fixseed = 42)
  expect_identical(res1$null_aucs, res2$null_aucs)
  expect_identical(res1$perm_p_value, res2$perm_p_value)
})

test_that("factor/character labels are coerced consistently with logical/0-1 labels", {
  sim <- .make_classifier_test_data(4, signal = "strong")
  res_numeric <- fit_patient_classifier(sim$features, sim$labels, cv_method = "loo")
  labels_factor <- factor(ifelse(sim$labels == 1, "Disease", "Healthy"))
  res_factor <- suppressMessages(fit_patient_classifier(sim$features, labels_factor, cv_method = "loo"))
  expect_equal(res_factor$auc, res_numeric$auc)
})

test_that("kfold and loo both run without error and kfold requires a valid k", {
  sim <- .make_classifier_test_data(5, signal = "strong")
  res_kfold <- fit_patient_classifier(sim$features, sim$labels, cv_method = "kfold", k = 5, fixseed = 1)
  expect_true(is.numeric(res_kfold$auc))
  expect_error(
    fit_patient_classifier(sim$features, sim$labels, cv_method = "kfold", k = 1),
    "k must be between"
  )
  expect_error(
    fit_patient_classifier(sim$features, sim$labels, cv_method = "kfold", k = 1000),
    "k must be between"
  )
})

test_that("clear errors for malformed inputs", {
  sim <- .make_classifier_test_data(6, signal = "strong")
  expect_error(fit_patient_classifier(sim$features[1:3, , drop = FALSE], sim$labels[1:3]), "at least 4")
  expect_error(fit_patient_classifier(sim$features, sim$labels[1:5]), "length")
  expect_error(fit_patient_classifier(sim$features, rep(1, nrow(sim$features))), "two classes")
  expect_error(fit_patient_classifier(sim$features, sim$labels, n_permutations = -1), "non-negative")
})

test_that("the AUC helper matches a hand-computed value on a tiny known example", {
  ## scores: 0.9, 0.8, 0.3, 0.2 for labels 1,1,0,0 -- perfect separation, AUC = 1
  expect_equal(.compute_auc(c(0.9, 0.8, 0.3, 0.2), c(1, 1, 0, 0)), 1)
  ## Completely reversed ranking -- AUC = 0
  expect_equal(.compute_auc(c(0.2, 0.3, 0.8, 0.9), c(1, 1, 0, 0)), 0)
  ## Ties split evenly -- AUC = 0.5
  expect_equal(.compute_auc(c(0.5, 0.5, 0.5, 0.5), c(1, 1, 0, 0)), 0.5)
})
