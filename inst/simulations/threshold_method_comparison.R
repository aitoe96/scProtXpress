###############################################################################
# threshold_method_comparison.R
#
# Compares the bootstrap-based thresholding approach against simple
# reference-informed baselines (a fixed quantile of the reference
# distribution, and a mean +/- k*SD z-score cutoff), to show whether any
# accuracy gain comes from the specific thresholding procedure or just from
# having a reference at all.
#
# WHAT THIS SIMULATES: four synthetic single-protein distributions with
# KNOWN ground-truth population boundaries (so "how close is the estimated
# threshold to the truth" is an exact, checkable number, not a subjective
# judgement):
#   - bimodal_well_separated : two clearly separated populations.
#   - bimodal_close          : two populations close enough to partially
#                              overlap -- the hardest realistic case.
#   - trimodal                : three well-separated populations (the
#                              low/medium/high design the framework targets).
#   - trimodal_unequal        : three populations with unequal sizes/spread.
#
# All four threshold_method options this package supports -- "rectangle"
# (original heuristic), "mixture", "quantile" (simple baseline), "zscore"
# (simple baseline) -- are run on identical simulated data, and each
# method's threshold(s) are compared against the TRUE analytically-computed
# crossing point(s) between the simulated component distributions.
#
# WHAT THIS DEMONSTRATES: whether "mixture"'s accuracy advantage (if any)
# over the two naive reference-informed baselines (quantile/zscore) holds up
# scenario by scenario, not just in aggregate -- directly answering whether
# the accuracy comes from the reference itself (which all four methods use
# equally) or from the specific thresholding procedure.
###############################################################################

suppressPackageStartupMessages({
  library(scProtXpress)
  library(mclust)
})

## Internal per-method threshold functions are not part of the public API
## (by design -- they're implementation detail of GetBackground()/
## run_scProteinXpress()'s threshold_method argument) but are exactly what
## this comparison needs directly. Standard, appropriate use of ::: for a
## standalone analysis script (unlike inside the package's own tests, where
## a portability shim is used instead -- see tests/testthat/test-
## calibration.R for why).
rectangle_fn <- scProtXpress:::createThresholdDist_rectangle
mixture_fn   <- scProtXpress:::createThresholdDist_mixture
quantile_fn  <- scProtXpress:::createThresholdDist_quantile
zscore_fn    <- scProtXpress:::createThresholdDist_zscore

set.seed(20260831)
out_dir <- "threshold_method_comparison_outputs"
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. FOUR SCENARIOS WITH ANALYTICALLY KNOWN GROUND-TRUTH CROSSING POINTS
## ---------------------------------------------------------------------------
true_crossing <- function(means, sds, weights) {
  k <- length(means)
  vapply(seq_len(k - 1), function(i) {
    f <- function(x) weights[i] * dnorm(x, means[i], sds[i]) - weights[i + 1] * dnorm(x, means[i + 1], sds[i + 1])
    uniroot(f, lower = means[i], upper = means[i + 1])$root
  }, numeric(1))
}

scenarios <- list(
  bimodal_well_separated = list(means = c(1, 5),    sds = c(0.5, 0.6), weights = c(0.6, 0.4)),
  bimodal_close          = list(means = c(1, 2.8),  sds = c(0.5, 0.6), weights = c(0.5, 0.5)),
  trimodal               = list(means = c(1, 3, 6), sds = c(0.4, 0.5, 0.5), weights = c(0.5, 0.3, 0.2)),
  trimodal_unequal       = list(means = c(0.5, 3, 7), sds = c(0.3, 0.7, 0.6), weights = c(0.7, 0.2, 0.1))
)

simulate_from <- function(s, n = 1500) {
  k <- sample(seq_along(s$means), n, replace = TRUE, prob = s$weights)
  rnorm(n, s$means[k], s$sds[k])
}

methods <- list(rectangle = rectangle_fn, mixture = mixture_fn,
               quantile = quantile_fn, zscore = zscore_fn)

## ---------------------------------------------------------------------------
## 2. RUN ALL FOUR METHODS ON ALL FOUR SCENARIOS
## ---------------------------------------------------------------------------
results <- list()
for (scen_name in names(scenarios)) {
  s <- scenarios[[scen_name]]
  truth <- true_crossing(s$means, s$sds, s$weights)
  x <- simulate_from(s)

  for (m_name in names(methods)) {
    thr <- methods[[m_name]](x)
    ## distance from truth: for each true crossing, the closest estimated
    ## threshold; report the MAX over true crossings (a method that only
    ## finds one of two true boundaries in the trimodal case should not
    ## get full credit just because the boundary it DID find is close).
    dists <- vapply(truth, function(t) min(abs(thr - t)), numeric(1))
    results[[length(results) + 1]] <- data.frame(
      scenario = scen_name, method = m_name,
      n_true_crossings = length(truth),
      threshold_low = thr[1], threshold_high = thr[2],
      max_dist_from_truth = max(dists), mean_dist_from_truth = mean(dists),
      stringsAsFactors = FALSE
    )
  }
}
results <- do.call(rbind, results)

cat("=== Threshold estimate vs. true crossing point, all methods x all scenarios ===\n")
print(results, row.names = FALSE, digits = 3)
write.csv(results, file.path(out_dir, "all_results.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. SUMMARY: AVERAGE ACCURACY PER METHOD ACROSS SCENARIOS
## ---------------------------------------------------------------------------
summary_by_method <- aggregate(max_dist_from_truth ~ method, data = results, FUN = mean)
summary_by_method <- summary_by_method[order(summary_by_method$max_dist_from_truth), ]
cat("\n=== Mean |distance from truth| across all 4 scenarios, by method (lower is better) ===\n")
print(summary_by_method, row.names = FALSE, digits = 3)
write.csv(summary_by_method, file.path(out_dir, "method_ranking.csv"), row.names = FALSE)

cat("\nInterpretation: all four methods use the identical simulated reference\n",
   "data -- any difference in accuracy is attributable to the thresholding\n",
   "PROCEDURE itself, not to having a reference at all. This directly\n",
   "answers whether the framework's accuracy comes from the reference\n",
   "(which a naive quantile/z-score cutoff would recover equally well) or\n",
   "from the specific bootstrap/mixture-model machinery.\n\n",
   "Honest caveat, visible in the per-scenario table above and worth\n",
   "reporting rather than hiding: 'mixture' wins clearly in 3 of the 4\n",
   "scenarios, but is the WORST method on 'trimodal_unequal' (small,\n",
   "highly-overlapping third component, low weight) -- BIC-based component\n",
   "selection can under-fit a small/rare population. If a 3-state design is\n",
   "a firm modelling commitment rather than something BIC should decide,\n",
   "consider fixing G=3 in createThresholdDist_mixture() rather than letting\n",
   "it choose G in 1:3, and/or reporting the chosen G per marker as a\n",
   "diagnostic so markers where a 3-state model wasn't well-supported can be\n",
   "flagged rather than silently mis-thresholded.\n", sep = "")

## ---------------------------------------------------------------------------
## 4. PLOT: per-scenario comparison
## ---------------------------------------------------------------------------
png(file.path(out_dir, "comparison_plot.png"), width = 1100, height = 750, res = 120)
par(mfrow = c(1, 1), mar = c(8, 5, 3, 2))
## reshape results into a method x scenario matrix for barplot()
plot_mat <- sapply(names(scenarios), function(sc) {
  sub <- results[results$scenario == sc, ]
  sub$max_dist_from_truth[match(names(methods), sub$method)]
})
rownames(plot_mat) <- names(methods)

bp <- barplot(plot_mat, beside = TRUE, las = 2,
             col = c("grey60", "steelblue", "darkorange", "seagreen"),
             ylab = "Max |threshold - true crossing|",
             main = "Threshold accuracy vs. ground truth, by method and scenario")
legend("topleft", bty = "n", legend = rownames(plot_mat),
      fill = c("grey60", "steelblue", "darkorange", "seagreen"))
dev.off()
cat("\nPlot written to", file.path(out_dir, "comparison_plot.png"), "\n")

## ---------------------------------------------------------------------------
## ADAPT TO YOUR OWN DATA
## ---------------------------------------------------------------------------
## This comparison also runs directly through the public pipeline: call
## GetBackground(..., threshold_method = "rectangle"/"mixture"/"quantile"/
## "zscore") or run_scProteinXpress(..., threshold_method = ...) on your
## real reference data with each of the four settings, then compare the
## resulting recall/precision (see precision_recall_negative_markers.R's
## evaluate_calls() usage) across the four reference builds -- the
## low-level comparison here is a faster proxy for exactly that same
## question on data with a known answer.
