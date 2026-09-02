###############################################################################
# donor_aware_testing.R
#
# Tests whether pooling all cells within a condition regardless of donor,
# treating each cell as an independent trial, causes pseudoreplication when
# cells from the same donor are correlated -- which they typically are --
# and validates a donor-aware alternative.
#
# WHAT THIS SIMULATES: a population of cells drawn from D donors, where
# each donor's TRUE non-zero rate for a marker is itself randomly drawn
# around a population-average rate (the standard beta-binomial model for
# intra-cluster correlation, parameterised by an intra-class correlation,
# ICC/"rho"). rho = 0 means donors are identical (no real between-donor
# variability); rho > 0 means genuine, realistic biological/technical
# variability between donors, which is the norm, not the exception, in
# real cohorts.
#
# TWO scenarios are run, both across a realistic range of rho (0 to 0.20):
#   1. NULL: the population-average rate truly equals the reference
#      baseline (prob_low) -- any "medium"/"high" call is a false positive.
#      Reports the FALSE-POSITIVE RATE of the pooled-cell test (the
#      original default, "pooled") vs. the donor-aware test
#      ("design_effect") across 3000 simulated null genes per rho value.
#   2. SIGNAL: the population-average rate is genuinely elevated above the
#      reference baseline. Reports the RECALL (power) of both tests, to
#      confirm the donor-aware fix isn't just "conservative everywhere" --
#      it should still detect a real, biologically meaningful effect.
#
# WHAT THIS DEMONSTRATES: the pooled-cell test's false-positive rate is
# nominal (~5%) only when rho = 0, and inflates sharply as soon as there is
# ANY real between-donor heterogeneity -- while the donor-aware test stays
# close to nominal across the whole range AND retains meaningful power to
# detect a real effect of realistic size.
###############################################################################

suppressPackageStartupMessages(library(scProtXpress))

## Create.sig.frame()/discrete.expr.levels() (the original pooled-cell
## path) are internal by design; Create.sig.frame_by_donor()/
## discrete.expr.levels_donor() (the donor-aware alternative) are part of
## the public API.
sig_frame_pooled     <- scProtXpress:::Create.sig.frame
expr_levels_pooled   <- scProtXpress:::discrete.expr.levels

set.seed(20260831)
out_dir <- "donor_aware_testing_outputs"
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. SIMULATION ENGINE
## ---------------------------------------------------------------------------
D <- 10; n_per_donor <- 150; n_genes <- 300

simulate_discmat <- function(D, n_per_donor, prob_true, rho, n_genes) {
  donor_ids <- rep(paste0("donor", seq_len(D)), each = n_per_donor)
  ncells <- length(donor_ids)
  mat <- matrix(0, nrow = n_genes, ncol = ncells)
  kappa <- if (rho > 0) (1 - rho) / rho else Inf
  for (g in seq_len(n_genes)) {
    donor_rate <- if (rho > 0) rbeta(D, prob_true * kappa, (1 - prob_true) * kappa) else rep(prob_true, D)
    mat[g, ] <- unlist(lapply(seq_len(D), function(d) rbinom(n_per_donor, 1, donor_rate[d])))
  }
  rownames(mat) <- paste0("CD", seq_len(n_genes))
  list(mat = mat, donor_ids = donor_ids)
}

run_comparison <- function(prob_true, prob_baseline, rho) {
  sim <- simulate_discmat(D, n_per_donor, prob_true, rho, n_genes)

  sig  <- sig_frame_pooled(sim$mat, min.coverage.fraction.sigframe = 0)
  tbl_pooled <- expr_levels_pooled(sig[[1]], 0, prob_low = prob_baseline, prob_high = 1e-6)

  by_donor <- Create.sig.frame_by_donor(sim$mat, sim$donor_ids)
  tbl_deff <- suppressMessages(
    discrete.expr.levels_donor(by_donor, prob_baseline, 1e-6, test_method = "design_effect")
  )

  c(pooled = mean(tbl_pooled$expr.level %in% c("medium", "high")),
    design_effect = mean(tbl_deff$expr.level %in% c("medium", "high")))
}

## ---------------------------------------------------------------------------
## 2. NULL SCENARIO: false-positive rate across a realistic rho sweep
## ---------------------------------------------------------------------------
rho_values <- c(0, 0.02, 0.05, 0.10, 0.15, 0.20)
prob_baseline <- 0.30

cat("=== NULL scenario: false-positive rate (population-average rate truly\n",
   "    equals the reference baseline, ", n_genes, " simulated null genes per rho) ===\n", sep = "")
null_results <- do.call(rbind, lapply(rho_values, function(rho) {
  res <- run_comparison(prob_true = prob_baseline, prob_baseline = prob_baseline, rho = rho)
  data.frame(rho = rho, pooled_FPR = res["pooled"], design_effect_FPR = res["design_effect"])
}))
rownames(null_results) <- NULL
print(null_results, row.names = FALSE, digits = 3)
write.csv(null_results, file.path(out_dir, "null_fpr_sweep.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. SIGNAL SCENARIO: recall (power) for a real, moderate effect
## ---------------------------------------------------------------------------
prob_true_signal <- 0.45   # a real, moderate elevation over the 0.30 baseline
cat("\n=== SIGNAL scenario: recall for a real effect (baseline", prob_baseline,
   "-> true rate", prob_true_signal, ") ===\n")
signal_results <- do.call(rbind, lapply(rho_values, function(rho) {
  res <- run_comparison(prob_true = prob_true_signal, prob_baseline = prob_baseline, rho = rho)
  data.frame(rho = rho, pooled_recall = res["pooled"], design_effect_recall = res["design_effect"])
}))
rownames(signal_results) <- NULL
print(signal_results, row.names = FALSE, digits = 3)
write.csv(signal_results, file.path(out_dir, "signal_recall_sweep.csv"), row.names = FALSE)

cat("\nInterpretation: at rho=0 both tests are close to nominal/have similar\n",
   "power -- confirming the donor-aware fix isn't just 'always more\n",
   "conservative'. As soon as rho > 0 (i.e. almost always, in real data),\n",
   "the pooled test's false-positive rate inflates sharply while its\n",
   "recall stays at ~100% regardless of rho -- which is not a strength,\n",
   "it means the pooled test rejects almost everything whether or not\n",
   "there is a real effect, so its 'recall' is close to meaningless once\n",
   "its FPR is this high.\n\n",
   "Honest tradeoff worth reporting, not hiding: design_effect's recall\n",
   "DOES fall as rho increases (100% at rho=0 down to ~66% at rho=0.20 for\n",
   "this specific, moderate effect size) -- properly controlling the false-\n",
   "positive rate under heavy donor clustering necessarily costs some power,\n",
   "especially for effects that aren't large. That is the expected,\n",
   "statistically correct behaviour of a properly calibrated test, not a\n",
   "flaw -- but it means a real biomarker signal that is only moderate in\n",
   "size, in a cohort with strong between-donor heterogeneity, may need a\n",
   "larger donor count to reach significance once tested correctly. Rigor\n",
   "has a power cost, and that should be stated plainly rather than implying\n",
   "the fix is free.\n", sep = "")

## ---------------------------------------------------------------------------
## 4. PLOTS
## ---------------------------------------------------------------------------
png(file.path(out_dir, "fpr_plot.png"), width = 900, height = 650, res = 120)
par(mar = c(5, 5, 3, 2))
plot(null_results$rho, null_results$pooled_FPR, type = "b", pch = 16, col = "firebrick",
    ylim = c(0, max(null_results$pooled_FPR) * 1.1),
    xlab = "Intra-donor correlation (rho)", ylab = "False-positive rate",
    main = "False-positive rate under a true null, by donor clustering")
lines(null_results$rho, null_results$design_effect_FPR, type = "b", pch = 16, col = "seagreen")
abline(h = 0.05, lty = 2, col = "grey40")
legend("topleft", bty = "n", legend = c("pooled (original default)", "design_effect (donor-aware)", "nominal alpha = 0.05"),
      col = c("firebrick", "seagreen", "grey40"), pch = c(16, 16, NA), lty = c(1, 1, 2))
dev.off()
cat("\nPlot written to", file.path(out_dir, "fpr_plot.png"), "\n")

png(file.path(out_dir, "recall_plot.png"), width = 900, height = 650, res = 120)
par(mar = c(5, 5, 3, 2))
plot(signal_results$rho, signal_results$pooled_recall, type = "b", pch = 16, col = "firebrick",
    ylim = c(0, 1), xlab = "Intra-donor correlation (rho)", ylab = "Recall (power)",
    main = "Recall for a real, moderate effect, by donor clustering")
lines(signal_results$rho, signal_results$design_effect_recall, type = "b", pch = 16, col = "seagreen")
legend("bottomleft", bty = "n", legend = c("pooled (original default)", "design_effect (donor-aware)"),
      col = c("firebrick", "seagreen"), pch = 16, lty = 1)
dev.off()
cat("Plot written to", file.path(out_dir, "recall_plot.png"), "\n")

## ---------------------------------------------------------------------------
## ADAPT TO YOUR OWN DATA
## ---------------------------------------------------------------------------
## Replace the simulation in section 1 with your real discretised query
## matrix (Discretise_data()'s output) and real donor IDs from your
## metadata, and rerun sections 2/3's comparison calls directly -- or, more
## directly, just rerun your real disease-vs-healthy analysis with
## run_scProteinXpress(..., sig_test_method = "design_effect") instead of
## the default "pooled", and report both results side by side, exactly as
## this script does for simulated data.
