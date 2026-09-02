###############################################################################
# patient_level_classifier.R
#
# Demonstrates a held-out, patient-level performance evaluation for
# clinical/diagnostic translation claims, rather than relying only on
# cell-population calls.
#
# WHAT THIS SIMULATES: a cohort of patients (default 30, a realistic
# disease-cohort scale) with a per-patient "biomarker proportion" feature
# -- e.g. the fraction of a population's cells called positive for a
# candidate marker, aggregated per patient, exactly what
# discrete.expr.levels_donor()'s by-donor breakdown would give you on real
# data. Two scenarios:
#   1. A genuine, moderate effect (disease patients have a higher true
#      proportion than healthy controls -- not a huge, unrealistic
#      separation).
#   2. A null scenario (no true difference between groups) at the SAME
#      cohort size, included specifically to demonstrate why the
#      permutation-null calibration matters (see fit_patient_classifier()'s
#      documentation for why raw small-sample cross-validated AUC is not
#      centered on 0.5 and should not be compared to 0.5 directly).
#
# WHAT THIS DEMONSTRATES: fit_patient_classifier() with leave-one-out
# cross-validation and n_permutations correctly distinguishes the real-
# effect scenario (significant permutation p-value) from the null scenario
# (non-significant), using a defensible, proof-of-concept-appropriate
# statistical framework rather than an unvalidated AUC number.
###############################################################################

suppressPackageStartupMessages(library(scProtXpress))

set.seed(20260831)
out_dir <- "patient_level_classifier_outputs"
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. SIMULATE A COHORT WITH A REAL, MODERATE EFFECT
## ---------------------------------------------------------------------------
n_patients <- 30
disease_status <- rep(c("Healthy", "Disease"), each = n_patients / 2)

## a realistic, moderate separation -- not the dramatic, unrealistic kind
## that would make any method look perfect.
biomarker_prop <- c(
  rbeta(n_patients / 2, 4, 16),   # Healthy: mean ~0.20
  rbeta(n_patients / 2, 7, 13)    # Disease: mean ~0.35 -- a real but modest shift
)
patient_features <- data.frame(biomarker_prop = biomarker_prop)

cat("=== Cohort summary ===\n")
cat(sprintf("n = %d patients (%d Healthy, %d Disease)\n", n_patients, n_patients/2, n_patients/2))
cat(sprintf("Mean biomarker proportion: Healthy = %.3f, Disease = %.3f\n",
           mean(biomarker_prop[disease_status == "Healthy"]),
           mean(biomarker_prop[disease_status == "Disease"])))

## ---------------------------------------------------------------------------
## 2. FIT WITH LEAVE-ONE-OUT CV + PERMUTATION-NULL CALIBRATION
## ---------------------------------------------------------------------------
cat("\nFitting leave-one-out cross-validated classifier with 300 label\n",
   "permutations for null calibration (this takes a little while)...\n", sep = "")

res_signal <- suppressMessages(fit_patient_classifier(
  patient_features, disease_status, cv_method = "loo",
  n_permutations = 300, fixseed = 1
))

cat("\n=== Real-effect scenario ===\n")
cat(sprintf("Observed AUC: %.3f\n", res_signal$auc))
cat(sprintf("Sensitivity: %.3f, Specificity: %.3f (at threshold 0.5)\n", res_signal$sensitivity, res_signal$specificity))
cat(sprintf("Null distribution (300 label permutations): mean = %.3f, 95th pct = %.3f\n",
           mean(res_signal$null_aucs), quantile(res_signal$null_aucs, 0.95)))
cat(sprintf("Permutation p-value: %.4f\n", res_signal$perm_p_value))

## ---------------------------------------------------------------------------
## 3. NULL SCENARIO AT THE SAME COHORT SIZE, FOR COMPARISON
## ---------------------------------------------------------------------------
biomarker_prop_null <- rbeta(n_patients, 5, 15)   # SAME distribution regardless of disease_status
patient_features_null <- data.frame(biomarker_prop = biomarker_prop_null)

res_null <- suppressMessages(fit_patient_classifier(
  patient_features_null, disease_status, cv_method = "loo",
  n_permutations = 300, fixseed = 2
))

cat("\n=== Null scenario (no true difference, same cohort size) ===\n")
cat(sprintf("Observed AUC: %.3f  (note: NOT necessarily close to 0.5 -- see below)\n", res_null$auc))
cat(sprintf("Null distribution: mean = %.3f, 95th pct = %.3f\n",
           mean(res_null$null_aucs), quantile(res_null$null_aucs, 0.95)))
cat(sprintf("Permutation p-value: %.4f\n", res_null$perm_p_value))

cat("\nInterpretation: comparing the two scenarios' RAW AUC values to each\n",
   "other, or to 0.5, would not reliably distinguish them -- both are small-\n",
   "cohort LOO-CV estimates and can land anywhere in a wide range under\n",
   "either the null or a real effect (see fit_patient_classifier()'s\n",
   "documented small-sample bias). The permutation p-value is the correct\n",
   "comparison: it is small for the real-effect scenario and large for the\n",
   "null scenario, because each is compared against ITS OWN empirical null\n",
   "distribution for its own data, not against an assumed reference point.\n",
   "Reporting the permutation p-value, with the AUC only as a secondary/\n",
   "descriptive number and the whole analysis framed as proof-of-concept, is\n",
   "what makes a patient-level clinical-translation claim defensible on a\n",
   "real cohort of this size.\n", sep = "")

results_summary <- data.frame(
  scenario = c("real_moderate_effect", "null_same_cohort_size"),
  n_patients = n_patients,
  observed_auc = c(res_signal$auc, res_null$auc),
  null_mean_auc = c(mean(res_signal$null_aucs), mean(res_null$null_aucs)),
  perm_p_value = c(res_signal$perm_p_value, res_null$perm_p_value)
)
write.csv(results_summary, file.path(out_dir, "results_summary.csv"), row.names = FALSE)
print(results_summary, row.names = FALSE, digits = 3)

## ---------------------------------------------------------------------------
## 4. PLOT: observed AUC against its own null distribution, both scenarios
## ---------------------------------------------------------------------------
png(file.path(out_dir, "null_distributions_plot.png"), width = 1000, height = 650, res = 120)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

hist(res_signal$null_aucs, breaks = 25, col = "grey80", border = "white",
    xlim = c(0, 1), xlab = "AUC", main = "Real-effect scenario",
    sub = sprintf("observed AUC = %.2f, p = %.3f", res_signal$auc, res_signal$perm_p_value), cex.sub = 0.9)
abline(v = res_signal$auc, col = "firebrick", lwd = 2)
legend("topleft", bty = "n", legend = "observed", col = "firebrick", lwd = 2, cex = 0.8)

hist(res_null$null_aucs, breaks = 25, col = "grey80", border = "white",
    xlim = c(0, 1), xlab = "AUC", main = "Null scenario",
    sub = sprintf("observed AUC = %.2f, p = %.3f", res_null$auc, res_null$perm_p_value), cex.sub = 0.9)
abline(v = res_null$auc, col = "firebrick", lwd = 2)
legend("topleft", bty = "n", legend = "observed", col = "firebrick", lwd = 2, cex = 0.8)

dev.off()
cat("\nPlot written to", file.path(out_dir, "null_distributions_plot.png"), "\n")

## ---------------------------------------------------------------------------
## ADAPT TO YOUR OWN DATA
## ---------------------------------------------------------------------------
## Replace `patient_features`/`disease_status` with your real per-patient
## data: aggregate discrete.expr.levels_donor()'s by-donor breakdown (or
## Create.sig.frame_by_donor()'s raw n0/n05/n1 counts) into one row per
## patient -- e.g. proportion of cells called "high" for your candidate
## biomarker(s) -- and your real disease/healthy labels. Multiple candidate
## biomarkers can go in as multiple columns of patient_features
## simultaneously. Always report n_permutations > 0 and the permutation
## p-value alongside the AUC for a cohort this size; do not report the raw
## AUC alone.
