###############################################################################
# precision_recall_negative_markers.R
#
# Demonstrates why a recall-only benchmark can be misleading, and why
# precision, false-positive rate, and F1 -- computed using known NEGATIVE
# marker/cell-type pairs, not just positive ones -- should also be reported.
#
# WHAT THIS SIMULATES: a marker panel benchmark of 100 markers x 5 cell
# types (500 marker/cell-type pairs total). Each marker has exactly one
# true "canonical" cell type (a true-positive pair); the other 4 cell
# types are true negatives for that marker (400 true-negative pairs in
# total -- the known negative pairs a recall-only report omits). A caller
# is simulated with realistic, plausible parameters: it correctly calls
# 92% of true positives (a good-looking recall) and incorrectly calls 12%
# of true negatives as positive (a realistic, non-zero false-positive rate
# -- e.g. from cross-reactivity, ambient signal, or borderline expression).
#
# WHAT THIS DEMONSTRATES: with those parameters, recall alone (0.90ish)
# looks like unambiguously strong performance. Precision, computed using
# the true-negative pairs the same benchmark already implies, tells a
# different story -- roughly a third of "positive" calls turn out to be
# false positives once negatives are actually counted. A recall-only
# report cannot distinguish "high performance" from "high recall achieved
# by calling almost everything positive."
###############################################################################

suppressPackageStartupMessages(library(scProtXpress))

set.seed(20260831)
out_dir <- "precision_recall_negative_markers_outputs"
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. SIMULATE THE BENCHMARK: markers, cell types, and a plausible caller
## ---------------------------------------------------------------------------
n_markers   <- 100
n_celltypes <- 5
true_positive_rate  <- 0.92   # the caller's real sensitivity
false_positive_rate <- 0.12   # the caller's real false-positive rate on true negatives

markers   <- paste0("CD", seq_len(n_markers))
celltypes <- paste0("CT", seq_len(n_celltypes))

## each marker has exactly one true ("canonical") cell type -- the other
## n_celltypes - 1 are true NEGATIVES for that marker.
true_celltype <- sample(celltypes, n_markers, replace = TRUE)
ground_truth <- data.frame(marker = markers, cell.type = true_celltype, stringsAsFactors = FALSE)

grid <- expand.grid(marker = markers, cell.type = celltypes, stringsAsFactors = FALSE)
is_tp_pair <- paste(grid$marker, grid$cell.type) %in% paste(ground_truth$marker, ground_truth$cell.type)
grid$called_positive <- NA_integer_
grid$called_positive[is_tp_pair]  <- rbinom(sum(is_tp_pair),  1, true_positive_rate)
grid$called_positive[!is_tp_pair] <- rbinom(sum(!is_tp_pair), 1, false_positive_rate)
grid$expr.level <- ifelse(grid$called_positive == 1, "high", "low")

cat("Simulated", nrow(grid), "marker x cell-type pairs (", sum(is_tp_pair),
   "true positives,", sum(!is_tp_pair), "true negatives).\n\n")

## ---------------------------------------------------------------------------
## 2. WHAT A RECALL-ONLY REPORT WOULD SHOW
## ---------------------------------------------------------------------------
recall_only <- mean(grid$called_positive[is_tp_pair])
cat(sprintf("Recall-only summary: recall = %.2f\n", recall_only))
cat("On its own, this looks like strong, unambiguous performance.\n\n")

## ---------------------------------------------------------------------------
## 3. THE FULL PICTURE: precision / recall / F1 / FPR via evaluate_calls()
## ---------------------------------------------------------------------------
full_metrics <- evaluate_calls(grid, ground_truth, positive_levels = "high")
cat("=== Full confusion-matrix metrics (using the true-negative pairs too) ===\n")
print(full_metrics, row.names = FALSE, digits = 3)

cat(sprintf(
  "\nInterpretation: recall alone (%.2f) looked strong. Precision (%.2f) shows\n",
  full_metrics$recall, full_metrics$precision))
cat(sprintf(
  "that %d of the %d calls flagged 'positive' were actually false positives --\n",
  full_metrics$FP, full_metrics$TP + full_metrics$FP))
cat(sprintf(
  "roughly %.0f%% of positive calls are wrong, at a false-positive rate of %.2f\n",
  100 * full_metrics$FP / (full_metrics$TP + full_metrics$FP), full_metrics$fpr))
cat("on the true-negative pairs. F1 (which balances both) is", round(full_metrics$f1, 2),
   "-- meaningfully lower than recall alone suggested.\n\n", sep = " ")

write.csv(full_metrics, file.path(out_dir, "full_metrics.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 4. SENSITIVITY CHECK: how does the recall-vs-precision gap change with
##    the true false-positive rate? (shows the demonstrated gap isn't an
##    artifact of one specific FPR choice)
## ---------------------------------------------------------------------------
fpr_sweep <- c(0.02, 0.05, 0.08, 0.12, 0.18, 0.25)
sweep_results <- do.call(rbind, lapply(fpr_sweep, function(fpr_true) {
  set.seed(which(fpr_sweep == fpr_true))
  g <- grid
  g$called_positive[is_tp_pair]  <- rbinom(sum(is_tp_pair),  1, true_positive_rate)
  g$called_positive[!is_tp_pair] <- rbinom(sum(!is_tp_pair), 1, fpr_true)
  g$expr.level <- ifelse(g$called_positive == 1, "high", "low")
  m <- evaluate_calls(g, ground_truth, positive_levels = "high")
  data.frame(true_fpr = fpr_true, recall = m$recall, precision = m$precision, f1 = m$f1)
}))
cat("=== How the gap changes with the caller's true false-positive rate ===\n")
print(sweep_results, row.names = FALSE, digits = 3)
write.csv(sweep_results, file.path(out_dir, "fpr_sweep.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. PLOT
## ---------------------------------------------------------------------------
png(file.path(out_dir, "metrics_plot.png"), width = 900, height = 650, res = 120)
par(mar = c(5, 5, 3, 2))
metrics_vec <- c(recall = full_metrics$recall, precision = full_metrics$precision,
                 f1 = full_metrics$f1, fpr = full_metrics$fpr)
bp <- barplot(metrics_vec, ylim = c(0, 1.1), col = c("seagreen", "firebrick", "steelblue", "grey50"),
             ylab = "Value", main = "Recall alone vs. the full picture")
text(bp, metrics_vec, labels = round(metrics_vec, 2), pos = 3)
dev.off()
cat("\nPlot written to", file.path(out_dir, "metrics_plot.png"), "\n")

png(file.path(out_dir, "fpr_sweep_plot.png"), width = 900, height = 650, res = 120)
par(mar = c(5, 5, 3, 2))
plot(sweep_results$true_fpr, sweep_results$recall, type = "b", pch = 16, col = "seagreen",
    ylim = c(0, 1), xlab = "True false-positive rate of the caller",
    ylab = "Value", main = "Recall stays flat; precision/F1 fall as FPR rises")
lines(sweep_results$true_fpr, sweep_results$precision, type = "b", pch = 16, col = "firebrick")
lines(sweep_results$true_fpr, sweep_results$f1, type = "b", pch = 16, col = "steelblue")
legend("bottomleft", bty = "n", legend = c("recall", "precision", "F1"),
      col = c("seagreen", "firebrick", "steelblue"), pch = 16, lty = 1)
dev.off()
cat("Plot written to", file.path(out_dir, "fpr_sweep_plot.png"), "\n")

## ---------------------------------------------------------------------------
## ADAPT TO YOUR OWN DATA
## ---------------------------------------------------------------------------
## Replace `ground_truth` with your real curated positive-marker list
## (data.frame(marker=, cell.type=)) and `grid`/`expr.level` with the real
## scProtXpress calls across all profiled cell types for the same markers
## (not just their canonical one -- you need calls for the true-negative
## pairs too, which is the part a recall-only report omits). evaluate_calls()
## then computes the same metrics on the real data unchanged. If some
## marker/cell-type pairs are genuinely untested/unknown rather than true
## negatives, exclude those rows before calling evaluate_calls() so they
## don't get miscounted as true negatives.
