###############################################################################
# reference_size_convergence.R
#
# Assesses how the framework's threshold estimates depend on the size of
# the reference, via a downsampling study showing how much reference data
# is actually needed.
#
# WHAT THIS SIMULATES: a single bimodal marker with 3000 cells pooled
# across 6 "datasets". downsample_reference_stability() rebuilds the
# bootstrap threshold at a sweep of subsample fractions (5% up to 100% of
# the cells), each fraction repeated 8 times to separate genuine
# convergence from single-draw noise, and compares every fraction's
# threshold estimate against the full-reference (100%) estimate.
#
# WHAT THIS DEMONSTRATES: a concrete, quantitative convergence curve --
# how much the threshold estimate typically differs from the "ground
# truth" full-reference value as a function of how much data went into it.
# This directly answers "how much reference is enough" with a number,
# rather than an assumption that more is always better without
# quantifying the return on additional data.
###############################################################################

suppressPackageStartupMessages(library(scProtXpress))

set.seed(20260831)
out_dir <- "reference_size_convergence_outputs"
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. SIMULATE A REALISTIC-SIZED REFERENCE (1 marker, 3000 cells, 6 datasets)
## ---------------------------------------------------------------------------
n_total <- 3000
mk <- c(rnorm(n_total * 0.6, 1, 0.6), rnorm(n_total * 0.4, 5, 0.7))
data_norm <- matrix(mk, nrow = 1, dimnames = list("CD_marker", NULL))
metadata <- data.frame(cell.id = paste0("c", seq_len(n_total)),
                       dataset = sample(paste0("ds", 1:6), n_total, replace = TRUE),
                       cell.type = "AllCells", stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 2. RUN THE DOWNSAMPLING STUDY
## ---------------------------------------------------------------------------
fractions <- c(0.05, 0.10, 0.25, 0.50, 0.75, 1.0)
cat("Running downsampling study:", length(fractions), "fractions x 8 replicates each",
   "(the 100% fraction has only 1 'replicate' -- it IS the reference)...\n\n")

result <- downsample_reference_stability(
  data_norm, metadata,
  fractions = fractions, subsample_by = "cells", n_replicates = 8,
  sample_n = 150, num.boot.sample = 100, fixseed = 1,
  bootstrap_stratify = "none", threshold_method = "rectangle", ncores = 1L
)

cat("=== Raw per-replicate results (first 10 rows) ===\n")
print(head(result, 10), row.names = FALSE, digits = 3)

## ---------------------------------------------------------------------------
## 3. SUMMARIZE: mean +/- SD shift from the full reference, per fraction
## ---------------------------------------------------------------------------
summary_tbl <- aggregate(cbind(abs_shift_lower, abs_shift_upper) ~ fraction + n_cells_used,
                         data = result, FUN = function(x) c(mean = mean(x), sd = sd(x)))
## aggregate() with a 2-column FUN returns a matrix column -- flatten it
summary_tbl <- do.call(data.frame, summary_tbl)
colnames(summary_tbl) <- c("fraction", "n_cells_used", "mean_shift_lower", "sd_shift_lower",
                          "mean_shift_upper", "sd_shift_upper")
summary_tbl <- summary_tbl[order(summary_tbl$fraction), ]

cat("\n=== Convergence summary: mean/SD |shift from full reference| by fraction ===\n")
print(summary_tbl, row.names = FALSE, digits = 3)
write.csv(summary_tbl, file.path(out_dir, "convergence_summary.csv"), row.names = FALSE)

## Practical takeaway: smallest fraction whose mean shift is within a
## chosen tolerance of the (necessarily zero, by construction) 100% row.
tolerance <- 0.10  # arbitrary but reasonable for illustration; set to whatever
                   # magnitude matters for your own marker's dynamic range
adequate <- summary_tbl[summary_tbl$fraction < 1.0 &
                        summary_tbl$mean_shift_lower < tolerance &
                        summary_tbl$mean_shift_upper < tolerance, ]
if (nrow(adequate) > 0) {
  best <- adequate[which.min(adequate$fraction), ]
  cat(sprintf("\nAt a tolerance of %.2f, fraction %.2f (%d cells) is already sufficient --\n",
             tolerance, best$fraction, best$n_cells_used))
  cat("additional data beyond that point yields diminishing returns for this marker.\n")
} else {
  cat(sprintf("\nNo fraction below 100%% meets the %.2f tolerance for this marker --\n", tolerance))
  cat("the full reference size tested here may still be adding real precision.\n")
}

## ---------------------------------------------------------------------------
## 4. PLOT: convergence curve with error bars
## ---------------------------------------------------------------------------
png(file.path(out_dir, "convergence_plot.png"), width = 950, height = 650, res = 120)
par(mar = c(5, 5, 3, 2))
plot(summary_tbl$fraction, summary_tbl$mean_shift_upper, type = "b", pch = 16, col = "steelblue",
    ylim = c(0, max(summary_tbl$mean_shift_upper + summary_tbl$sd_shift_upper, na.rm = TRUE) * 1.1),
    xlab = "Fraction of full reference used", ylab = "Mean |threshold shift| from full reference",
    main = "Reference-size convergence (Upper threshold)")
arrows(summary_tbl$fraction, summary_tbl$mean_shift_upper - summary_tbl$sd_shift_upper,
      summary_tbl$fraction, summary_tbl$mean_shift_upper + summary_tbl$sd_shift_upper,
      angle = 90, code = 3, length = 0.05, col = "steelblue")
abline(h = tolerance, lty = 2, col = "grey40")
text(0.15, tolerance, "illustrative tolerance", pos = 3, cex = 0.8, col = "grey40")
dev.off()
cat("\nPlot written to", file.path(out_dir, "convergence_plot.png"), "\n")

## ---------------------------------------------------------------------------
## ADAPT TO YOUR OWN DATA
## ---------------------------------------------------------------------------
## Replace section 1 with your real reference NormData/Metadata (needs a
## 'dataset' column only if you want subsample_by = "datasets" instead of
## "cells" -- that variant tests a different question: does the reference
## need many independent STUDIES, not just many cells). Run this once per
## marker of interest, or loop over a representative subset of markers if
## running it for all ~700 would be too slow -- markers that already
## converge by a low fraction in this analysis are the ones for which
## additional reference data would matter least; markers that don't
## converge even at 100% may be under-supported regardless of fraction
## (compare against leave_one_dataset_out_stability()'s per-marker results
## for those specifically).
