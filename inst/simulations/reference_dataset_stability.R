###############################################################################
# reference_dataset_stability.R
#
# Assesses how sensitive protein-specific reference thresholds are to the
# composition of the reference compendium, via leave-one-dataset-out
# analysis.
#
# WHAT THIS SIMULATES: a reference compendium of 6 "datasets" (studies)
# contributing cells to a shared reference. Three markers with different,
# known ground-truth behaviour:
#   - CD_stable   : the SAME true bimodal distribution in every dataset
#                   (an idealised, well-behaved marker).
#   - CD_batch    : normal in 5 datasets, but one dataset ("ds6") has a
#                   substantial positive-population shift -- simulating a
#                   failed antibody titration, a different clone, or a
#                   genuine batch effect in that one study.
#   - CD_rare     : only measured in 2 of the 6 datasets (a marker that
#                   happens to be in only a couple of the contributing
#                   panels) -- included to show the analysis correctly
#                   excludes markers without enough independent support
#                   rather than reporting a falsely-confident number for
#                   them.
#
# WHAT THIS DEMONSTRATES: leave_one_dataset_out_stability() correctly (a)
# excludes CD_rare, (b) reports CD_stable as materially more stable than
# CD_batch, and (c) the per-dataset breakdown pinpoints ds6 specifically as
# the source of CD_batch's instability. Swap in your own real, already-
# normalized reference data (see "ADAPT TO YOUR OWN DATA" below) to run
# this on your own reference.
#
# Output: a summary table (printed and written to
# stability_summary.csv) plus a bar plot of each marker's largest
# single-dataset threshold shift.
###############################################################################

suppressPackageStartupMessages(library(scProtXpress))

set.seed(20260831)
out_dir <- "reference_dataset_stability_outputs"
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. SIMULATE THE REFERENCE (replace this block with your own data -- see
##    "ADAPT TO YOUR OWN DATA" at the bottom)
## ---------------------------------------------------------------------------
datasets   <- paste0("ds", 1:6)
n_per_ds   <- 250
dataset_vec <- rep(datasets, each = n_per_ds)
n_total    <- length(dataset_vec)
cell_ids   <- paste0("cell_", seq_len(n_total))

simulate_bimodal <- function(n, neg_mean, pos_mean, pos_frac = 0.35, sd = 0.55) {
  is_pos <- runif(n) < pos_frac
  ifelse(is_pos, rnorm(n, pos_mean, sd), rnorm(n, neg_mean, sd))
}

CD_stable <- unlist(lapply(datasets, function(d) simulate_bimodal(n_per_ds, 1.0, 5.0)))

CD_batch <- unlist(lapply(datasets, function(d) {
  if (d == "ds6") simulate_bimodal(n_per_ds, 1.0, 8.5)   # antibody/batch shift in ds6 only
  else            simulate_bimodal(n_per_ds, 1.0, 5.0)
}))

CD_rare <- rep(NA_real_, n_total)
in_two_datasets <- dataset_vec %in% c("ds1", "ds2")
CD_rare[in_two_datasets] <- simulate_bimodal(sum(in_two_datasets), 0.8, 3.2)

data_norm <- rbind(CD_stable = CD_stable, CD_batch = CD_batch, CD_rare = CD_rare)
colnames(data_norm) <- cell_ids
metadata <- data.frame(cell.id = cell_ids, dataset = dataset_vec,
                       cell.type = "AllCells", stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 2. RUN THE STABILITY ANALYSIS
## ---------------------------------------------------------------------------
cat("Running leave-one-dataset-out stability analysis (this rebuilds the\n",
   "bootstrap threshold distribution once per dataset held out)...\n\n", sep = "")

result <- leave_one_dataset_out_stability(
  data_norm, metadata,
  min_datasets_per_marker = 3,
  sample_n = 150, num.boot.sample = 300, fixseed = 1,
  ## "rectangle" (the original heuristic) kept deliberately simple here so
  ## the printed threshold VALUES are easy to sanity-check by eye against
  ## the simulated means -- this script is about DATASET stability, not
  ## about which threshold_method is best (see threshold_method_comparison.R
  ## for that). Swap to "mixture" if that's what your real reference build
  ## uses.
  bootstrap_stratify = "none", threshold_method = "rectangle", ncores = 1L
)

summary_tbl <- result$stability_summary
summary_tbl <- summary_tbl[order(-summary_tbl$max_abs_shift_upper), ]

cat("=== Stability summary (sorted by largest single-dataset threshold shift) ===\n")
print(summary_tbl, row.names = FALSE, digits = 3)

write.csv(summary_tbl, file.path(out_dir, "stability_summary.csv"), row.names = FALSE)

## Which dataset is responsible for CD_batch's instability? Pull it out of
## the per-dataset results directly.
if ("CD_batch" %in% names(result$per_dataset[[1]]) || TRUE) {
  per_ds_upper <- vapply(result$per_dataset, function(th) {
    if (is.null(th[["CD_batch"]])) NA_real_ else mean(th[["CD_batch"]][, "Upper"])
  }, numeric(1))
  full_upper <- mean(result$full_reference[["CD_batch"]][, "Upper"])
  cat("\n=== CD_batch's Upper threshold when EACH dataset is held out ===\n")
  print(round(sort(per_ds_upper), 3))
  cat(sprintf("\nFull-reference (all 6 datasets) Upper threshold: %.3f\n", full_upper))
  cat("Interpretation: ds6 is the one dataset with a simulated batch effect on\n",
     "this marker. It correctly shows up as the clear outlier above -- holding\n",
     "it out gives a noticeably different (here, lower) threshold estimate than\n",
     "holding out any of the other 5 datasets, which barely move the estimate\n",
     "at all. That asymmetry is the diagnostic signal: it correctly localizes\n",
     "the instability to one specific dataset rather than just flagging the\n",
     "marker as 'unstable' in the abstract.\n", sep = "")
}

## ---------------------------------------------------------------------------
## 3. PLOT
## ---------------------------------------------------------------------------
png(file.path(out_dir, "stability_plot.png"), width = 1000, height = 650, res = 120)
par(mar = c(7, 5, 3, 1))
ylim_max <- max(summary_tbl$max_abs_shift_upper) * 1.18
bp <- barplot(summary_tbl$max_abs_shift_upper, names.arg = summary_tbl$marker,
             las = 2, col = c("firebrick", "grey60")[1 + (summary_tbl$marker != "CD_batch")],
             ylim = c(0, ylim_max),
             ylab = "Max |threshold shift| across leave-one-out runs (Upper)",
             main = "Largest single-dataset threshold shift per marker")
text(bp, summary_tbl$max_abs_shift_upper, labels = round(summary_tbl$max_abs_shift_upper, 2),
    pos = 3, cex = 0.9)
dev.off()
cat("\nPlot written to", file.path(out_dir, "stability_plot.png"), "\n")

## ---------------------------------------------------------------------------
## ADAPT TO YOUR OWN DATA
## ---------------------------------------------------------------------------
## Replace section 1 above with:
##   data_norm <- readRDS("path/to/your_reference_NormData.rds")   # markers x cells
##   metadata  <- readRDS("path/to/your_reference_Metadata.rds")   # needs a 'dataset' column
## then rerun section 2 with sample_n/num.boot.sample/bootstrap_stratify/
## threshold_method matching whatever your real reference build uses (check
## the reference's *_provenance.json if available). With ~45 real datasets
## you will likely want ncores > 1.
