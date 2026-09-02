#' From raw counts to a donor-aware disease-vs-healthy sig.frame
#'
#' Full pipeline example: QC + ADTnorm normalization + background
#' calculation + discretization + population-level significance testing,
#' with the population-level test done ACROSS DONORS rather than pooled
#' cells.
#'
#' Two things this specifically demonstrates, both easy to get wrong
#' silently:
#'   1. sig_test_method must be set explicitly -- the default ("pooled")
#'      reproduces the original pseudoreplication-prone behaviour exactly,
#'      by design, for backward compatibility.
#'   2. donor_col must be set explicitly whenever "one sample" != "one
#'      donor" for your data. This example uses a deliberately realistic
#'      case: 3 donors, each contributing 2 samples (e.g. two timepoints,
#'      or two tissues) -- if you leave donor_col unset here, the pipeline
#'      silently treats each of the 6 samples as its own "donor", which
#'      reintroduces pseudoreplication one level up from the problem
#'      sig_test_method exists to fix. Verified in this exact scenario:
#'      unset -> resolves to (wrongly) 6 "donors"; donor_col set -> 3
#'      (right).

suppressPackageStartupMessages(library(scProtXpress))

# ----------------------------------------------------------------------------
# 1. LOAD THE LIST (raw counts, same shape as before)
# ----------------------------------------------------------------------------
raw_data <- readRDS("path/to/your_raw_counts_list.rds")
unpacked <- unpack_list_citeseq(raw_data)

out_dir <- "path/to/your_output_dir"   # <- EDIT

# ----------------------------------------------------------------------------
# 2. QC + NORMALIZE + BACKGROUND, WITH DISCRETIZATION + SIG.FRAMES TURNED ON
# ----------------------------------------------------------------------------
run_scProteinXpress(
  data     = unpacked$data,
  metadata = unpacked$metadata,
  dir      = out_dir,
  file.name = "my_reference",
  modality  = "CITE-seq",

  cell_id_col   = "cell_id",     # <- EDIT to your columns
  cell_type_col = "cell_type",
  dataset_col   = "dataset_id",
  sample_col    = "sample_id",   # one sample per timepoint/tissue in this example
  batch_col     = "batch",
  tissue_col    = "tissue_id",

  # THE column that actually identifies the donor/patient -- NOT the same
  # as sample_col here (this example has 2 samples per donor). If you
  # don't set this and donor really does equal sample for your data,
  # leave it out entirely -- it defaults to sample_col's resolution, which
  # is exactly right in that case and needs no extra configuration.
  donor_col = "patient_ref",     # <- EDIT (or delete this line if donor == sample for you)

  is_transformed = unpacked$is_transformed,
  qc_outlier_cells = "no",       # <- set "auto" if you have the scuttle package installed

  discretize = TRUE,
  sig.frames = TRUE,

  num.boot.sample  = 1000L,
  num.cells.sample = 100L,
  ncores           = 4L,         # <- EDIT
  bootstrap_stratify = "cell.type_tissue",   # stratifies the REFERENCE (Rev. 1 Major 1)
  threshold_method   = "mixture",

  # stratifies the DISEASE-TEST (Rev. 1 Major 5) -- a completely separate
  # stage and a completely separate fix from bootstrap_stratify above.
  # "pooled" (the default if you omit this) reproduces the original,
  # pseudoreplication-prone behaviour -- you have to opt into the fix.
  sig_test_method      = "design_effect",
  min_donors_for_call  = 2L
)

# ----------------------------------------------------------------------------
# 3. RESULT
# ----------------------------------------------------------------------------
# out_dir/Ref_my_reference/sig_frames/<celltype>_sig_tbl.txt, one per
# population, now computed across donors (patient_ref) rather than pooled
# cells. Columns match discrete.expr.levels_donor()'s output: gene, n.0,
# n.0.5, n.1, n_donors, p_low, p_high, expr.level.
sig_dir <- file.path(out_dir, "Ref_my_reference", "sig_frames")
print(list.files(sig_dir, pattern = "_sig_tbl.txt$"))

example_tbl <- read.delim(list.files(sig_dir, pattern = "_sig_tbl.txt$", full.names = TRUE)[1])
head(example_tbl)

# ----------------------------------------------------------------------------
# WHAT HAPPENS IF YOU GET donor_col WRONG (for reference -- not run here)
# ----------------------------------------------------------------------------
# Leaving donor_col unset when donor != sample doesn't error -- it silently
# resolves "donor" to whatever sample_col gave you, so every population's
# sig-table would be tested across the WRONG grouping (6 "donors" instead
# of the true 3 in this example). There is no warning for this specific
# mistake because the pipeline has no way to know your samples and donors
# differ -- only you know that from how the data was collected. If you're
# not sure whether this applies to you, ask: "can the same person/animal
# contribute more than one row of my sample_id column?" If yes, set
# donor_col to whatever column tracks that identity instead.
