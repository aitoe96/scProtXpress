# Changelog

All notable changes to scProtXpress are documented here.

## [0.2.0]

### Added

**Mixture-model threshold estimation.** New `threshold_method` option on
`run_scProteinXpress()`/`GetBackground()`: `"mixture"` fits a 1-3 component
Gaussian mixture (via `mclust`) to each bootstrap draw and sets the
low/high cut at the Bayes decision boundary between fitted components,
rather than the `"rectangle"` heuristic's cumulative-distribution
approximation. Falls back to `"rectangle"` automatically when `mclust`
isn't installed, the fit is unstable, or BIC selects a single component.
On synthetic data with known ground-truth crossing points, `"mixture"`
recovers the true boundary more accurately than `"rectangle"` in
well-separated, close, and standard trimodal scenarios, but can
under-perform when a true component is small and heavily overlapping
(BIC-based selection tends to under-fit rare components in that regime).
Two additional simple reference-informed baselines are also available:
`"quantile"` (fixed percentile cutoffs) and `"zscore"` (mean +/- k*SD).
See `inst/simulations/threshold_method_comparison.R` for a full comparison.

**Donor-aware population-level testing.** The population-level
significance test previously pooled every cell in a population regardless
of donor, which inflates the false-positive rate whenever cells from the
same donor are correlated (a realistic and common condition, not an edge
case). New `sig_test_method` parameter on `run_scProteinXpress()`:
`"pooled"` (original behaviour, default, kept for backward compatibility),
`"donor_t"`, `"donor_wilcox"`, or `"design_effect"` to test across donors
instead. New `donor_col` parameter to specify the donor identifier
explicitly when it differs from the sample identifier (e.g. longitudinal
designs with multiple timepoints per donor, multi-tissue-per-donor
designs, or multiplexed/hashed runs). New `min_donors_for_call` parameter
(default 2) forces populations with too few contributing donors to
`"not.significant"` rather than a potentially spurious confident call.
Underlying new functions: `Create.sig.frame_by_donor()`,
`discrete.expr.levels_donor()`, `.estimate_icc()`. On simulated data with
realistic between-donor heterogeneity, the pooled test's false-positive
rate inflates several-fold above nominal as donor clustering increases,
while the donor-aware tests stay close to nominal across the same range.
See `inst/simulations/donor_aware_testing.R`.

**Reference robustness and validation tooling**, in `R/validation_tools.R`
and `R/internalFunctions.R`:
- `leave_one_dataset_out_stability()`: rebuilds a reference's bootstrap
  thresholds with each contributing dataset held out in turn, and reports
  how much each marker's threshold shifts, to identify markers
  disproportionately determined by a single dataset.
- `downsample_reference_stability()`: rebuilds thresholds at a sweep of
  reference-size fractions (by cells or by whole datasets), with multiple
  replicates per fraction, to quantify how threshold estimates converge as
  reference size grows.
- `tag_marker_specificity()`: classifies markers in a long-format
  marker-to-cell-type map as "broad" (canonical for more than one cell
  type) or "specific" (canonical for exactly one).
- `evaluate_calls()`: computes precision, recall, F1, and false-positive
  rate for a table of predicted calls against a ground-truth set,
  optionally stratified by a grouping column (e.g. marker specificity
  tier).
- `fit_patient_classifier()`: cross-validated (leave-one-out or k-fold)
  logistic regression on per-patient summary features, reporting AUC,
  sensitivity, and specificity, explicitly framed as a proof-of-concept
  tool rather than a validated diagnostic. Includes an `n_permutations`
  option that builds an empirical null distribution by re-running the
  entire cross-validation procedure with shuffled labels; this exists
  because cross-validated AUC on small, weakly- or non-predictive cohorts
  is not centered on 0.5 (the bias worsens as the method approaches
  leave-one-out), so a raw AUC should not be compared to 0.5 directly for
  small cohorts -- the permutation p-value is the statistically
  appropriate comparison.

**Incremental reference updates.** New `add_reference_dataset()` in
`R/reference_update.R`: given an existing reference and a new dataset,
splits the new dataset's markers into already-known (query-aligned to
existing landmark parameters) and genuinely new (reference-aligned from
scratch), then re-bootstraps only the markers the new dataset touched,
leaving every other marker's threshold distribution unaffected. Untouched
markers are verified `identical()` before and after; touched markers
converge to the same estimate as a full rebuild within ordinary bootstrap
sampling noise. Reads prior build settings from the provenance sidecar
(below) when not given explicitly, backs up existing reference files
before writing, and writes a per-marker shift audit CSV. See
`inst/simulations/incremental_reference_update.R`.

**Reference provenance tracking.** `GetBackground()` and
`add_reference_dataset()` now call `.write_provenance()`, which writes a
`<file.id>_provenance.json` (or `.txt` if `jsonlite` isn't installed)
alongside every reference's outputs, recording the package version, git
commit (best effort), and every parameter that determined that reference
(`bootstrap_stratify`, `threshold_method`, `num.boot.sample`, `sample_n`,
`fixseed`), so the effective build settings are self-documenting.

**External population baseline.** New `ComputeExternalBaseline()` and
`apply_isotype_floor()` in `R/internalFunctions.R`: computes the
population-level baseline probabilities (`prob_low`/`prob_high`) from a
caller-supplied independent discretised matrix rather than the same
population being tested, avoiding a self-referential baseline.
`apply_isotype_floor()` optionally floors `prob_low` using isotype/
negative-control ADTs' own non-zero rate when the panel includes them.

**Multiple-testing correction options.** `discrete.expr.levels()` gained
`fdr_method = "BH"` (Benjamini-Hochberg correction within a population's
whole-marker-panel scan), `min_n_for_call` (forces `"not.significant"`
below a minimum cell count rather than only on an exact tie), and
`effect_size` (adds `effect_size_low`/`effect_size_high` columns). Raw
`p_low`/`p_high` are now always returned as columns. All are additive and
off-by-default or unchanged-by-default; default-argument output is
byte-for-byte identical to the previous behaviour on shared columns. New
`correct_across_scan()` applies BH correction to `p_low`/`p_high` jointly
across an entire list of per-population tables (e.g. a full disease-panel
scan spanning many proteins, cell types, and diseases), which
within-table correction does not do. `Discretise_data()`'s previously
hardcoded per-cell correction behaviour is now an explicit `p_correct`
argument (default unchanged, `FALSE`).

**Input-scale sanity check.** New `.sanity_check_transformed_scale()`
(used by `.detect_input_mode()`) warns when input declared or detected as
already-transformed has a substantial fraction of negative values (default
trigger: more than 5%), which is inconsistent with arcsinh-transformed
non-negative counts and more consistent with DSB- or CLR-normalised input
being scored on the wrong scale.

**SQS before/after comparison script.** New
`inst/analysis/sqs_before_after_comparison.R` with `compare_sqs_runs()`:
builds a reference twice (with and without Stain Quality Score filtering)
and quantifies how many (sample, marker) combinations were excluded and
how much each affected marker's threshold moved as a result.

**Validation and simulation scripts**, in `inst/simulations/`: each is
self-contained, simulates illustrative data with known ground truth, runs
the relevant package function(s), and reports a summary table, CSV, and
plot: `reference_dataset_stability.R`, `marker_specificity_benchmark.R`,
`precision_recall_negative_markers.R`, `threshold_method_comparison.R`,
`donor_aware_testing.R`, `incremental_reference_update.R`,
`reference_size_convergence.R`, `patient_level_classifier.R`.

**Package metadata.** `Authors@R` now lists all four manuscript authors.
`inst/CITATION` includes both a manuscript entry and a software entry with
correct author names. `mclust` and `jsonlite` added to `Suggests` (used by
the mixture threshold estimator and the provenance sidecar respectively,
both optional with graceful fallback if absent). `NAMESPACE` is now
generated by roxygen2 from source annotations rather than hand-maintained,
which also produced previously-missing `man/*.Rd` documentation files.

### Changed

- `inst/drivers/citeseq_driver.R` and `inst/drivers/cytof_driver.R` now
  explicitly set `bootstrap_stratify`, `apply_sqs`, and `threshold_method`
  (env-var overridable, defaulting to `"cell.type_tissue"`, `"ref_only"`,
  and `"mixture"` respectively) instead of leaving them unset and
  silently resolving to `"none"`/`"none"`/`"rectangle"`. A `[CONFIG]`
  message prints the resolved values at the start of every run.
- `run_scProteinXpress()`'s `sig.frames` step now dispatches between
  pooled and donor-aware population testing via `sig_test_method` and a
  new internal `.compute_sig_table_for_population()`, rather than always
  running the original pooled-cell test unconditionally. Each
  population's donor/sample metadata is now persisted alongside its
  discretised data so this can be used even when discretization and
  significance testing happen in separate calls.

### Fixed

- `downsample_reference_stability()`: replicate draws within the same
  subsample fraction were not statistically independent, because
  `GetBackground()` resets the RNG to a `fixseed`-determined state on
  every call, and successive replicates in the loop called it repeatedly
  without an intervening seed change. Replicates now use a distinct,
  deterministic seed per (fraction, replicate) combination, while the
  function remains exactly reproducible given the same `fixseed`.
- `run_scProteinXpress()` returned a visibly-printed `NULL` at the end of
  the function instead of `invisible(NULL)` as documented (`return(x)` is
  always visible in R regardless of `x`; only `invisible(x)` suppresses
  auto-printing).
- `NAMESPACE` was missing an `importFrom(magrittr, "%>%")` entry despite
  `get_peak_mode_patched.R` using `%>%` in several places; resolved by
  regenerating `NAMESPACE` with roxygen2 and fully qualifying the
  remaining base/stats calls that don't need an import
  (`stats::dnorm`, `stats::uniroot`, `stats::t.test`,
  `stats::wilcox.test`).
- `Create.sig.frame_by_donor()`/`discrete.expr.levels_donor()` existed as
  tested, standalone functions but were not called anywhere in
  `run_scProteinXpress()`'s pipeline; `sig.frames = TRUE` always ran the
  original pooled-cell test regardless of intent. Now wired in via
  `sig_test_method` (see Changed, above).
- `sqs_before_after_comparison.R`'s SQS log reader: the log's `action`
  column uses the literal string `"NA"` to mean "excluded", which
  `read.csv()`'s default `na.strings` silently converted to a real `NA`,
  breaking any filter on `action == "NA"`. Fixed by reading with
  `na.strings = character(0)`.
- A flaky threshold-estimator test asserted a bound on a single stochastic
  draw, which occasionally failed depending on execution context. Now
  averages over 8 independent draws with an appropriately wider bound.

### Known limitations

- `add_reference_dataset()`'s and `compare_sqs_runs()`'s ADTnorm-dependent
  alignment calls are exercised by the package's existing smoke tests and
  used exactly as `run_scProteinXpress()` already uses them, but have not
  been independently run end-to-end against a real CITE-seq/CyTOF dataset
  outside that shared path. Validate the full pipeline against your own
  data before relying on either in production.
- `ComputeExternalBaseline()`/`apply_isotype_floor()` are validated on
  synthetic data (return well-formed probabilities, floor behaves
  monotonically); validating their effect on real disease-cohort biomarker
  calls requires real data and is left to the user.
