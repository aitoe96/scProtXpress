# Changelog

All notable changes to scProtXpress are documented here.

## [0.3.1]

### Added

- **FDR for the population-level calls (`sig.frames`).** New argument
  `fdr_method` of `run_scProteinXpress()` (default `"BH"`), also in
  `discrete.expr.levels_donor()` (default `NULL`); the pooled
  `discrete.expr.levels()` already had it. The correction is done within
  each population, over its gene panel, in two stages: `p_low` is
  Benjamini-Hochberg adjusted over every gene (`q_low`); then `p_high` over
  the genes that reach the medium/high test only, i.e. not low after `q_low`
  and not tied (`q_high`). `expr.level` is called on the q-values at 0.05.
  The raw `p_low`/`p_high` and the call on them (`expr.level_raw`) are kept
  in the table. `fdr_method = NULL` gives the calls of 0.3.0.
- **Positive peak of the background.** When no sample of the background
  detects the positive peak of a marker, it is imputed in every sample with
  the median of the per-sample maxima (all cells of the sample, before
  alignment); before, such a marker was aligned with two landmarks and had
  no positive target. A positive peak at or below the valley of its sample
  is replaced by the maximum of that sample. Markers where at least one
  sample detects it are imputed as before (`ADTnorm::landmark_fill_na()`).
- **Single-peak query samples.** No positive peak is imputed in a query. A
  sample with both peaks is aligned with its three landmarks; a sample with
  a single peak with that peak and its valley onto the matching background
  targets. New argument `single_peak` of `run_scProteinXpress()`,
  `PrepQueryData_ADT()` and `ADTprot_norm()`: `"nearest"` (default) takes
  the peak as the positive one when it is nearer to the positive target of
  the background than to the negative one (valley searched below the peak,
  as ADTnorm does for `positive_peak`; warp onto the valley and positive
  targets), and as the negative one otherwise (warp onto the negative and
  valley targets); `"negative"` always takes it as the negative one.
- **Cell-type balancing of the queries is optional.** New argument
  `balance_celltypes_query` (default `FALSE`) of the same functions: query
  samples are not balanced unless it is `TRUE`, which uses the balancing
  settings saved in the background. The background is still balanced by
  default (`balance_celltypes = TRUE`), and queries still use its bimodal
  markers and `min_cells_align`.
- **Debug files.** New argument `debug` (default `FALSE`) of the same
  functions. For the background and the query it writes
  `<file.id>_debug_landmarks.csv` (one row per marker and sample: values,
  status of the sample and why, detected peaks and valley, landmarks used by
  the warp, which were imputed and how, targets), `<file.id>_debug/<marker>.rds`
  (detection and alignment objects of the marker; the cells of the balanced
  copy and how many times each was drawn, without cell values) and
  `<file.id>_debug_settings.csv`.

### Changed

- `discrete.expr.levels_donor()` computes `p_high` for every gene with
  non-zero cells and no 0.5/1 tie (before: only when `p_low` was
  significant), so that both the raw and the FDR calls can be made; the
  calls with `fdr_method = NULL` are unchanged.
- In the pooled `discrete.expr.levels()` with `fdr_method = "BH"`, `p_high`
  is now adjusted over the genes that reach the medium/high test only
  (before: over every gene with non-zero cells).
- The bimodal peak detector stays on by default (`bimodal_marker = "auto"`,
  as in 0.3.0).
- Version 0.3.1.

## [0.3.0]

### Added

**Cell-type balanced landmark detection (on by default).** The peaks and
valleys of every sample are now detected on a copy of the sample in which
its cell types are balanced, and the landmarks found there warp the
original cells (the copy is discarded). Each sample's copy has `balance_n`
cells (default 2000) shared equally among the cell types with at least
`balance_min_cells` cells (default 1): a smaller type keeps all its cells
and is completed with draws with replacement, a larger one is down-sampled
without replacement (`balance_replace = TRUE`: with replacement). Samples
with fewer than two eligible types are used as they are. `balance_n` also
accepts `"sample"` (as many cells as the sample) and `"maxK"` (every type as
large as the largest). New arguments of `run_scProteinXpress()` and
`ADTprot_norm()`: `balance_celltypes`, `balance_n`, `balance_min_cells`,
`balance_replace`. The defaults are the configuration chosen in the
evaluation of the method (2000-5000 cells per sample performed best; 2000
is the cheaper).

**Bimodal peak detector on by default.** `bimodal_marker` now defaults to
`"auto"`: every marker when `balance_celltypes = TRUE`, none otherwise
(`"all"`, `NULL` or marker names are also accepted). It is resolved to
marker names before it is saved in the reference.

**Thresholds per group of cell types (optional).** New arguments of
`run_scProteinXpress()`: `background_groups` and `query_groups` (data.frames
with columns `group` and `label`, labels of `cell_type_col`; a label can be
in several groups), `min_cells_group` (default 200) and `min_donors_group`
(default 2).

- **Background.** With `background_groups`, after the global thresholds one
  bootstrap per group is run on the aligned background cells of its labels,
  with the same settings. A group is built when the background has more
  than `min_cells_group` of its cells, from at least `min_donors_group`
  donors (dataset::donor). The results are saved in
  `<file.name>_th_dist_groups.rds`, `<file.name>_groups.csv` (cells, donors,
  built or why not) and `groups/<group>/`.
- **Queries.** With `query_groups` and `precBack`, each population is also
  discretised with the thresholds of its groups, in
  `groups/<group>/Query_<label>/` and `groups/<group>/sig_frames/`; the
  self-query of a background uses `background_groups`. Populations without a
  group, or whose groups were not built, keep only the global outputs and
  are listed in `groups_unmatched.csv` and in a warning.
- **What does not change.** The alignment and the global outputs are the
  same as without groups. A background built before without groups is
  reused, and only the outputs per group are added.

Without these arguments nothing changes.

**Alignment log.** Samples set to NA for a marker (`too_small`,
`not_valid`, `warp_failed`) or whose valley was imputed
(`no_valley_imputed`) are reported in one message per marker and listed in
`<file.id>_alignment_log.csv` (reference and query).

### Changed

- **Queries follow their background.** Balancing, bimodal markers and
  `min_cells_align` are saved in `<file.name>_bg_prot_pars.rds`
  (`global_settings`, with `alignment_version = "0.3.0"`), and queries always
  use them; explicit arguments of `run_scProteinXpress()` that differ are
  ignored with a warning. A background built with scProtXpress <= 0.2.0 is
  followed as it was built (no balancing, valley call of that version), with
  a message recommending to rebuild it.
- **Minimum cells per sample.** A sample is aligned for a marker only with
  at least `min_cells_align` finite values (default 25; was 5), >= 3
  distinct values and non-zero variance. Samples that are not aligned are
  set to NA for that marker; before, they kept their unaligned values.
- **Seeds.** The landmark detection of every marker sets `set.seed(seed)`
  (balanced copy and the jitter of the peak detector), in reference and
  query; `seed` is `fixseed` of `run_scProteinXpress()` and a new argument of
  `ADTprot_norm()` and `PrepQueryData_ADT()`. Results no longer depend on the
  number of workers. The reference subsampling also uses this seed (it was
  fixed at 1234, the default of `fixseed`).
- **Collation.** `run_scProteinXpress()`, `ADTprot_norm()`,
  `PrepQueryData_ADT()`, `GetBackground()` and `RunBootstrap()` run in C
  collation (`LC_COLLATE`) and restore the caller's on exit. The reference
  subsampling and the stratified bootstrap order cell types, tissues and
  datasets before drawing cells, so with the same seed they drew different
  cells in different locales (e.g. `C` vs `en_US.UTF-8`); the result is now
  the same on any machine.
- Cell-type labels with "/" or "\\" are written with "_" in file and folder
  names (`Query_<label>/...`); the labels are not changed.
- `num.cells.sample = NULL` now resolves to 1000L for CITE-seq too (was
  100L).
- A background can be built with a single cell type (was >= 2).
- Version 0.3.0.

### Fixed

- **Bootstrap replicates (bug C).** `RunBootstrap()` fork workers all
  started from the parent's RNG state, so most iterations repeated the same
  draws (about 35 distinct of 1000). Every iteration now sets its own seed,
  drawn from `fixseed`; the result is the same with any number of cores.
  This also affects `GetBackground()`, `leave_one_dataset_out_stability()`,
  `downsample_reference_stability()` and `add_reference_dataset()`.
- **Warp failures (bug B).** When `fda::landmarkreg` fails, every sample is
  warped on its own with the settings of the whole call; only the samples
  that still fail are set to NA (`warp_failed`), instead of losing the whole
  marker. If every sample fails the marker is not aligned.
- **Missing valley.** `ADTnorm::get_valley_location()` stopped the whole
  marker when a single-peak sample had neither valley nor shoulder. A local
  copy, `.get_valley_location_patched()`, leaves that valley NA, imputed by
  `ADTnorm::landmark_fill_na()`.
- **Valley call arguments (bug D).** `get_valley_location()` was called by
  position and skipped `lower_peak_thres`, so `neg_candidate_thres` was
  used as `lower_peak_thres`, the arcsinh flag as `neg_candidate_thres`
  (1 or 0) and `arcsinh_transform_flag` was always TRUE. The arguments are
  now passed by name, as ADTnorm does: the valley floor of single-peak
  samples with > 80 % negative cells is `neg_candidate_thres` again.
- **Bimodal / trimodal markers in queries.** Queries looked up
  `bimodal_marker_index` / `trimodal_marker_index`, which the reference
  never saved, so the bimodal and trimodal detectors were never applied to
  queries. They are now taken from the marker names saved in the reference.
- **Alignment without the ADTnorm fork.** `.peak_alignment_scprot()`, an
  internal copy of the fixed-mode `peak_alignment()` of the scProt fork of
  ADTnorm, is used for reference and query alignment (published ADTnorm
  left queries unaligned without error); `fda` added to `Imports`.
- **Samples without peaks (bug A).** `.get_peak_mode_patched()` records a
  sample with no peak as NA, imputed by `landmark_fill_na()`, instead of
  failing the whole marker.
- **Markers that cannot be aligned** are set to NA (reference: left out of
  the background; query: not discretised) and reported with the failing
  step, also from parallel workers.

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
