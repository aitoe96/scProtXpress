# scProtXpress

`scProtXpress` is a single-cell **protein** expression quantification tool with
its own identity: it aligns per-marker sample distributions to landmarks
(negative peak, valley, positive peak) and then **discretises** the resulting
matrix and emits per-population expression-level calls via a binomial test
on the discretised counts. The same pipeline handles **CyTOF** and **CITE-seq**
inputs through a single `run_scProteinXpress()` entry point that branches
only on the ingest step.

## Why a single pipeline

CyTOF and CITE-seq deliver the same biological readout (per-cell antibody
intensity) through different instruments and containers. Historically the
group maintained two pipelines that did the same downstream work twice
with slight numerical drift. `scProtXpress` collapses both into one
codepath, with two thin ingest helpers (`unpack_sce_cytof()`,
`unpack_list_citeseq()`) that hand the canonical matrix +
data.frame to the unified API.

## Method outline

1. **Reference build.** Per-marker landmark detection
   (`ADTnorm::get_peak_mode` patched, `get_valley_location`,
   `landmark_fill_na`, `peak_alignment`).
2. **Mode C global QC.** A single `cleanData_ADT()` pass over the
   NA-padded union matrix with marker-level coverage filtering and an
   exclusion audit CSV.
3. **Stratified bootstrap.** Per-marker threshold distributions
   sampled stratified by `cell.type [x tissue]`, with
   SLURM_CPUS_PER_TASK-aware ncores and BLAS pinning.
4. **Discretisation.** Each cell receives a {0, 0.5, 1, NA} label
   based on the ECDF of the reference bootstrap thresholds.
5. **Signature frame + binomial call.** `Create.sig.frame()` derives
   per-population baselines (`prob_low`, `prob_high`) on a coverage-
   filtered subset of markers, then `discrete.expr.levels()` calls
   each marker as `low` / `medium` / `high` / `not.significant` with
   two one-sided binomial tests.

## Installation

`scProtXpress` `Imports: ADTnorm`, which is currently a GitHub-only
package. Install ADTnorm first:

```r
install.packages("remotes")
remotes::install_github("yezhengSTAT/ADTnorm")
```

Then install `scProtXpress`:

```r
remotes::install_github("aitoe96/scProtXpress")
```

Optional dependencies (only required for specific paths):

| Use                                   | Package(s)                                |
| ------------------------------------- | ----------------------------------------- |
| Unpack a CyTOF SingleCellExperiment   | `SingleCellExperiment`, `SummarizedExperiment` |
| Cell-outlier filter (`qc_outlier_cells != "no"`) | `scuttle`                            |
| Faster BLAS pinning in workers        | `RhpcBLASctl`                             |
| SQS diagnostic plots                  | `ggplot2`                                 |

## Minimum CITE-seq use

```r
library(scProtXpress)
qobj     <- readRDS("data/citeseq_reference.rds")
unpacked <- unpack_list_citeseq(qobj)

run_scProteinXpress(
  data             = unpacked$data,
  metadata         = unpacked$metadata,
  dir              = "results",
  file.name        = "CITEseq_run",
  modality         = "CITE-seq",
  cell_id_col      = "cell.id",
  cell_type_col    = "cell.type_harmonised",
  num.boot.sample  = 1000L,
  ncores           = 4L,
  qc_outlier_cells = "auto"   # scuttle filter active for CITE-seq
)
```

## Minimum CyTOF use

```r
library(scProtXpress)
sce      <- readRDS("data/cytof_reference.rds")
unpacked <- unpack_sce_cytof(sce, assay_name = "exprs")

run_scProteinXpress(
  data            = unpacked$data,
  metadata        = unpacked$metadata,
  dir             = "results",
  file.name       = "CyTOF_run",
  modality        = "CyTOF",
  cell_id_col     = "cell_id",
  cell_type_col   = "cell_type",
  num.boot.sample = 1000L,
  ncores          = 4L,
  is_transformed  = unpacked$is_transformed   # TRUE for assay "exprs"
)
```

## Activating SQS (Stain Quality Score)

SQS is off by default. To turn it on for the reference build only:

```r
run_scProteinXpress(..., apply_sqs = "ref_only", sqs_threshold = 1.0)
```

Score logs land at `<Ref_dir>/<file.name>_SQS_ref_log.csv`. Optional
diagnostic PNG via the bundled reporter:

```r
source(system.file("sqs/sqs_report.R", package = "scProtXpress"))
df <- load_sqs_logs(ref_dir = "results/Ref_CyTOF_run")
save_sqs_report(df, out_dir = "results/SQS_report")
```

## Output files

```
Ref_<file.name>/
  <file.name>_NormData.rds          # normalised proteins x cells matrix
  <file.name>_Metadata.rds          # metadata for retained cells
  <file.name>_bg_prot_pars.rds      # per-marker fitted parameters
  <file.name>_th_dist.rds           # per-marker bootstrap thresholds
  <file.name>_sf.txt                # scaling factors (sf = 1 by default)
  <file.name>_BackCounts.rds        # sanitised reference counts (P1)
  <file.name>_BackMetadata.rds      # sanitised reference metadata (P1)
  <file.name>_excluded_markers.csv  # Mode C exclusion audit
  Query_<pop>/                      # per-cell-type discretisation
    <pop>_NormData.rds
    <pop>_pvals.rds
    <pop>_DiscretisedData.rds
  sig_frames/
    <pop>_sig_tbl.txt                       # binomial expression level
    audit_excluded_genes/
      <pop>_excluded_from_baseline.csv      # P3 coverage-filter audit
```

## SLURM template

A reference SLURM submission script lives at:

```
inst/slurm/scProtXpress.sh
```

It dispatches to either `inst/drivers/cytof_driver.R` or
`inst/drivers/citeseq_driver.R` based on its first argument:

```
sbatch inst/slurm/scProtXpress.sh CyTOF
sbatch inst/slurm/scProtXpress.sh CITE-seq
```

## Credits

`scProtXpress` builds on two GPL-licensed projects:

* **ADTnorm** (GPL-3) by Ye Zheng et al. ([yezhengSTAT/ADTnorm](https://github.com/yezhengSTAT/ADTnorm)).
  Imported as a regular dependency. `R/get_peak_mode_patched.R` is a
  local fork of `ADTnorm::get_peak_mode` that re-introduces
  `, drop = FALSE` at 12 internal subset operations; everything else
  delegates to upstream ADTnorm. The patch is enumerated in the file
  header and machine-readable as `inst/extdata/get_peak_mode.diff`
  (also shipped at `cytof_scProt/get_peak_mode.diff` in the source
  history).
* **scGeneXpress** (GPL-2) by Sascha Jung et al.
  ([saschajung/scGeneXpress](https://github.com/saschajung/scGeneXpress)).
  No code from scGeneXpress is vendored; its discretisation +
  binomial-test methodology informs the design but the implementation
  here is independent.

This is **not** "an adaptation of scGeneXpress to proteomics". It is a
purpose-built proteomic-discretisation pipeline that draws on the same
statistical idea.

## Licence

`scProtXpress` is distributed under GPL-3 (see `LICENSE.md`). GPL-3 is
required because the package vendors a single function derived from
GPL-3-licensed ADTnorm (`R/get_peak_mode_patched.R`).
