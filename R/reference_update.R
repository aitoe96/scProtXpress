###############################################################################
# reference_update.R -- incremental reference expansion
#
# run_scProteinXpress()'s default bootstrap ("none") pools every cell's
# normalised value per marker into one vector before resampling from it
# (see RunBootstrap()'s docstring). That means adding a new dataset to the
# reference does NOT require rebuilding peak alignment for markers the new
# dataset doesn't touch, and does not require re-bootstrapping markers whose
# pooled values haven't changed. add_reference_dataset() exploits exactly
# that: it aligns the new dataset's cells (reusing the *existing* per-marker
# landmark parameters for markers already in the reference, via
# PrepQueryData_ADT -- precisely what happens to any query), fits fresh
# landmark parameters only for genuinely new markers (via ADTprot_norm in
# "Ref" mode, restricted to those markers), re-bootstraps ONLY the markers
# touched by the new dataset, and leaves every other marker's threshold
# distribution byte-for-byte untouched.
#
# The bookkeeping in this file (which markers count as "touched", the merge
# of old/new threshold lists, the audit report) is unit-tested in
# tests/testthat/test-add-reference-dataset.R against mocked threshold-list
# objects and does not require ADTnorm. The ADTnorm-dependent alignment
# calls (PrepQueryData_ADT / ADTprot_norm) are exercised by the package's
# existing smoke tests and are used here exactly as run_scProteinXpress()
# already uses them. Validate the full incremental-update path against your
# own reference before relying on it in production.
###############################################################################

#' Incrementally add a new dataset to an existing scProtXpress reference
#'
#' @param ref_dir Path to the existing `Ref_<file.name>` directory (as
#'   produced by `run_scProteinXpress()` / `GetBackground()`).
#' @param file.name The `file.name` used to build that reference (used to
#'   locate `<file.name>_NormData.rds`, `_Metadata.rds`, `_bg_prot_pars.rds`,
#'   `_th_dist.rds`, and, if present, `_provenance.json`/`.txt`).
#' @param new_data,new_metadata The new dataset only (not concatenated with
#'   the old data) -- proteins x cells matrix/Matrix and a metadata
#'   data.frame, same contract as `run_scProteinXpress(data=, metadata=)`.
#' @param modality "CyTOF" or "CITE-seq".
#' @param cell_id_col,cell_type_col,dataset_col,sample_col,batch_col,tissue_col
#'   As in `run_scProteinXpress()`.
#' @param is_transformed As in `run_scProteinXpress()`.
#' @param out_dir Where to write the updated reference. Default: overwrite
#'   `ref_dir` in place (after backing up; see `backup`).
#' @param backup logical, default `TRUE`. Copies the four existing `.rds`
#'   files to `<ref_dir>/backup_<timestamp>/` before writing anything.
#' @param num.boot.sample,num.cells.sample,ncores,fixseed,bootstrap_stratify,threshold_method
#'   Bootstrap parameters for the *touched-marker* re-run. If a
#'   `_provenance.json`/`.txt` sidecar exists for the old reference and an
#'   argument here is left `NULL`, its value is read from the sidecar so the
#'   update uses the same settings as the original build unless you
#'   deliberately override them (a mismatch here -- e.g. rectangle vs.
#'   mixture thresholds between the old and new marker pools -- would
#'   silently make the reference internally inconsistent).
#' @param min.cells.feature,min.coverage.fraction,qc_outlier_cells,qc_nmads
#'   QC parameters applied to `new_data` (same semantics as
#'   `run_scProteinXpress()`); ideally should match what was used to build
#'   `ref_dir` originally.
#' @param sqs_threshold,apply_sqs SQS parameters for aligning `new_data`.
#' @return Invisibly, a list with `touched_markers`, `new_markers`,
#'   `unchanged_markers`, and `threshold_shift` (a data.frame of per-touched-
#'   marker summary shift in the bootstrap Lower/Upper means, old vs new).
#'   Also written to `<out_dir>/<file.name>_update_audit.csv`.
#' @export
add_reference_dataset <- function(
    ref_dir, file.name, new_data, new_metadata,
    modality = c("CyTOF", "CITE-seq"),
    cell_id_col   = "cell_id",
    cell_type_col = "cell_type",
    dataset_col   = "dataset_id",
    sample_col    = "sample_id",
    batch_col     = "batch",
    tissue_col    = "tissue",
    tissue_default = "PBMC",
    is_transformed = NULL,
    out_dir = ref_dir,
    backup = TRUE,
    num.boot.sample = NULL, num.cells.sample = NULL, ncores = 1L, fixseed = NULL,
    bootstrap_stratify = NULL, threshold_method = NULL,
    min.cells.feature = 10, min.coverage.fraction = 0.10,
    qc_outlier_cells = c("auto", "yes", "no"), qc_nmads = 3,
    sqs_threshold = 1.0, apply_sqs = c("none", "ref_only", "both")
){
  modality         <- match.arg(modality)
  qc_outlier_cells <- match.arg(qc_outlier_cells)
  apply_sqs        <- match.arg(apply_sqs)

  ## ---- 0. locate + load existing reference artefacts ---------------------
  p_norm  <- file.path(ref_dir, paste0(file.name, "_NormData.rds"))
  p_meta  <- file.path(ref_dir, paste0(file.name, "_Metadata.rds"))
  p_pars  <- file.path(ref_dir, paste0(file.name, "_bg_prot_pars.rds"))
  p_th    <- file.path(ref_dir, paste0(file.name, "_th_dist.rds"))
  for (p in c(p_norm, p_meta, p_pars, p_th)) {
    if (!file.exists(p)) stop("add_reference_dataset: expected file not found: ", p)
  }
  old_norm <- readRDS(p_norm)
  old_meta <- readRDS(p_meta)
  bg.prot.pars <- readRDS(p_pars)
  old_th   <- readRDS(p_th)

  ## ---- resolve bootstrap params from the provenance sidecar when left NULL
  prov <- .read_provenance_if_exists(ref_dir, file.name)
  `%||%` <- function(a, b) if (is.null(a)) b else a
  num.boot.sample    <- num.boot.sample    %||% prov$num.boot.sample    %||% 1000L
  num.cells.sample   <- num.cells.sample   %||% prov$sample_n           %||% (if (identical(modality, "CyTOF")) 1000L else 100L)
  fixseed            <- fixseed            %||% prov$fixseed           %||% 1234L
  bootstrap_stratify <- bootstrap_stratify %||% prov$bootstrap_stratify %||% "none"
  threshold_method   <- threshold_method   %||% prov$threshold_method  %||% "rectangle"
  if (is.null(prov)) {
    warning("add_reference_dataset: no _provenance.json/.txt found next to ", ref_dir,
            " -- falling back to defaults (bootstrap_stratify='", bootstrap_stratify,
            "', threshold_method='", threshold_method, "'). If this reference predates ",
            "the provenance sidecar, double-check these match what was actually used ",
            "to build it before trusting the merged result.", call. = FALSE)
  }

  ## ---- 1. QC + resolve metadata columns on the NEW data only -------------
  new_meta_resolved <- .resolve_new_metadata_cols(
    new_metadata, cell_id_col, cell_type_col, dataset_col, sample_col, batch_col,
    tissue_col, tissue_default
  )
  qc_outlier_resolved <- if (identical(qc_outlier_cells, "auto")) {
    if (identical(modality, "CITE-seq")) "yes" else "no"
  } else qc_outlier_cells

  qc <- cleanData_ADT(new_data, nmads = qc_nmads,
                      min.cells.feature = min.cells.feature,
                      min.coverage.fraction = min.coverage.fraction,
                      qc_outlier_cells = qc_outlier_resolved)
  new_counts_qc <- qc$Counts
  new_meta_qc   <- new_meta_resolved[new_meta_resolved$cell.id %in% colnames(new_counts_qc), , drop = FALSE]
  new_meta_qc   <- new_meta_qc[match(colnames(new_counts_qc), new_meta_qc$cell.id), , drop = FALSE]

  ## ---- 2. split markers: already-known (query-align) vs. brand-new (ref-align)
  new_markers      <- rownames(new_counts_qc)
  existing_markers <- names(bg.prot.pars)
  markers_known    <- intersect(new_markers, existing_markers)
  markers_brand_new <- setdiff(new_markers, existing_markers)

  message(sprintf(
    "[add_reference_dataset] New dataset contributes %d markers: %d already in the reference (will be QUERY-aligned to existing landmarks), %d genuinely new (will be REF-aligned from scratch on the new data alone).",
    length(new_markers), length(markers_known), length(markers_brand_new)))

  tmp_dir <- file.path(tempdir(), paste0("add_ref_", file.name, "_", as.integer(Sys.time())))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  aligned_known <- NULL
  if (length(markers_known) > 0) {
    aligned_known <- PrepQueryData_ADT(
      query_dat      = new_counts_qc[markers_known, , drop = FALSE],
      query_metadata = new_meta_qc,
      Ref_data       = ref_dir,
      dir            = tmp_dir,
      file.id        = "incoming_known",
      bg.prot.pars   = bg.prot.pars[markers_known],
      normalise      = TRUE, scale_data = TRUE, saveData = FALSE, return.data = TRUE,
      is_transformed = is_transformed,
      sqs_threshold  = sqs_threshold, apply_sqs = apply_sqs,
      min.cells.feature = 1, min.coverage.fraction = 0,  # already QC'd above; don't filter twice
      qc_outlier_cells  = "no",                          # already applied above
      ncores = ncores
    )
  }

  aligned_new <- NULL
  new_pars    <- NULL
  if (length(markers_brand_new) > 0) {
    ref_res <- ADTprot_norm(
      data_counts = new_counts_qc[markers_brand_new, , drop = FALSE],
      metadata    = new_meta_qc,
      data.type   = "Ref",
      file.id     = "incoming_new_markers",
      saveData    = FALSE, return = TRUE, outdir = tmp_dir,
      sample_col = "sample", batch_col = "batch",
      marker_to_process = markers_brand_new,
      is_transformed = is_transformed,
      sqs_threshold = sqs_threshold, apply_sqs = apply_sqs
    )
    aligned_new <- .extract_align_field(ref_res, "NormData", pos = 1)
    new_pars    <- .extract_align_field(ref_res, "bg.prot.pars", pos = 2)
  }

  ## ---- 3. assemble the updated NormData / Metadata / bg.prot.pars --------
  all_markers <- union(rownames(old_norm), new_markers)
  pad <- function(mat, markers, cells) {
    out <- matrix(NA_real_, nrow = length(markers), ncol = length(cells),
                  dimnames = list(markers, cells))
    common_m <- intersect(rownames(mat), markers)
    out[common_m, colnames(mat)] <- as.matrix(mat[common_m, , drop = FALSE])
    out
  }
  new_aligned_combined <- if (!is.null(aligned_known) && !is.null(aligned_new)) {
    rbind(aligned_known, aligned_new)[, colnames(new_counts_qc), drop = FALSE]
  } else if (!is.null(aligned_known)) aligned_known else aligned_new

  old_padded <- pad(old_norm, all_markers, colnames(old_norm))
  new_padded <- pad(new_aligned_combined, all_markers, colnames(new_aligned_combined))
  updated_norm <- cbind(old_padded, new_padded)

  updated_meta <- rbind(old_meta, new_meta_qc[, colnames(old_meta), drop = FALSE])

  updated_pars <- bg.prot.pars
  if (!is.null(new_pars)) updated_pars[markers_brand_new] <- new_pars[markers_brand_new]

  ## ---- 4. re-bootstrap ONLY touched markers -------------------------------
  touched_markers <- union(markers_known, markers_brand_new)   # every marker the new data measured
  unchanged_markers <- setdiff(names(old_th), touched_markers)

  th_touched <- GetBackground(
    data_norm = updated_norm[touched_markers, , drop = FALSE],
    dir_name  = tmp_dir, file.id = "touched",
    metadata  = updated_meta,
    sample_n  = num.cells.sample, parallelize = (ncores > 1), ncores = ncores,
    num.boot.sample = num.boot.sample, fixseed = fixseed,
    bootstrap_stratify = bootstrap_stratify, threshold_method = threshold_method
  )

  updated_th <- .merge_threshold_lists(old_th, th_touched, touched_markers)
  audit      <- .build_update_audit(old_th, updated_th, touched_markers, markers_brand_new)

  ## ---- 6. backup + write -----------------------------------------------
  if (isTRUE(backup)) {
    bdir <- file.path(ref_dir, paste0("backup_", format(Sys.time(), "%Y%m%dT%H%M%S")))
    dir.create(bdir, showWarnings = FALSE, recursive = TRUE)
    for (p in c(p_norm, p_meta, p_pars, p_th)) file.copy(p, file.path(bdir, basename(p)))
    message("[add_reference_dataset] Backed up previous reference files to: ", bdir)
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(updated_norm, file.path(out_dir, paste0(file.name, "_NormData.rds")))
  saveRDS(updated_meta, file.path(out_dir, paste0(file.name, "_Metadata.rds")))
  saveRDS(updated_pars, file.path(out_dir, paste0(file.name, "_bg_prot_pars.rds")))
  saveRDS(updated_th,   file.path(out_dir, paste0(file.name, "_th_dist.rds")))
  utils::write.csv(audit, file.path(out_dir, paste0(file.name, "_update_audit.csv")), row.names = FALSE)
  .write_provenance(out_dir, file.name,
                    bootstrap_stratify = bootstrap_stratify, threshold_method = threshold_method,
                    num.boot.sample = num.boot.sample, sample_n = num.cells.sample, fixseed = fixseed,
                    incremental_update = TRUE, n_touched_markers = length(touched_markers),
                    n_new_markers = length(markers_brand_new), n_unchanged_markers = length(unchanged_markers))

  message(sprintf(
    "[add_reference_dataset] Done. %d markers touched (%d brand new), %d markers left byte-for-byte unchanged. Audit: %s",
    length(touched_markers), length(markers_brand_new), length(unchanged_markers),
    file.path(out_dir, paste0(file.name, "_update_audit.csv"))))

  invisible(list(touched_markers = touched_markers, new_markers = markers_brand_new,
                unchanged_markers = unchanged_markers, threshold_shift = audit))
}

#' Merge an old threshold-distribution list with freshly-bootstrapped
#' results for just the touched markers, leaving every other marker
#' byte-for-byte identical to `old_th`. Pulled out of
#' `add_reference_dataset()` so this bookkeeping -- the actual "incremental"
#' guarantee -- can be unit-tested without ADTnorm
#' (`tests/testthat/test-add-reference-dataset.R`).
#' @param old_th Named list, as read from `<file.name>_th_dist.rds`.
#' @param th_touched Named list covering (at least) `touched_markers`, as
#'   returned by a `GetBackground()` call restricted to those markers.
#' @param touched_markers Character vector of marker names to overwrite/add.
#' @return Named list: same names/values as `old_th` for every marker NOT in
#'   `touched_markers`; `th_touched`'s value for every marker in
#'   `touched_markers` (including ones absent from `old_th`, i.e. brand new).
#' @noRd
.merge_threshold_lists <- function(old_th, th_touched, touched_markers) {
  missing_touched <- setdiff(touched_markers, names(th_touched))
  if (length(missing_touched) > 0) {
    stop(".merge_threshold_lists: th_touched is missing results for: ",
         paste(missing_touched, collapse = ", "))
  }
  updated_th <- old_th
  updated_th[touched_markers] <- th_touched[touched_markers]
  updated_th
}

#' Per-marker summary of how much a threshold moved, for markers touched by
#' an `add_reference_dataset()` call. `status` distinguishes markers that
#' existed before (`"existing_marker_updated"`) from ones the new dataset
#' introduced (`"new_marker"`, for which `old_mean_*`/`shift_*` are `NA`).
#' @noRd
.build_update_audit <- function(old_th, updated_th, touched_markers, markers_brand_new) {
  .mean_col <- function(th_entry, col) if (is.null(th_entry)) NA_real_ else mean(th_entry[, col])
  do.call(rbind, lapply(touched_markers, function(m) {
    is_new <- m %in% markers_brand_new
    old_l <- if (is_new) NA_real_ else .mean_col(old_th[[m]], "Lower")
    old_u <- if (is_new) NA_real_ else .mean_col(old_th[[m]], "Upper")
    new_l <- .mean_col(updated_th[[m]], "Lower")
    new_u <- .mean_col(updated_th[[m]], "Upper")
    data.frame(
      marker = m,
      status = if (is_new) "new_marker" else "existing_marker_updated",
      old_mean_lower = old_l, new_mean_lower = new_l,
      shift_lower    = if (is_new) NA_real_ else new_l - old_l,
      old_mean_upper = old_u, new_mean_upper = new_u,
      shift_upper    = if (is_new) NA_real_ else new_u - old_u,
      stringsAsFactors = FALSE
    )
  }))
}

#' @noRd
.resolve_new_metadata_cols <- function(metadata, cell_id_col, cell_type_col, dataset_col,
                                       sample_col, batch_col, tissue_col, tissue_default) {
  out <- metadata
  ren <- function(df, want, fallbacks, canon) {
    cand <- c(want, fallbacks)
    hit <- cand[cand %in% colnames(df)][1]
    if (is.na(hit)) return(df)
    colnames(df)[colnames(df) == hit] <- canon
    df
  }
  out <- ren(out, cell_id_col,   c("cell_id", "cell.id"), "cell.id")
  out <- ren(out, cell_type_col, c("cell_type", "cell.type", "cell.type_harmonised"), "cell.type")
  out <- ren(out, dataset_col,   c("dataset_id", "dataset"), "dataset")
  out <- ren(out, sample_col,    c("sample_id", "sample"), "sample")
  out <- ren(out, batch_col,     c("batch"), "batch")
  out <- ren(out, tissue_col,    c("tissue", "tissue_id"), "tissue")
  if (!"tissue" %in% colnames(out)) out$tissue <- tissue_default
  if (!"batch" %in% colnames(out)) out$batch <- out$sample
  out
}

#' @noRd
.read_provenance_if_exists <- function(dir_name, file.id) {
  p_json <- file.path(dir_name, paste0(file.id, "_provenance.json"))
  p_txt  <- file.path(dir_name, paste0(file.id, "_provenance.txt"))
  if (file.exists(p_json) && requireNamespace("jsonlite", quietly = TRUE)) {
    return(tryCatch(jsonlite::fromJSON(p_json), error = function(e) NULL))
  }
  if (file.exists(p_txt)) {
    lines <- readLines(p_txt)
    kv <- strsplit(lines, ": ", fixed = TRUE)
    vals <- setNames(vapply(kv, function(x) if (length(x) > 1) x[2] else NA_character_, character(1)),
                     vapply(kv, `[`, character(1), 1))
    return(as.list(vals))
  }
  NULL
}
