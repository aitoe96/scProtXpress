###############################################################################
# Each function here is deliberately generic: it takes already-normalized
# data + metadata (the same contract GetBackground() uses) so it works on
# ANY scProtXpress reference, not just the one in the associated manuscript.
###############################################################################

#' Leave-one-dataset-out reference threshold stability
#'
#' For every marker measured in at least `min_datasets_per_marker` distinct
#' datasets, rebuilds the bootstrap threshold distribution with each
#' contributing dataset held out in turn (reusing `GetBackground()`
#' unchanged), and reports how much the threshold moves. A marker whose
#' threshold barely changes regardless of which dataset is removed is
#' "underpinned by the whole reference"; a marker whose threshold swings a
#' lot when one specific dataset is removed is disproportionately
#' determined by that one dataset -- a useful diagnostic for reference
#' generalizability.
#'
#' @param data_norm proteins x cells matrix, ADTnorm-normalized (same
#'   contract as `GetBackground(data_norm = )`).
#' @param metadata data.frame, one row per cell, ALIGNED POSITIONALLY to
#'   `data_norm`'s columns (row i describes column i -- see
#'   `GetBackground()`'s own requirement), with a `dataset` column.
#' @param min_datasets_per_marker Integer, default 3. Markers measured in
#'   fewer datasets than this are excluded -- leave-one-out on a marker
#'   seen in only 1-2 datasets isn't a meaningful stability statement, and
#'   removing its only dataset would leave nothing to rebuild from.
#' @param markers Optional character vector to restrict the analysis to;
#'   default `NULL` (test every eligible marker in `data_norm`).
#' @param sample_n,num.boot.sample,bootstrap_stratify,threshold_method,ncores,fixseed
#'   Forwarded to every `GetBackground()` call (the full-reference one and
#'   every leave-one-out one) -- keep these identical to whatever your real
#'   reference build uses, or the comparison isn't apples-to-apples.
#' @return list(full_reference, per_dataset, stability_summary).
#'   `stability_summary` is a data.frame with one row per eligible marker:
#'   `n_datasets_loo` (how many leave-one-out runs contributed), the full-
#'   reference threshold, the mean/SD of the leave-one-out thresholds, the
#'   largest absolute shift from the full-reference threshold
#'   (`max_abs_shift_lower/upper`), and a coefficient of variation
#'   (`cv_lower/upper`) across the leave-one-out runs. Sort by
#'   `max_abs_shift_*` descending to find the markers most sensitive to a
#'   single dataset's inclusion.
#' @export
leave_one_dataset_out_stability <- function(
    data_norm, metadata,
    min_datasets_per_marker = 3,
    markers = NULL,
    sample_n = 100, num.boot.sample = 200,
    bootstrap_stratify = c("none", "cell.type", "cell.type_tissue"),
    threshold_method = c("rectangle", "mixture", "quantile", "zscore"),
    ncores = 1L, fixseed = 1234L
){
  bootstrap_stratify <- match.arg(bootstrap_stratify)
  threshold_method   <- match.arg(threshold_method)

  if (!"dataset" %in% colnames(metadata)) {
    stop("leave_one_dataset_out_stability: metadata must have a 'dataset' column ",
         "(the canonical name resolved by run_scProteinXpress()'s dataset_col).")
  }
  data_norm <- .as_dense_matrix(data_norm)
  if (ncol(data_norm) != nrow(metadata)) {
    stop("leave_one_dataset_out_stability: ncol(data_norm) must equal nrow(metadata) ",
         "(rows of metadata are matched POSITIONALLY to columns of data_norm, exactly ",
         "as GetBackground() requires).")
  }

  dataset_vec <- as.character(metadata$dataset)
  datasets <- unique(dataset_vec)
  if (length(datasets) < 3) {
    stop("leave_one_dataset_out_stability: need at least 3 distinct datasets for a ",
         "meaningful leave-one-out analysis; found ", length(datasets), ".")
  }

  marker_dataset_counts <- vapply(rownames(data_norm), function(m) {
    length(unique(dataset_vec[!is.na(data_norm[m, ])]))
  }, integer(1))
  eligible_markers <- names(marker_dataset_counts)[marker_dataset_counts >= min_datasets_per_marker]
  if (!is.null(markers)) eligible_markers <- intersect(eligible_markers, markers)
  if (length(eligible_markers) == 0) {
    stop("leave_one_dataset_out_stability: no markers are measured in >= ",
         min_datasets_per_marker, " distinct datasets (after intersecting with `markers`, if given).")
  }
  data_sub <- data_norm[eligible_markers, , drop = FALSE]

  .run_gb <- function(d_norm, meta, label) {
    GetBackground(d_norm, tempfile(paste0("loo_", label, "_")), label, meta,
                  sample_n = sample_n, parallelize = (ncores > 1), ncores = ncores,
                  num.boot.sample = num.boot.sample, fixseed = fixseed,
                  bootstrap_stratify = bootstrap_stratify, threshold_method = threshold_method)
  }

  full_th <- .run_gb(data_sub, metadata, "full")

  per_dataset <- list()
  for (d in datasets) {
    keep <- dataset_vec != d
    if (sum(keep) == 0) next
    sub_meta <- metadata[keep, , drop = FALSE]
    sub_data <- data_sub[, keep, drop = FALSE]
    still_present <- rownames(sub_data)[rowSums(!is.na(sub_data)) > 0]
    if (length(still_present) == 0) next
    per_dataset[[d]] <- .run_gb(sub_data[still_present, , drop = FALSE], sub_meta,
                                paste0("held_out_", make.names(d)))
  }
  if (length(per_dataset) == 0) {
    stop("leave_one_dataset_out_stability: every leave-one-out run ended up with zero ",
         "eligible markers -- check that datasets overlap in which markers they measure.")
  }

  stability_summary <- do.call(rbind, lapply(eligible_markers, function(m) {
    full_lower <- mean(full_th[[m]][, "Lower"]); full_upper <- mean(full_th[[m]][, "Upper"])
    loo_lower <- vapply(per_dataset, function(th) {
      if (is.null(th[[m]])) NA_real_ else mean(th[[m]][, "Lower"])
    }, numeric(1))
    loo_upper <- vapply(per_dataset, function(th) {
      if (is.null(th[[m]])) NA_real_ else mean(th[[m]][, "Upper"])
    }, numeric(1))
    loo_lower <- loo_lower[!is.na(loo_lower)]; loo_upper <- loo_upper[!is.na(loo_upper)]
    n_loo <- length(loo_lower)
    data.frame(
      marker = m,
      n_datasets_loo   = n_loo,
      full_lower       = full_lower,
      full_upper       = full_upper,
      mean_loo_lower   = if (n_loo > 0) mean(loo_lower) else NA_real_,
      sd_loo_lower     = if (n_loo > 1) sd(loo_lower) else NA_real_,
      mean_loo_upper   = if (n_loo > 0) mean(loo_upper) else NA_real_,
      sd_loo_upper     = if (n_loo > 1) sd(loo_upper) else NA_real_,
      max_abs_shift_lower = if (n_loo > 0) max(abs(loo_lower - full_lower)) else NA_real_,
      max_abs_shift_upper = if (n_loo > 0) max(abs(loo_upper - full_upper)) else NA_real_,
      cv_lower = if (n_loo > 1 && full_lower != 0) sd(loo_lower) / abs(full_lower) else NA_real_,
      cv_upper = if (n_loo > 1 && full_upper != 0) sd(loo_upper) / abs(full_upper) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(full_reference = full_th, per_dataset = per_dataset, stability_summary = stability_summary)
}
