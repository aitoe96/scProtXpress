###############################################################################
# internalFunctions.R -- scProtXpress
#
# ADTnorm-based protein normalisation backend shared by both modalities.
#
# Behaviour (preserved from the pre-package CyTOF reference implementation):
#   * .ADTprot_norm_ref  : per-marker foreach %dopar% (fork/PSOCK + BLAS pin).
#   * .ADTprot_norm_query: per-marker foreach %dopar% (fork/PSOCK + BLAS pin).
#   * RunBootstrap       : one-shot precompute + .RunBackComputations_fast
#                          inside foreach %dopar% (cell.type x tissue strata).
#   * resolve_ncores honours SLURM_CPUS_PER_TASK.
#   * Global QC with audit CSV (cleanData_ADT).
#   * Binomial sig_frame (Create.sig.frame + discrete.expr.levels).
#   * Self-query persistence handled by run_scProteinXpress().
###############################################################################

#' @import foreach
#' @importFrom Matrix Matrix
#' @importFrom stats density ecdf median optimize p.adjust pbinom quantile sd var rnorm setNames
#' @importFrom utils head read.csv read.delim write.csv write.table
NULL

if (!requireNamespace("tictoc", quietly = TRUE)) {
  tictoc <- list(tic = function(...) invisible(NULL), toc = function(...) invisible(NULL))
}

`%dopar%` <- foreach::`%dopar%`
`%do%`    <- foreach::`%do%`

# -----------------------------------------------------------------------------
# SQS (Stain Quality Score) - static helpers
# Return 0 if data is missing or variance is zero; never a fatal error.
# -----------------------------------------------------------------------------
.sqs_area_under_curve <- function(densityObj, threshold, peak_num){
  if (is.null(densityObj) || !is.finite(threshold)) return(0)
  if (peak_num > 1) {
    index <- which(densityObj$x >= threshold)
    if (length(index) == 0) return(0)
    auc <- sum(densityObj$y[index] * diff(c(densityObj$x[index], max(densityObj$x))))
    if (!is.finite(auc)) return(0)
    return(auc)
  } else {
    index <- which(densityObj$x >= threshold)
    if (length(index) == 0) return(0)
    auc <- sum(densityObj$y[index] * diff(c(densityObj$x[index], max(densityObj$x))))
    if (!is.finite(auc)) return(0)
    return(min(auc, 1 - auc))
  }
}

.sqs_valley_deep_scaler <- function(adt_tmp, peak_row, valley_row, adjust = 3, peak_num = NA){
  if (length(adt_tmp) < 2) return(0)
  adt_tmp <- adt_tmp[is.finite(adt_tmp)]
  if (length(adt_tmp) < 2 || stats::var(adt_tmp) <= 0) return(0)
  density_res <- tryCatch(stats::density(adt_tmp, adjust = adjust, na.rm = TRUE),
                          error = function(e) NULL)
  if (is.null(density_res)) return(0)
  x <- density_res$x; y <- density_res$y
  if (peak_num > 1) {
    peak_last <- peak_row[length(peak_row)]
    if (!is.finite(peak_last)) return(0)
    x_peakR <- which(abs(x - peak_last) == min(abs(x - peak_last)))[1]
    y_peakR <- y[x_peakR]
  } else {
    y_peakR <- 0
  }
  if (peak_num <= 2) {
    v_ref <- valley_row[1]
  } else {
    v_ref <- valley_row[length(valley_row)]
  }
  if (!is.finite(v_ref)) return(0)
  x_valley <- which(abs(x - v_ref) == min(abs(x - v_ref)))[1]
  y_valley <- y[x_valley]
  auc_scaler <- .sqs_area_under_curve(density_res, threshold = valley_row[1], peak_num = peak_num)
  deep_scaler <- (1 + y_peakR - y_valley) * (1 + auc_scaler)
  if (!is.finite(deep_scaler)) return(0)
  deep_scaler
}

.sqs_two_peak <- function(adt_tmp, peak_row, valley_row, peak_num = 2){
  adt_tmp <- adt_tmp[is.finite(adt_tmp)]
  if (length(adt_tmp) < 2) return(0)
  peak_first <- peak_row[1]; peak_last <- peak_row[length(peak_row)]
  v1 <- valley_row[1]
  if (!all(is.finite(c(peak_first, peak_last, v1)))) return(0)
  mean_diff <- peak_last - peak_first
  left_idx  <- which(adt_tmp < v1)
  right_idx <- which(adt_tmp > v1)
  ss <- 0
  if (length(left_idx)  > 0) ss <- ss + sum((adt_tmp[left_idx]  - peak_first)^2)
  if (length(right_idx) > 0) ss <- ss + sum((adt_tmp[right_idx] - peak_last )^2)
  within_peak_sd <- sqrt(ss / length(adt_tmp))
  if (!is.finite(within_peak_sd) || within_peak_sd <= 0) return(0)
  deep_scaler <- .sqs_valley_deep_scaler(adt_tmp, peak_row, valley_row, peak_num = peak_num)
  score <- (mean_diff * deep_scaler) / within_peak_sd
  if (!is.finite(score)) return(0)
  score
}

.sqs_multi_peak <- function(adt_tmp, peak_row, valley_row, peak_num = 3){
  adt_tmp <- adt_tmp[is.finite(adt_tmp)]
  if (length(adt_tmp) < 2) return(0)
  peak_first <- peak_row[1]; peak_last <- peak_row[length(peak_row)]
  if (!all(is.finite(c(peak_first, peak_last)))) return(0)
  mean_diff <- peak_last - peak_first
  within_peak_var <- 0
  ncols <- length(peak_row)
  for (col_index in seq_len(ncols)) {
    pk <- peak_row[col_index]
    if (!is.finite(pk)) next
    if (col_index == 1) {
      v_right <- valley_row[col_index]
      if (!is.finite(v_right)) next
      sel <- which(adt_tmp < v_right)
    } else if (col_index == ncols) {
      v_left <- valley_row[col_index - 1]
      if (!is.finite(v_left)) next
      sel <- which(adt_tmp > v_left)
    } else {
      v_left <- valley_row[col_index - 1]; v_right <- valley_row[col_index]
      if (!is.finite(v_left) || !is.finite(v_right)) next
      sel <- which(adt_tmp < v_right & adt_tmp > v_left)
    }
    if (length(sel) > 0) within_peak_var <- within_peak_var + sum((adt_tmp[sel] - pk)^2)
  }
  within_peak_sd <- sqrt(within_peak_var / length(adt_tmp))
  if (!is.finite(within_peak_sd) || within_peak_sd <= 0) return(0)
  deep_scaler <- .sqs_valley_deep_scaler(adt_tmp, peak_row, valley_row, peak_num = peak_num)
  score <- (mean_diff * deep_scaler) / within_peak_sd
  if (!is.finite(score)) return(0)
  score
}

.sqs_one_peak <- function(adt_tmp, peak_row, valley_row, peak_num = 1){
  adt_tmp <- adt_tmp[is.finite(adt_tmp)]
  if (length(adt_tmp) < 2) return(0)
  v1 <- valley_row[1]
  if (!is.finite(v1)) return(0)
  ncols <- length(peak_row)
  if (ncols == 1) {
    p_ref <- peak_row[1]
  } else if (is.na(peak_row[1])) {
    p_ref <- peak_row[ncols]
  } else {
    p_ref <- peak_row[1]
  }
  if (!is.finite(p_ref)) return(0)
  mean_diff <- abs(v1 - p_ref)
  sd_val <- stats::sd(adt_tmp)
  if (!is.finite(sd_val) || sd_val <= 0) return(0)
  within_peak_sd <- sqrt(sd_val)
  if (!is.finite(within_peak_sd) || within_peak_sd <= 0) return(0)
  deep_scaler <- .sqs_valley_deep_scaler(adt_tmp, peak_row, valley_row, peak_num = peak_num)
  score <- (mean_diff * deep_scaler) / within_peak_sd
  if (!is.finite(score)) return(0)
  score
}

# Dispatcher: computes SQS per sample for a marker, returns a named numeric.
.sqs_compute_per_sample <- function(marker_vec, sample_vec, valid_samples, peak_info, valley_info){
  out <- setNames(rep(0, length(valid_samples)), valid_samples)
  if (is.null(peak_info) || is.null(valley_info)) return(out)
  for (s in valid_samples) {
    if (!(s %in% rownames(peak_info)) || !(s %in% rownames(valley_info))) { out[s] <- 0; next }
    peak_row   <- peak_info[s, ]
    valley_row <- valley_info[s, ]
    peak_num   <- sum(!is.na(peak_row))
    idx_s <- which(sample_vec == s)
    if (length(idx_s) == 0) { out[s] <- 0; next }
    adt_tmp <- marker_vec[idx_s]
    if (peak_num <= 0) { out[s] <- 0; next }
    score <- tryCatch({
      if (peak_num == 1) {
        .sqs_one_peak(adt_tmp, peak_row, valley_row, peak_num = peak_num)
      } else if (peak_num == 2) {
        .sqs_two_peak(adt_tmp, peak_row, valley_row, peak_num = peak_num)
      } else {
        .sqs_multi_peak(adt_tmp, peak_row, valley_row, peak_num = peak_num)
      }
    }, error = function(e) 0)
    if (!is.finite(score)) score <- 0
    out[s] <- score
  }
  out
}

# Append-only CSV logger (no dependencies; thread-unsafe but used sequentially).
.sqs_log_append <- function(log_path, df_row){
  if (is.null(log_path) || !nzchar(log_path)) return(invisible(NULL))
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  write_header <- !file.exists(log_path)
  utils::write.table(df_row, file = log_path, sep = ",",
                     row.names = FALSE, col.names = write_header,
                     append = !write_header, quote = FALSE)
  invisible(NULL)
}

# Internal timestamped logger.
# Returns invisibly so it can be embedded in pipelines without affecting the
# return values of the wrapping block. The format is ISO-8601 down to the
# second to allow log greps and downstream parsing.
#
# Intentionally NOT a replacement for the existing message()/cat() calls
# elsewhere in this file, which keep their wording and order unchanged.
# `.log_step()` is opt-in scaffolding for new instrumentation only.
.log_step <- function(label, ..., level = c("info", "warn"), .stream = stderr()) {
  level <- match.arg(level)
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  msg <- sprintf("[%s][%s] %s", ts, toupper(level),
                 paste0(c(label, ...), collapse = " "))
  if (identical(level, "warn")) {
    warning(msg, call. = FALSE)
  } else {
    message(msg)
  }
  invisible(msg)
}

# -----------------------------------------------------------------------------
# Modality dispatch helpers.
#
# `.ingest_<modality>()` is invoked from run_scProteinXpress() at the very
# top of the pipeline. Its sole responsibility is to coerce `data` into the
# storage backend that the rest of the pipeline expects for that modality,
# and to return `metadata` as a base data.frame. Both helpers leave
# metadata column translation and validation to run_scProteinXpress().
#
#   CyTOF    : dense base::matrix (panels are small; arcsinh values dense).
#   CITE-seq : sparse dgCMatrix   (panels can be large + raw ADT counts are
#              sparse with many zeros; preserving sparsity matters for RAM).
#
# Downstream helpers (`.as_dense_matrix`, `make_ADTnorm_input`,
# `cleanData_ADT`) accept both backends because `.as_dense_matrix()`
# already handles `Matrix` objects via `as.matrix()`.
# -----------------------------------------------------------------------------
.ingest_cytof <- function(data, metadata) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (inherits(data, "Matrix")) data <- as.matrix(data)
  if (!is.matrix(data)) {
    stop(".ingest_cytof: data must be matrix, Matrix, or data.frame; got ",
         paste(class(data), collapse = ","))
  }
  storage.mode(data) <- "double"
  list(data = data, metadata = as.data.frame(metadata, stringsAsFactors = FALSE))
}

.ingest_citeseq <- function(data, metadata) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (is.matrix(data)) {
    data <- methods::as(Matrix::Matrix(data, sparse = TRUE), "dgCMatrix")
  } else if (inherits(data, "Matrix") && !inherits(data, "dgCMatrix")) {
    data <- methods::as(data, "dgCMatrix")
  }
  if (!inherits(data, "dgCMatrix")) {
    stop(".ingest_citeseq: data must be matrix, Matrix, dgCMatrix, or data.frame; got ",
         paste(class(data), collapse = ","))
  }
  list(data = data, metadata = as.data.frame(metadata, stringsAsFactors = FALSE))
}

# `.resolve_meta_col()` standardises metadata column resolution:
# - Tries `primary` first; then each entry in `fallbacks` in order.
# - First column found is renamed in place to `internal_name`.
# - If none match: either fall back to a constant (`fallback_value`),
#   to another already-existing column (`fallback_from`), or stop().
.resolve_meta_col <- function(metadata, primary, fallbacks, internal_name,
                              fallback_value = NA, fallback_from = NULL) {
  candidates <- unique(c(primary, fallbacks))
  hit <- candidates[candidates %in% colnames(metadata)][1]
  if (!is.na(hit) && length(hit) == 1L) {
    if (hit != internal_name) {
      ## If `internal_name` already exists too, drop it before renaming so we
      ## do not create duplicated colnames.
      if (internal_name %in% colnames(metadata)) {
        metadata[[internal_name]] <- NULL
      }
      colnames(metadata)[colnames(metadata) == hit] <- internal_name
    }
    return(metadata)
  }
  if (!is.null(fallback_from) && fallback_from %in% colnames(metadata)) {
    metadata[[internal_name]] <- metadata[[fallback_from]]
    return(metadata)
  }
  if (!isTRUE(is.na(fallback_value))) {
    metadata[[internal_name]] <- fallback_value
    return(metadata)
  }
  stop(sprintf(
    "metadata: could not resolve column for '%s'. Tried: %s",
    internal_name, paste(candidates, collapse = ", ")))
}

# --- Helpers ---
.as_dense_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (inherits(x, "Matrix")) return(as.matrix(x))
  if (is.data.frame(x)) return(as.matrix(x))
  stop("Cannot coerce object to dense matrix", call. = FALSE)
}

.pad_to_features <- function(mat, feats) {
  mat <- .as_dense_matrix(mat)
  miss <- setdiff(feats, rownames(mat))
  if (length(miss) > 0) {
    zero <- matrix(NA_real_, nrow = length(miss), ncol = ncol(mat))
    rownames(zero) <- miss
    colnames(zero) <- colnames(mat)
    mat <- rbind(mat, zero)
  }
  mat <- mat[feats, , drop = FALSE]
  .as_dense_matrix(mat)
}

.safe_to_matrix <- function(res, rnames) {
  if (is.null(res)) return(matrix(NA_real_, nrow=length(rnames), ncol=1, dimnames=list(rnames, NULL)))
  if (is.data.frame(res)) res <- as.matrix(res)
  if (!is.matrix(res)) res <- matrix(res, nrow = length(rnames), byrow = TRUE)
  rownames(res) <- rnames
  res
}

make_ADTnorm_input <- function(counts_fxC, metadata, sample_col = "sample", batch_col = "batch") {
  counts_fxC <- .as_dense_matrix(counts_fxC)
  metadata <- as.data.frame(metadata)
  common.cells <- intersect(colnames(counts_fxC), rownames(metadata))
  counts_fxC <- counts_fxC[, common.cells, drop = FALSE]
  metadata <- metadata[common.cells, , drop = FALSE]
  
  cell_x_adt <- t(counts_fxC)
  cell_x_feature <- metadata
  cell_x_feature$sample <- as.character(cell_x_feature[[sample_col]])
  
  if (batch_col %in% colnames(cell_x_feature)) {
    cell_x_feature$batch <- as.character(cell_x_feature[[batch_col]])
  } else {
    cell_x_feature$batch <- cell_x_feature$sample
  }
  
  if (!is.factor(cell_x_feature$sample)) cell_x_feature$sample <- factor(cell_x_feature$sample, levels = unique(cell_x_feature$sample))
  cell_x_feature <- cell_x_feature[rownames(cell_x_adt), , drop = FALSE]
  list(cell_x_adt = cell_x_adt, cell_x_feature = cell_x_feature)
}

## .detect_input_mode() -- three-valued is_transformed semantics.
##
##   is_transformed = NULL   (default)
##     "decide for me": run the integer-coercibility heuristic. Treat the
##     matrix as already transformed iff at least one usable value differs
##     from its integer cast (i.e. there is a non-integer real value).
##
##   is_transformed = TRUE
##     "force: already transformed". Skip arcsinh regardless of values.
##
##   is_transformed = FALSE
##     "force: NOT transformed". Always apply arcsinh. This bypasses the
##     heuristic, so callers that *know* the matrix is raw can override
##     ADTnorm's auto-detection of non-integer noise.
##
## Returns:
##   list(arcsinh_transform_flag, neg_candidate_thres)
##
## Note on backwards compatibility: the pre-unification API treated FALSE
## as "autodetect" because the heuristic ran via `is_transformed || any(...)`.
## With the new semantics, FALSE means force-apply. The CyTOF driver
## passes `is_transformed = TRUE` explicitly for the arcsinh-transformed
## "exprs" assay; it is unaffected. Drivers that previously omitted the
## argument (CITE-seq `codebase/run_scProteinXpress.R`) now get NULL via
## the unified API, which preserves the prior autodetect behaviour.
#' Distributional sanity check for the "already transformed" branch
#'
#' `.detect_input_mode()`'s autodetection is purely syntactic (is the matrix
#' integer-coercible?), and an explicit `is_transformed = TRUE` is taken on
#' faith. Neither checks that "already transformed" actually means
#' "arcsinh-transformed", which is what the reference thresholds assume.
#' Public CITE-seq deposits frequently ship DSB- or CLR-normalised ADT
#' values instead (roughly Gaussian, substantial negative mass) -- those
#' pass this gate silently and get scored against arcsinh-based reference
#' thresholds on a mismatched scale with no warning anywhere else in the
#' pipeline. arcsinh(raw_count/5 + 1)-style transforms of non-negative raw
#' input are essentially always >= 0, so a simple, cheap proxy is the
#' fraction of negative values.
#' @noRd
.sanity_check_transformed_scale <- function(matrix_in, source = c("autodetect", "forced")){
  source <- match.arg(source)
  vals <- matrix_in[!is.na(matrix_in)]
  if (length(vals) == 0) return(invisible(NA_real_))
  frac_negative <- mean(vals < 0)
  if (frac_negative > 0.05) {
    warning(sprintf(
      paste0(
        ".detect_input_mode: input was %s as already-transformed, but %.1f%% ",
        "of finite values are negative. arcsinh-transformed raw counts (what this ",
        "package's reference thresholds assume) are essentially always >= 0 for ",
        "non-negative raw input -- this looks more consistent with a DSB- or ",
        "CLR-normalised (or z-scored) input than an arcsinh transform. If so, ",
        "query-time comparisons against the arcsinh-based reference thresholds will ",
        "silently run on a mismatched scale. Re-derive from raw counts + arcsinh, or ",
        "build/query with data normalised the same way throughout."
      ), source, 100 * frac_negative), call. = FALSE)
  }
  invisible(frac_negative)
}

.detect_input_mode <- function(cell_x_adt, neg_candidate_thres = NULL,
                               is_transformed = NULL, verbose = TRUE) {
  matrix_in <- as.matrix(cell_x_adt)

  if (is.null(is_transformed)) {
    as_int <- matrix_in
    mode(as_int) <- "integer"
    is_already_transformed <- any(
      as_int[!is.na(as_int)] != matrix_in[!is.na(matrix_in)],
      na.rm = TRUE
    )
    detect_source <- "autodetect"
  } else if (isTRUE(is_transformed)) {
    is_already_transformed <- TRUE
    detect_source <- "forced"
  } else if (isFALSE(is_transformed)) {
    is_already_transformed <- FALSE
    detect_source <- "forced"
  } else {
    stop(".detect_input_mode: is_transformed must be NULL, TRUE, or FALSE.")
  }

  if (is_already_transformed) {
    arcsinh_transform_flag <- FALSE
    if (is.null(neg_candidate_thres)) {
      neg_candidate_thres <- asinh(1/5 + 1)
    }
    .sanity_check_transformed_scale(matrix_in, source = detect_source)
  } else {
    if (any(matrix_in[!is.na(matrix_in)] < 0)) {
      warning("ADTnorm detects integer input ADT count matrix and assumes the input is raw count. Please check why there is any negative values in the raw count matrix.")
    }
    arcsinh_transform_flag <- TRUE
    if (is.null(neg_candidate_thres)) neg_candidate_thres <- asinh(4/5 + 1)
  }
  list(arcsinh_transform_flag = arcsinh_transform_flag,
       neg_candidate_thres    = neg_candidate_thres)
}

.build_target_landmark <- function(target_landmark_location, landmark_align_type, peak_mode_res){
  if (is.null(target_landmark_location)) return(NULL)
  if (length(target_landmark_location) == 1 && identical(target_landmark_location, "fixed")) target_landmark_location <- c(1, 5)
  if (landmark_align_type == "negPeak_valley_posPeak") {
    if (ncol(peak_mode_res) == 1) return(c(target_landmark_location[1], round(mean(target_landmark_location), 1)))
    else if (ncol(peak_mode_res) == 2) return(c(target_landmark_location[1], round(mean(target_landmark_location), 1), target_landmark_location[2]))
    else return(c(target_landmark_location[1], (target_landmark_location[1] + round(mean(target_landmark_location), 1)) / 2, target_landmark_location[2]))
  }
  return(NULL)
}

.coerce_peak_mode_to_reference <- function(peak_mode_res, ref_n) {
  if (is.null(dim(peak_mode_res))) peak_mode_res <- matrix(peak_mode_res, nrow = length(peak_mode_res))
  if (ref_n == 2) peak_mode_res <- peak_mode_res[, 1, drop = FALSE]
  if (ref_n == 3) {
    if (ncol(peak_mode_res) == 1) peak_mode_res <- cbind(peak_mode_res, matrix(NA, nrow = nrow(peak_mode_res), ncol = 1))
    else if (ncol(peak_mode_res) > 2) peak_mode_res <- peak_mode_res[, c(1, ncol(peak_mode_res)), drop = FALSE]
  }
  peak_mode_res
}

.coerce_landmark_matrix_to_reference <- function(landmark_matrix, fixed_target_landmark, neg_candidate_thres = NULL){
  if (is.null(dim(landmark_matrix))) landmark_matrix <- matrix(landmark_matrix, nrow = length(landmark_matrix))
  ref_n <- length(fixed_target_landmark)
  if (ncol(landmark_matrix) < ref_n) landmark_matrix <- cbind(landmark_matrix, matrix(NA, nrow = nrow(landmark_matrix), ncol = ref_n - ncol(landmark_matrix)))
  else if (ncol(landmark_matrix) > ref_n) landmark_matrix <- landmark_matrix[, seq_len(ref_n), drop = FALSE]
  
  for (j in seq_len(ncol(landmark_matrix))) {
    idx_na <- is.na(landmark_matrix[, j])
    if (any(idx_na)) {
      col_non_na <- landmark_matrix[!idx_na, j]
      fill_val <- if (length(col_non_na) > 0) stats::median(col_non_na, na.rm = TRUE) else fixed_target_landmark[min(j, length(fixed_target_landmark))]
      landmark_matrix[idx_na, j] <- fill_val
    }
  }
  if (ncol(landmark_matrix) >= 2) {
    bad <- landmark_matrix[, 2] <= landmark_matrix[, 1]
    landmark_matrix[bad, 2] <- landmark_matrix[bad, 1] + 1e-3
  }
  landmark_matrix
}

.extract_align_field <- function(align_res, field, fallback = NULL, pos = NULL) {
  if (is.list(align_res) && field %in% names(align_res)) return(align_res[[field]])
  if (!is.null(pos) && length(align_res) >= pos) return(align_res[[pos]])
  fallback
}

# -----------------------------------------------------------------------------
# .run_one_marker_ref  
# -----------------------------------------------------------------------------
.run_one_marker_ref <- function(
    adt_marker_index,
    cell_x_adt,
    cell_x_feature,
    original_marker_name,
    quantile_clip,
    bw_smallest_bi,
    bw_smallest_tri,
    bw_smallest_adjustments,
    bimodal_marker_index,
    trimodal_marker_index,
    positive_peak,
    peak_type,
    neg_candidate_thres,
    lower_peak_thres,
    arcsinh_transform_flag,
    shoulder_valley,
    multi_sample_per_batch,
    valley_density_adjust,
    shoulder_valley_slope,
    detect_outlier_valley,
    sqs_threshold,
    apply_sqs,
    file.id,
    landmark_align_type,
    midpoint_type,
    target_landmark_location
){
  adt_marker_select      <- colnames(cell_x_adt)[adt_marker_index]      
  adt_marker_select_name <- original_marker_name[adt_marker_index]      
  
  out_skip_pre_sqs <- list(
    marker_name       = adt_marker_select_name,
    tmp_name          = adt_marker_select,
    norm_col_override = NULL,
    pars              = NULL,
    sqs_rows          = list()
  )
  
  if (quantile_clip < 1) {
    quant <- stats::quantile(cell_x_adt[, adt_marker_select], quantile_clip, na.rm = TRUE)
    clip_idx <- !is.na(cell_x_adt[, adt_marker_select]) & (cell_x_adt[, adt_marker_select] > quant)
    cell_x_adt[clip_idx, adt_marker_select] <- NA
  }
  
  valid_samples <- c()
  for (s in levels(cell_x_feature$sample)) {
    vals <- cell_x_adt[cell_x_feature$sample == s, adt_marker_select]
    valid_vals <- vals[is.finite(vals)]
    if (length(valid_vals) >= 5 && length(unique(valid_vals)) >= 3 && stats::var(valid_vals) > 1e-6) {
      valid_samples <- c(valid_samples, s)
    }
  }
  if (length(valid_samples) == 0) return(out_skip_pre_sqs)
  
  valid_idx <- cell_x_feature$sample %in% valid_samples
  cell_x_adt_sub     <- cell_x_adt[valid_idx, , drop = FALSE]
  cell_x_feature_sub <- cell_x_feature[valid_idx, , drop = FALSE]
  cell_x_feature_sub$sample <- factor(cell_x_feature_sub$sample, levels = valid_samples)
  
  bwFac_smallest <- bw_smallest_bi
  
  peak_mode_res <- tryCatch({
    if (peak_type == "mode") {
      .get_peak_mode_patched(cell_x_adt_sub, cell_x_feature_sub, adt_marker_select, adt_marker_index,
                    bwFac_smallest, bimodal_marker_index, trimodal_marker_index, positive_peak,
                    neg_candidate_thres, lower_peak_thres, arcsinh_transform_flag)
    } else {
      ADTnorm::get_peak_midpoint(cell_x_adt_sub, cell_x_feature_sub, adt_marker_select, adt_marker_index,
                        bwFac_smallest, bimodal_marker_index, trimodal_marker_index, positive_peak,
                        neg_candidate_thres, lower_peak_thres, arcsinh_transform_flag)
    }
  }, error = function(e) {
    message(paste0("\n[!] Skipping normalization for ", adt_marker_select_name,
                   ": Peak detection failed -> ", conditionMessage(e)))
    NULL
  })
  if (is.null(peak_mode_res)) return(out_skip_pre_sqs)
  
  peak_valley_list <- tryCatch({
    ADTnorm::get_valley_location(cell_x_adt_sub, cell_x_feature_sub, adt_marker_select, peak_mode_res,
                        shoulder_valley, positive_peak, multi_sample_per_batch,
                        valley_density_adjust, 20, shoulder_valley_slope,
                        neg_candidate_thres, arcsinh_transform_flag)
  }, error = function(e) {
    message(paste0("\n[!] Skipping normalization for ", adt_marker_select_name, ": ",
                   conditionMessage(e)))
    NULL
  })
  if (is.null(peak_valley_list)) return(out_skip_pre_sqs)
  
  valley_location_res <- peak_valley_list$valley_location_list
  peak_landmark_list  <- peak_valley_list$peak_landmark_list
  
  if (detect_outlier_valley != FALSE) {
    valley_location_res <- tryCatch({
      ADTnorm::detect_impute_outlier_valley(valley_location_res, adt_marker_select, cell_x_adt_sub,
                                   cell_x_feature_sub, 3, detect_outlier_valley, 3, 0.75)
    }, error = function(e) valley_location_res)
  }
  
  sqs_rows <- list()
  norm_col_override <- cell_x_adt[, adt_marker_select]
  
  peak_info_raw   <- .safe_to_matrix(peak_landmark_list,  valid_samples)
  valley_info_raw <- .safe_to_matrix(valley_location_res, valid_samples)
  
  sqs_scores <- .sqs_compute_per_sample(
    marker_vec    = cell_x_adt_sub[, adt_marker_select],
    sample_vec    = as.character(cell_x_feature_sub$sample),
    valid_samples = valid_samples,
    peak_info     = peak_info_raw,
    valley_info   = valley_info_raw
  )
  
  for (s in valid_samples) {
    sqs_rows[[length(sqs_rows) + 1L]] <- data.frame(
      stage      = "ref",
      file.id    = file.id,
      marker     = adt_marker_select_name,
      sample     = s,
      sqs        = unname(sqs_scores[s]),
      threshold  = sqs_threshold,
      apply_sqs  = apply_sqs,
      action     = if (apply_sqs %in% c("ref_only", "both") && sqs_scores[s] < sqs_threshold) "NA" else "keep",
      stringsAsFactors = FALSE
    )
  }
  
  if (apply_sqs %in% c("ref_only", "both")) {
    fail_samples <- names(sqs_scores)[sqs_scores < sqs_threshold]
    if (length(fail_samples) > 0) {
      fail_idx <- cell_x_feature$sample %in% fail_samples
      norm_col_override[fail_idx] <- NA_real_
      
      keep_mask           <- !(valid_samples %in% fail_samples)
      valid_samples       <- valid_samples[keep_mask]
      peak_landmark_list  <- peak_info_raw[keep_mask, , drop = FALSE]
      valley_location_res <- valley_info_raw[keep_mask, , drop = FALSE]
      
      if (length(valid_samples) == 0) {
        return(list(
          marker_name       = adt_marker_select_name,
          tmp_name          = adt_marker_select,
          norm_col_override = norm_col_override,
          pars              = NULL,
          sqs_rows          = sqs_rows
        ))
      }
      
      valid_idx <- cell_x_feature$sample %in% valid_samples
      cell_x_adt_sub     <- cell_x_adt[valid_idx, , drop = FALSE]
      cell_x_feature_sub <- cell_x_feature[valid_idx, , drop = FALSE]
      cell_x_feature_sub$sample <- factor(cell_x_feature_sub$sample, levels = valid_samples)
    }
  }
  
  full_samples <- levels(cell_x_feature$sample)
  peak_landmark_list <- .safe_to_matrix(peak_landmark_list, valid_samples)
  pm_full <- matrix(NA_real_, nrow = length(full_samples), ncol = ncol(peak_landmark_list))
  rownames(pm_full) <- full_samples
  pm_full[valid_samples, ] <- peak_landmark_list
  peak_mode_res_final <- pm_full
  
  valley_location_res <- .safe_to_matrix(valley_location_res, valid_samples)
  vl_full <- matrix(NA_real_, nrow = length(full_samples), ncol = ncol(valley_location_res))
  rownames(vl_full) <- full_samples
  vl_full[valid_samples, ] <- valley_location_res
  valley_location_res_final <- vl_full
  
  landmark_matrix <- tryCatch({
    ADTnorm::landmark_fill_na(peak_mode_res_final, valley_location_res_final,
                     landmark_align_type, midpoint_type, neg_candidate_thres)
  }, error = function(e) {
    message(paste0("\n[!] Skipping normalization for ", adt_marker_select_name,
                   ": Landmark alignment failed -> ", conditionMessage(e)))
    NULL
  })
  if (is.null(landmark_matrix)) {
    return(list(
      marker_name       = adt_marker_select_name,
      tmp_name          = adt_marker_select,
      norm_col_override = norm_col_override,
      pars              = NULL,
      sqs_rows          = sqs_rows
    ))
  }
  
  target_landmark <- .build_target_landmark(target_landmark_location,
                                            landmark_align_type, peak_mode_res_final)
  
  align_res <- tryCatch({
    ADTnorm::peak_alignment(cell_x_adt_sub[, adt_marker_select], cell_x_feature_sub,
                   landmark_matrix[valid_samples, , drop = FALSE], target_landmark)
  }, error = function(e) {
    message(paste0("\n[!] Skipping normalization for ", adt_marker_select_name,
                   ": Alignment failed -> ", conditionMessage(e)))
    NULL
  })
  if (is.null(align_res)) {
    return(list(
      marker_name       = adt_marker_select_name,
      tmp_name          = adt_marker_select,
      norm_col_override = norm_col_override,
      pars              = NULL,
      sqs_rows          = sqs_rows
    ))
  }
  
  norm_col_override[valid_idx] <- .extract_align_field(align_res, "cell_x_adt_norm", pos = 1)
  
  pars <- list(
    tmp_name              = adt_marker_select,
    fixed_target_landmark = .extract_align_field(align_res, "target_landmark_used",
                                                 fallback = target_landmark),
    warp_parameters       = .extract_align_field(align_res, "warp_parameters"),
    bwFac_smallest        = bwFac_smallest,
    peak_type             = peak_type,
    lower_peak_thres      = lower_peak_thres,
    neg_candidate_thres   = neg_candidate_thres,
    shoulder_valley       = shoulder_valley,
    shoulder_valley_slope = shoulder_valley_slope,
    valley_density_adjust = valley_density_adjust,
    multi_sample_per_batch= multi_sample_per_batch,
    midpoint_type         = midpoint_type,
    detect_outlier_valley = detect_outlier_valley,
    positive_peak         = positive_peak
  )
  
  list(
    marker_name       = adt_marker_select_name,
    tmp_name          = adt_marker_select,
    norm_col_override = norm_col_override,
    pars              = pars,
    sqs_rows          = sqs_rows
  )
}
# -----------------------------------------------------------------------------
# REFERENCE MODE
# -----------------------------------------------------------------------------
.ADTprot_norm_ref <- function(
    data_counts, 
    metadata, 
    file.id, 
    outdir, 
    sample_col, 
    batch_col, 
    marker_to_process,
    exclude_zeroes, 
    bimodal_marker, 
    trimodal_marker, 
    positive_peak, 
    bw_smallest_bi,
    bw_smallest_tri, 
    bw_smallest_adjustments, 
    quantile_clip, 
    peak_type, 
    multi_sample_per_batch,
    shoulder_valley, 
    shoulder_valley_slope, 
    valley_density_adjust, 
    landmark_align_type,
    midpoint_type, 
    neg_candidate_thres, 
    lower_peak_thres, 
    detect_outlier_valley,
    target_landmark_location, 
    clean_adt_name, 
    override_landmark, 
    verbose, 
    saveData,
    is_transformed,
    sqs_threshold = 1.0,
    apply_sqs = "none",
    do_subsample = TRUE,
    target_ct = 2000L,
    # NOTE: these defaults reflect the internal names that run_scProteinXpress
    # already enforces via the "internal metadata translation" (cell.type / tissue /
    # dataset). If .ADTprot_norm_ref is called with metadata that still uses the
    # external names (cell_type / dataset_id), they must be passed explicitly.
    ct_col = "cell.type",
    tissue_col = "tissue",
    dataset_col = "dataset",
    seed = 1234,
    ncores = 1L){
  
  data_counts <- .as_dense_matrix(data_counts)
  metadata <- as.data.frame(metadata)
  
  adt_in <- make_ADTnorm_input(data_counts, metadata, sample_col, batch_col)
  cell_x_adt <- adt_in$cell_x_adt
  cell_x_feature <- adt_in$cell_x_feature
  
  mode_info <- .detect_input_mode(cell_x_adt, neg_candidate_thres, is_transformed, verbose)
  arcsinh_transform_flag <- mode_info$arcsinh_transform_flag
  neg_candidate_thres <- mode_info$neg_candidate_thres
  
  all_marker_name <- colnames(cell_x_adt)
  col_sums <- colSums(cell_x_adt, na.rm = TRUE)
  if (any(col_sums == 0)) {
    keep <- col_sums > 0
    cell_x_adt <- cell_x_adt[, keep, drop = FALSE]
    all_marker_name <- colnames(cell_x_adt)
  }
  
  if (arcsinh_transform_flag) {
    if(verbose) message("Applying arcsinh transformation...")
    cell_x_adt <- ADTnorm::arcsinh_transform(cell_x_adt = cell_x_adt)
  }
  
  if (!is.factor(cell_x_feature$sample)) cell_x_feature$sample <- factor(cell_x_feature$sample, levels = unique(cell_x_feature$sample))
  
  bimodal_marker_index <- which(all_marker_name %in% bimodal_marker)
  trimodal_marker_index <- if (!is.null(trimodal_marker)) which(all_marker_name %in% trimodal_marker) else NULL
  
  original_marker_name <- all_marker_name
  colnames(cell_x_adt) <- paste0("tmpName", seq_len(ncol(cell_x_adt)))
  adt_marker_index_list <- seq_len(ncol(cell_x_adt))
  
  cell_x_adt_norm <- cell_x_adt
  
  bg.prot.pars <- list(
    method = "ADTnorm-compact", file.id = file.id, sample_col = sample_col, batch_col = batch_col,
    original_marker_name = original_marker_name,
    global_settings = list(
      exclude_zeroes = exclude_zeroes, quantile_clip = quantile_clip, peak_type = peak_type,
      multi_sample_per_batch = multi_sample_per_batch, shoulder_valley = shoulder_valley,
      shoulder_valley_slope = shoulder_valley_slope, valley_density_adjust = valley_density_adjust,
      landmark_align_type = landmark_align_type, midpoint_type = midpoint_type,
      neg_candidate_thres = neg_candidate_thres, lower_peak_thres = lower_peak_thres,
      detect_outlier_valley = detect_outlier_valley, target_landmark_location = target_landmark_location,
      bimodal_marker = bimodal_marker, trimodal_marker = trimodal_marker, positive_peak = positive_peak,
      bw_smallest_bi = bw_smallest_bi, bw_smallest_tri = bw_smallest_tri, bw_smallest_adjustments = bw_smallest_adjustments,
      is_transformed_input = is_transformed, 
      sqs_threshold = sqs_threshold, apply_sqs = apply_sqs
    ),
    markers = list()
  )
  
  sqs_log_path <- file.path(outdir, paste0(file.id, "_SQS_ref_log.csv"))
  
  # Cap ncores to the number of markers: no sense in requesting more workers than tasks.
  effective_ncores <- if (is.numeric(ncores) && length(ncores) == 1L)
    min(as.integer(ncores), length(adt_marker_index_list)) else 1L
  
  cl <- NULL
  use_parallel <- effective_ncores >= 2L && length(adt_marker_index_list) >= 2L
  
  # --- Pin BLAS to 1 thread in the master process and in each worker (KEY) -----
  # Without this, OpenBLAS/MKL launch implicit multi-threading in density/splines/etc.
  # With N workers * M BLAS threads the oversubscription kills the speedup.
  # We save the previous value to restore it on exit.
  .prev_blas_threads <- NA_integer_
  .prev_omp_threads  <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  .prev_mkl_threads  <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
  .prev_openblas_t   <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  if (use_parallel) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      .prev_blas_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
                                     error = function(e) NA_integer_)
      try(RhpcBLASctl::blas_set_num_threads(1L),  silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(1L),   silent = TRUE)
    }
    # Fallback via env-vars (affects fork workers; PSOCK inherits them).
    Sys.setenv(OMP_NUM_THREADS      = "1")
    Sys.setenv(MKL_NUM_THREADS      = "1")
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  }
  
  if (use_parallel) {
    try(unregister_dopar(), silent = TRUE)
    cl <- tryCatch({
      if (.Platform$OS.type == "unix") {
        parallel::makeForkCluster(nnodes = effective_ncores)
      } else {
        parallel::makePSOCKcluster(effective_ncores)
      }
    }, error = function(e) {
      message("[parallel] Could not create cluster (", conditionMessage(e),
              "); running sequentially.")
      NULL
    })
    if (!is.null(cl)) {
      doParallel::registerDoParallel(cl)
      # Ensure in PSOCK workers as well (fork inherits env-vars from the parent).
      if (.Platform$OS.type != "unix") {
        try(parallel::clusterEvalQ(cl, {
          Sys.setenv(OMP_NUM_THREADS = "1",
                     MKL_NUM_THREADS = "1",
                     OPENBLAS_NUM_THREADS = "1")
          if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
            try(RhpcBLASctl::omp_set_num_threads(1L),  silent = TRUE)
          }
        }), silent = TRUE)
      }
    }
  }
  
  # Explicit log: leave no doubt about whether parallelization is active.
  if (use_parallel && !is.null(cl)) {
    message(sprintf("[PARALLEL][ref] Cluster created: %d workers (%s) for %d markers.",
                    effective_ncores,
                    if (.Platform$OS.type == "unix") "fork" else "PSOCK",
                    length(adt_marker_index_list)))
  } else {
    message(sprintf("[PARALLEL][ref] SEQUENTIAL execution (ncores=%s, use_parallel=%s).",
                    as.character(ncores), as.character(use_parallel)))
  }
  
  on.exit({
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
      try(foreach::registerDoSEQ(), silent = TRUE)
      try(unregister_dopar(), silent = TRUE)
    }
    # Restaurar hilos BLAS/OMP/MKL previos.
    if (!is.na(.prev_blas_threads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(.prev_blas_threads), silent = TRUE)
    }
    if (is.na(.prev_omp_threads))  Sys.unsetenv("OMP_NUM_THREADS")      else Sys.setenv(OMP_NUM_THREADS      = .prev_omp_threads)
    if (is.na(.prev_mkl_threads))  Sys.unsetenv("MKL_NUM_THREADS")      else Sys.setenv(MKL_NUM_THREADS      = .prev_mkl_threads)
    if (is.na(.prev_openblas_t))   Sys.unsetenv("OPENBLAS_NUM_THREADS") else Sys.setenv(OPENBLAS_NUM_THREADS = .prev_openblas_t)
  }, add = TRUE)
  
  .pkgs_needed <- c("Matrix", "stats", "utils", "ADTnorm")
  .exports_needed <- c(
    ".run_one_marker_ref", ".safe_to_matrix", ".as_dense_matrix",
    ".build_target_landmark", ".extract_align_field",
    ".sqs_compute_per_sample", ".sqs_one_peak", ".sqs_two_peak",
    ".sqs_multi_peak", ".sqs_valley_deep_scaler", ".sqs_area_under_curve",
    ".get_peak_mode_patched"
  )

  `%OP%` <- if (use_parallel && !is.null(cl)) `%dopar%` else `%do%`
  
  # -------------------------------------------------------------------------
  # REDUCED PAYLOAD PER WORKER (KEY for large matrices)
  # ---------------------------------------------------------------------
  # In the previous design each worker received the full cell_x_adt and
  # cell_x_feature. With fork this uses COW but, in R, `cell_x_adt[, j] <- NA`
  # inside the worker triggers copy-on-modify at object level -> each
  # worker duplicated the entire matrix in RAM (oversubscription -> thrashing).
  # Here we pre-extract ONLY the marker column + a slim version of the
  # feature-frame (sample, batch). The helper remains functional because
  # its internal calls (get_peak_mode, get_valley_location, etc.)
  # operate on cell_x_adt_sub[, adt_marker_select], which is exactly
  # that column.
  # -------------------------------------------------------------------------
  sample_vec_full <- cell_x_feature$sample
  feature_slim    <- cell_x_feature[, intersect(c("sample", "batch"),
                                                colnames(cell_x_feature)),
                                    drop = FALSE]
  # Build the per-marker data list up front (ONLY ONCE):
  # with fork this is shared via COW; with PSOCK only the minimum is serialized.
  marker_payload <- lapply(adt_marker_index_list, function(k) {
    tmp_name <- colnames(cell_x_adt)[k]
    list(
      adt_marker_index = k,
      tmp_name         = tmp_name,
      marker_col       = cell_x_adt[, tmp_name, drop = TRUE]
    )
  })
  
  t_foreach <- proc.time()
  marker_results <- foreach::foreach(
    payload          = marker_payload,
    .packages        = .pkgs_needed,
    .export          = .exports_needed,
    .errorhandling   = "pass"
  ) %OP% {
    # Reconstruct a ONE-column matrix with the same rowname that
    # get_peak_mode / get_valley_location / peak_alignment expect.
    local_mat <- matrix(payload$marker_col, ncol = 1L,
                        dimnames = list(names(payload$marker_col), payload$tmp_name))
    
    .run_one_marker_ref(
      adt_marker_index         = 1L,                     # always 1 in the submatrix
      cell_x_adt               = local_mat,
      cell_x_feature           = feature_slim,
      original_marker_name     = original_marker_name[payload$adt_marker_index],
      quantile_clip            = quantile_clip,
      bw_smallest_bi           = bw_smallest_bi,
      bw_smallest_tri          = bw_smallest_tri,
      bw_smallest_adjustments  = bw_smallest_adjustments,
      bimodal_marker_index     = if (payload$adt_marker_index %in% bimodal_marker_index) 1L else integer(0),
      trimodal_marker_index    = if (payload$adt_marker_index %in% trimodal_marker_index) 1L else integer(0),
      positive_peak            = positive_peak,
      peak_type                = peak_type,
      neg_candidate_thres      = neg_candidate_thres,
      lower_peak_thres         = lower_peak_thres,
      arcsinh_transform_flag   = arcsinh_transform_flag,
      shoulder_valley          = shoulder_valley,
      multi_sample_per_batch   = multi_sample_per_batch,
      valley_density_adjust    = valley_density_adjust,
      shoulder_valley_slope    = shoulder_valley_slope,
      detect_outlier_valley    = detect_outlier_valley,
      sqs_threshold            = sqs_threshold,
      apply_sqs                = apply_sqs,
      file.id                  = file.id,
      landmark_align_type      = landmark_align_type,
      midpoint_type            = midpoint_type,
      target_landmark_location = target_landmark_location
    )
  }
  message(sprintf("[TIMING][ref] foreach elapsed: %.2fs",
                  (proc.time() - t_foreach)[["elapsed"]]))
  
  # Reassign the real tmp_name (the helper returned it as the submatrix's name).
  for (i in seq_along(marker_results)) {
    if (!inherits(marker_results[[i]], "error") &&
        !is.null(marker_results[[i]]) &&
        !is.null(marker_results[[i]]$tmp_name)) {
      marker_results[[i]]$tmp_name <- marker_payload[[i]]$tmp_name
    }
  }
  
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
    try(unregister_dopar(), silent = TRUE)
    cl <- NULL
  }
  
  for (res in marker_results) {
    if (inherits(res, "error") || is.null(res) || is.null(res$tmp_name)) {
      msg <- if (inherits(res, "error")) conditionMessage(res) else "invalid result"
      message("[parallel] A worker failed: ", msg, " (skipping the marker).")
      next
    }
    if (!is.null(res$norm_col_override)) {
      cell_x_adt_norm[, res$tmp_name] <- res$norm_col_override
    }
    if (!is.null(res$pars)) {
      bg.prot.pars$markers[[res$marker_name]] <- res$pars
    }
    for (row in res$sqs_rows) .sqs_log_append(sqs_log_path, row)
  }
  
  colnames(cell_x_adt_norm) <- original_marker_name
  norm_fxC <- t(as.matrix(cell_x_adt_norm))
  
  if (do_subsample) {
    if (verbose) message("Applying hierarchical subsampling...")
    set.seed(seed)
    
    if (!"cell_id" %in% colnames(metadata)) {
      metadata$cell_id <- rownames(metadata)
    }
    
    all_na_cells <- rowSums(!is.na(cell_x_adt_norm)) == 0
    if (any(all_na_cells)) {
      drop_ids <- rownames(cell_x_adt_norm)[all_na_cells]
      if (verbose) message("SQS: excluding ", length(drop_ids),
                           " cells with 100% of markers in NA from the subsampling pool.")
      metadata <- metadata[!(metadata$cell_id %in% drop_ids), , drop = FALSE]
    }
    
    alloc_bal <- function(tgt, caps) {
      caps <- caps[caps > 0]; if(length(caps) == 0) return(NULL)
      if(sum(caps) < tgt) return(caps)
      al <- setNames(integer(length(caps)), names(caps))
      while(tgt > 0) {
        act <- names(caps)[caps > al]; act <- act[order(al[act])]
        tk <- act[seq_len(min(tgt, length(act)))]; al[tk] <- al[tk] + 1L; tgt <- tgt - length(tk)
      }
      return(al)
    }
    
    selected_ids <- character()
    
    for(c in unique(metadata[[ct_col]])) {
      if (is.na(c)) next
      df_c <- metadata[metadata[[ct_col]] == c, ]
      
      t_al <- alloc_bal(target_ct, table(df_c[[tissue_col]]))
      if(is.null(t_al)) next
      
      for(tn in names(t_al)[t_al > 0]) {
        df_t <- df_c[df_c[[tissue_col]] == tn, ]
        
        d_al <- alloc_bal(t_al[[tn]], table(df_t[[dataset_col]]))
        if(is.null(d_al)) next
        
        for(dn in names(d_al)[d_al > 0]) {
          ids <- df_t$cell_id[df_t[[dataset_col]] == dn]
          selected_ids <- c(selected_ids, sample(ids, d_al[[dn]], replace = FALSE))
        }
      }
    }
    
    valid_ids <- intersect(selected_ids, colnames(norm_fxC))
    
    if (length(valid_ids) == 0) {
      warning("Subsampling returned 0 cells or the IDs do not intersect with norm_fxC. The original matrix is returned.")
    } else {
      norm_fxC <- norm_fxC[, valid_ids, drop = FALSE]
      if (verbose) message("Subsampling finished. Cells retained: ", ncol(norm_fxC))
    }
  }
  
  if ("cell_id" %in% colnames(metadata)) {
    metadata <- metadata[metadata$cell_id %in% colnames(norm_fxC), , drop = FALSE]
  } else {
    metadata <- metadata[rownames(metadata) %in% colnames(norm_fxC), , drop = FALSE]
  }
  
  if (saveData) {
    saveRDS(norm_fxC, file = file.path(outdir, paste0(file.id, "_NormData.rds")))
    saveRDS(bg.prot.pars, file = file.path(outdir, paste0(file.id, "_bg_prot_pars.rds")))
    saveRDS(metadata, file = file.path(outdir, paste0(file.id, "_Metadata.rds")))
  }
  
  norm_fxC
}
# -----------------------------------------------------------------------------
# .run_one_marker_query  [NUEVO]
# -----------------------------------------------------------------------------
.run_one_marker_query <- function(marker_name, adt_marker_select, adt_marker_index,
                                  cell_x_adt, cell_x_feature, marker_pars, gs, mode_info,
                                  file.id, sqs_threshold, apply_sqs) {
  
  out_skip <- list(
    marker_name       = marker_name,
    tmp_name          = adt_marker_select,
    norm_col_override = NULL,  # Keep the philosophy: NULL = ignore in master
    sqs_rows          = list()
  )
  
  marker_vec <- cell_x_adt[, adt_marker_select]
  if (gs$quantile_clip < 1) {
    quant <- stats::quantile(marker_vec, gs$quantile_clip, na.rm = TRUE)
    clip_idx <- !is.na(marker_vec) & (marker_vec > quant)
    marker_vec[clip_idx] <- NA
  }
  
  valid_samples <- c()
  for (s in levels(cell_x_feature$sample)) {
    vals <- marker_vec[cell_x_feature$sample == s]
    valid_vals <- vals[is.finite(vals)]
    if (length(valid_vals) >= 5 && length(unique(valid_vals)) >= 3 && stats::var(valid_vals) > 1e-6) {
      valid_samples <- c(valid_samples, s)
    }
  }
  
  # If there are no valid samples, return the clipped vector and skip (as in the sequential branch)
  if (length(valid_samples) == 0) {
    return(list(marker_name = marker_name, tmp_name = adt_marker_select,
                norm_col_override = marker_vec, sqs_rows = list()))
  }
  
  valid_idx <- cell_x_feature$sample %in% valid_samples
  cell_x_adt_sub <- cell_x_adt[valid_idx, , drop = FALSE]
  cell_x_adt_sub[, adt_marker_select] <- marker_vec[valid_idx]
  cell_x_feature_sub <- cell_x_feature[valid_idx, , drop = FALSE]
  cell_x_feature_sub$sample <- factor(cell_x_feature_sub$sample, levels = valid_samples)
  
  peak_mode_res <- tryCatch({
    if (marker_pars$peak_type == "mode") {
      .get_peak_mode_patched(cell_x_adt_sub, cell_x_feature_sub, adt_marker_select, adt_marker_index, marker_pars$bwFac_smallest, gs$bimodal_marker_index, gs$trimodal_marker_index, marker_pars$positive_peak, marker_pars$neg_candidate_thres, marker_pars$lower_peak_thres, mode_info$arcsinh_transform_flag)
    } else {
      ADTnorm::get_peak_midpoint(cell_x_adt_sub, cell_x_feature_sub, adt_marker_select, adt_marker_index, marker_pars$bwFac_smallest, gs$bimodal_marker_index, gs$trimodal_marker_index, marker_pars$positive_peak, marker_pars$neg_candidate_thres, marker_pars$lower_peak_thres, mode_info$arcsinh_transform_flag)
    }
  }, error = function(e) NULL)
  
  if (is.null(peak_mode_res)) return(out_skip)
  
  peak_mode_res_raw <- peak_mode_res
  ref_n <- length(marker_pars$fixed_target_landmark)
  peak_mode_res <- .coerce_peak_mode_to_reference(peak_mode_res, ref_n)
  
  peak_valley_list <- tryCatch(
    ADTnorm::get_valley_location(cell_x_adt_sub, cell_x_feature_sub, adt_marker_select, peak_mode_res, marker_pars$shoulder_valley, marker_pars$positive_peak, marker_pars$multi_sample_per_batch, marker_pars$valley_density_adjust, 20, marker_pars$shoulder_valley_slope, marker_pars$neg_candidate_thres, mode_info$arcsinh_transform_flag),
    error = function(e) NULL
  )
  
  if (is.null(peak_valley_list)) return(out_skip)
  
  valley_location_res_raw   <- peak_valley_list$valley_location_list
  peak_landmark_list_raw    <- peak_valley_list$peak_landmark_list
  valley_location_res <- peak_valley_list$valley_location_list
  peak_landmark_list <- .coerce_peak_mode_to_reference(peak_valley_list$peak_landmark_list, ref_n)
  
  if (marker_pars$detect_outlier_valley != FALSE) {
    valley_location_res <- tryCatch(
      ADTnorm::detect_impute_outlier_valley(valley_location_res, adt_marker_select, cell_x_adt_sub, cell_x_feature_sub, 3, marker_pars$detect_outlier_valley, 3, 0.75),
      error = function(e) NULL
    )
  }
  
  sqs_rows <- list()
  norm_col_override <- cell_x_adt[, adt_marker_select] 
  
  peak_info_raw_mat   <- .safe_to_matrix(peak_landmark_list_raw,  valid_samples)
  valley_info_raw_mat <- .safe_to_matrix(valley_location_res_raw, valid_samples)
  
  sqs_scores_q <- .sqs_compute_per_sample(
    marker_vec    = cell_x_adt_sub[, adt_marker_select],
    sample_vec    = as.character(cell_x_feature_sub$sample),
    valid_samples = valid_samples,
    peak_info     = peak_info_raw_mat,
    valley_info   = valley_info_raw_mat
  )
  
  for (s in valid_samples) {
    sqs_rows[[length(sqs_rows) + 1L]] <- data.frame(
      stage      = "query",
      file.id    = file.id,
      marker     = marker_name,
      sample     = s,
      sqs        = unname(sqs_scores_q[s]),
      threshold  = sqs_threshold,
      apply_sqs  = apply_sqs,
      action     = if (apply_sqs == "both" && sqs_scores_q[s] < sqs_threshold) "NA" else "keep",
      stringsAsFactors = FALSE
    )
  }
  
  if (apply_sqs == "both") {
    fail_samples_q <- names(sqs_scores_q)[sqs_scores_q < sqs_threshold]
    if (length(fail_samples_q) > 0) {
      fail_idx_q <- cell_x_feature$sample %in% fail_samples_q
      norm_col_override[fail_idx_q] <- NA_real_
      
      keep_mask_q     <- !(valid_samples %in% fail_samples_q)
      valid_samples   <- valid_samples[keep_mask_q]
      peak_landmark_list  <- peak_landmark_list[keep_mask_q, , drop = FALSE]
      if (!is.null(valley_location_res)) {
        valley_location_res <- .safe_to_matrix(valley_location_res, rownames(peak_info_raw_mat))
        valley_location_res <- valley_location_res[keep_mask_q, , drop = FALSE]
      }
      
      if (length(valid_samples) == 0) {
        return(list(marker_name = marker_name, tmp_name = adt_marker_select, norm_col_override = norm_col_override, sqs_rows = sqs_rows))
      }
      
      valid_idx <- cell_x_feature$sample %in% valid_samples
      cell_x_adt_sub     <- cell_x_adt[valid_idx, , drop = FALSE]
      cell_x_adt_sub[, adt_marker_select] <- marker_vec[valid_idx]
      cell_x_feature_sub <- cell_x_feature[valid_idx, , drop = FALSE]
      cell_x_feature_sub$sample <- factor(cell_x_feature_sub$sample, levels = valid_samples)
    }
  }
  
  full_samples <- levels(cell_x_feature$sample)
  peak_landmark_list <- .safe_to_matrix(peak_landmark_list, valid_samples)
  pm_full <- matrix(NA_real_, nrow = length(full_samples), ncol = ncol(peak_landmark_list))
  rownames(pm_full) <- full_samples
  pm_full[valid_samples, ] <- peak_landmark_list
  peak_mode_res_final <- pm_full
  
  valley_location_res <- .safe_to_matrix(valley_location_res, valid_samples)
  vl_full <- matrix(NA_real_, nrow = length(full_samples), ncol = ncol(valley_location_res))
  rownames(vl_full) <- full_samples
  vl_full[valid_samples, ] <- valley_location_res
  valley_location_res_final <- vl_full
  
  landmark_matrix <- tryCatch(
    ADTnorm::landmark_fill_na(peak_mode_res_final, valley_location_res_final, gs$landmark_align_type, marker_pars$midpoint_type, marker_pars$neg_candidate_thres),
    error = function(e) NULL
  )
  
  if (is.null(landmark_matrix)) {
    return(list(marker_name = marker_name, tmp_name = adt_marker_select, norm_col_override = norm_col_override, sqs_rows = sqs_rows))
  }
  
  fixed_target_landmark <- marker_pars$fixed_target_landmark
  landmark_matrix <- .coerce_landmark_matrix_to_reference(landmark_matrix, fixed_target_landmark, marker_pars$neg_candidate_thres)
  
  align_res <- tryCatch(
    ADTnorm::peak_alignment(cell_x_adt_sub[, adt_marker_select], cell_x_feature_sub, landmark_matrix[valid_samples, , drop=FALSE], fixed_target_landmark, "fixed", marker_pars$warp_parameters),
    error = function(e) NULL
  )
  
  if (!is.null(align_res)) {
    norm_col_override[valid_idx] <- .extract_align_field(align_res, "cell_x_adt_norm", pos = 1)
  }
  
  list(marker_name = marker_name, tmp_name = adt_marker_select, norm_col_override = norm_col_override, sqs_rows = sqs_rows)
}

# -----------------------------------------------------------------------------
# QUERY MODE
# -----------------------------------------------------------------------------
.ADTprot_norm_query <- function(data_counts, metadata, bg.prot.pars, file.id, outdir,
                                sample_col, batch_col, marker_to_process, verbose,
                                saveData, is_transformed,
                                ## Defaults aligned across the chain so SQS is
                                ## off-by-default. Callers downstream of run_scProteinXpress
                                ## still inherit gs$apply_sqs from bg.prot.pars when
                                ## available (see ADTprot_norm wrapper).
                                sqs_threshold = 1.0, apply_sqs = "none",
                                ncores = 1L){
  ## Validate bg.prot.pars before touching anything else. An empty or
  ## malformed reference fitted-parameters list would otherwise surface as a
  ## confusing "subscript out of bounds" deep inside the foreach loop.
  if (is.null(bg.prot.pars) || !is.list(bg.prot.pars) ||
      is.null(bg.prot.pars$markers) || length(bg.prot.pars$markers) == 0L) {
    stop(".ADTprot_norm_query: bg.prot.pars must contain a non-empty `markers` list ",
         "(loaded from <Ref>/<file.name>_bg_prot_pars.rds).")
  }
  if (is.null(bg.prot.pars$global_settings)) {
    stop(".ADTprot_norm_query: bg.prot.pars$global_settings is missing.")
  }

  data_counts <- .as_dense_matrix(data_counts)
  metadata <- as.data.frame(metadata)

  ref_markers <- names(bg.prot.pars$markers)
  ## Explicit Ref/Query marker intersection check. .pad_to_features pads
  ## with NA, so the matrix would survive even with zero overlap. Detect that
  ## case explicitly because it is almost always a metadata/panel mismatch.
  query_markers <- rownames(data_counts)
  if (length(intersect(ref_markers, query_markers)) == 0L) {
    stop(sprintf(
      ".ADTprot_norm_query: no overlap between query markers (%d) and Ref markers (%d). ",
      length(query_markers), length(ref_markers)),
      "Examples query: ", paste(utils::head(query_markers, 5), collapse = ", "),
      "; examples Ref: ", paste(utils::head(ref_markers, 5), collapse = ", "))
  }
  data_counts <- .pad_to_features(data_counts, ref_markers)
  
  adt_in <- make_ADTnorm_input(data_counts, metadata, sample_col, batch_col)
  cell_x_adt <- adt_in$cell_x_adt
  cell_x_feature <- adt_in$cell_x_feature
  gs <- bg.prot.pars$global_settings
  
  mode_info <- .detect_input_mode(cell_x_adt, gs$neg_candidate_thres, is_transformed, verbose)
  if (mode_info$arcsinh_transform_flag) {
    if(verbose) message("Applying arcsinh transformation to query...")
    cell_x_adt <- ADTnorm::arcsinh_transform(cell_x_adt = cell_x_adt)
  }
  
  if (!is.factor(cell_x_feature$sample)) cell_x_feature$sample <- factor(cell_x_feature$sample, levels = unique(cell_x_feature$sample))
  
  original_marker_name <- colnames(cell_x_adt)
  colnames(cell_x_adt) <- paste0("tmpName", seq_len(ncol(cell_x_adt)))
  name_map <- stats::setNames(colnames(cell_x_adt), original_marker_name)
  
  markers_to_run <- ref_markers
  cell_x_adt_norm <- cell_x_adt[, name_map[markers_to_run], drop = FALSE]
  
  sqs_log_path_q <- file.path(outdir, paste0(file.id, "_SQS_query_log.csv"))
  
  effective_ncores <- if (is.numeric(ncores) && length(ncores) == 1L)
    min(as.integer(ncores), length(markers_to_run)) else 1L
  
  cl <- NULL
  use_parallel <- effective_ncores >= 2L && length(markers_to_run) >= 2L
  
  # --- Pin BLAS to 1 thread per worker: avoids oversubscription ---------------
  .prev_blas_threads <- NA_integer_
  .prev_omp_threads  <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  .prev_mkl_threads  <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
  .prev_openblas_t   <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  if (use_parallel) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      .prev_blas_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
                                     error = function(e) NA_integer_)
      try(RhpcBLASctl::blas_set_num_threads(1L),  silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(1L),   silent = TRUE)
    }
    Sys.setenv(OMP_NUM_THREADS      = "1")
    Sys.setenv(MKL_NUM_THREADS      = "1")
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  }
  
  if (use_parallel) {
    try(unregister_dopar(), silent = TRUE)
    cl <- tryCatch({
      if (.Platform$OS.type == "unix") {
        parallel::makeForkCluster(nnodes = effective_ncores)
      } else {
        parallel::makePSOCKcluster(effective_ncores)
      }
    }, error = function(e) {
      message("[parallel] Could not create cluster (", conditionMessage(e), "); running sequentially.")
      NULL
    })
    if (!is.null(cl)) {
      doParallel::registerDoParallel(cl)
      if (.Platform$OS.type != "unix") {
        try(parallel::clusterEvalQ(cl, {
          Sys.setenv(OMP_NUM_THREADS = "1",
                     MKL_NUM_THREADS = "1",
                     OPENBLAS_NUM_THREADS = "1")
          if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
            try(RhpcBLASctl::omp_set_num_threads(1L),  silent = TRUE)
          }
        }), silent = TRUE)
      }
    }
  }
  
  if (use_parallel && !is.null(cl)) {
    message(sprintf("[PARALLEL][query] Cluster created: %d workers (%s) for %d markers.",
                    effective_ncores,
                    if (.Platform$OS.type == "unix") "fork" else "PSOCK",
                    length(markers_to_run)))
  } else {
    message(sprintf("[PARALLEL][query] SEQUENTIAL execution (ncores=%s, use_parallel=%s).",
                    as.character(ncores), as.character(use_parallel)))
  }
  
  on.exit({
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
      try(foreach::registerDoSEQ(), silent = TRUE)
      try(unregister_dopar(), silent = TRUE)
    }
    if (!is.na(.prev_blas_threads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(.prev_blas_threads), silent = TRUE)
    }
    if (is.na(.prev_omp_threads))  Sys.unsetenv("OMP_NUM_THREADS")      else Sys.setenv(OMP_NUM_THREADS      = .prev_omp_threads)
    if (is.na(.prev_mkl_threads))  Sys.unsetenv("MKL_NUM_THREADS")      else Sys.setenv(MKL_NUM_THREADS      = .prev_mkl_threads)
    if (is.na(.prev_openblas_t))   Sys.unsetenv("OPENBLAS_NUM_THREADS") else Sys.setenv(OPENBLAS_NUM_THREADS = .prev_openblas_t)
  }, add = TRUE)
  
  .pkgs_needed <- c("Matrix", "stats", "utils", "ADTnorm")
  .exports_needed <- c(
    ".run_one_marker_query", ".safe_to_matrix", ".as_dense_matrix",
    ".build_target_landmark", ".extract_align_field",
    ".sqs_compute_per_sample", ".sqs_one_peak", ".sqs_two_peak",
    ".sqs_multi_peak", ".sqs_valley_deep_scaler", ".sqs_area_under_curve",
    ".coerce_peak_mode_to_reference", ".coerce_landmark_matrix_to_reference",
    ".get_peak_mode_patched"
  )
  
  `%OP%` <- if (use_parallel && !is.null(cl)) `%dopar%` else `%do%`
  
  t_foreach_q <- proc.time()
  marker_results_q <- foreach::foreach(
    marker_name      = markers_to_run,
    .packages        = .pkgs_needed,
    .export          = .exports_needed,
    .errorhandling   = "pass"
  ) %OP% {
    adt_marker_select <- name_map[[marker_name]]
    adt_marker_index <- which(colnames(cell_x_adt) == adt_marker_select)
    marker_pars <- bg.prot.pars$markers[[marker_name]]
    
    .run_one_marker_query(
      marker_name = marker_name,
      adt_marker_select = adt_marker_select,
      adt_marker_index = adt_marker_index,
      cell_x_adt = cell_x_adt,
      cell_x_feature = cell_x_feature,
      marker_pars = marker_pars,
      gs = gs,
      mode_info = mode_info,
      file.id = file.id,
      sqs_threshold = sqs_threshold,
      apply_sqs = apply_sqs
    )
  }
  message(sprintf("[TIMING][query] foreach elapsed: %.2fs",
                  (proc.time() - t_foreach_q)[["elapsed"]]))
  
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
    try(unregister_dopar(), silent = TRUE)
    cl <- NULL
  }
  
  for (res in marker_results_q) {
    if (inherits(res, "error") || is.null(res) || is.null(res$tmp_name)) {
      msg <- if (inherits(res, "error")) conditionMessage(res) else "invalid result"
      message("[parallel] A worker failed: ", msg, " (skipping the marker in query).")
      next
    }
    if (!is.null(res$norm_col_override)) {
      cell_x_adt_norm[, res$tmp_name] <- res$norm_col_override
    }
    for (row in res$sqs_rows) .sqs_log_append(sqs_log_path_q, row)
  }
  
  colnames(cell_x_adt_norm) <- markers_to_run
  norm_fxC <- t(as.matrix(cell_x_adt_norm))
  
  if (saveData) saveRDS(norm_fxC, file = file.path(outdir, paste0(file.id, "_NormData.rds")))
  norm_fxC
}
###############################################################################
# Public wrapper
###############################################################################
#' Reference / query ADT normalisation
#'
#' Dispatcher around `.ADTprot_norm_ref()` / `.ADTprot_norm_query()` that
#' enforces SLURM-aware parallelisation and the documented QC filters when
#' invoked through `run_scProteinXpress()`.
#'
#' @param data_counts proteins x cells matrix or Matrix.
#' @param metadata data.frame; must contain `cell.id`, `sample`, `batch`.
#' @param data.type "Ref" or "query".
#' @param file.id character(1) used as filename prefix.
#' @param bg.prot.pars list returned by Ref; ignored for Ref builds.
#' @param saveData logical; whether to save the normalised matrix to disk.
#' @param return logical; whether to return the normalised matrix.
#' @param outdir character(1); output directory for saved files.
#' @param sample_col character(1); name of the metadata column holding sample ids.
#' @param batch_col character(1); name of the metadata column holding batch ids.
#' @param marker_to_process optional character vector; restrict normalisation to these markers.
#' @param exclude_zeroes logical; ADTnorm parameter; passed through.
#' @param bimodal_marker optional character vector; markers known to be bimodal.
#' @param trimodal_marker optional character vector; markers known to be trimodal.
#' @param positive_peak optional list; per-marker positive peak hints.
#' @param bw_smallest_bi numeric; minimum bandwidth for bimodal markers.
#' @param bw_smallest_tri numeric; minimum bandwidth for trimodal markers.
#' @param bw_smallest_adjustments named list; per-marker bandwidth adjustments.
#' @param quantile_clip numeric in the range 0 to 1; quantile used to clip extreme values.
#' @param peak_type character; ADTnorm peak detection mode ("mode" or "midpoint").
#' @param multi_sample_per_batch logical; whether multiple samples per batch are expected.
#' @param shoulder_valley logical; enable shoulder-valley detection.
#' @param shoulder_valley_slope numeric; slope threshold for shoulder detection.
#' @param valley_density_adjust numeric; KDE bandwidth adjust for valley detection.
#' @param landmark_align_type character; landmark alignment mode (e.g. "negPeak_valley_posPeak").
#' @param midpoint_type character; midpoint definition (typically "valley").
#' @param neg_candidate_thres optional numeric; threshold for the negative-peak candidate.
#' @param lower_peak_thres numeric; lower bound for peak height to be considered.
#' @param detect_outlier_valley logical; enable outlier-valley imputation.
#' @param target_landmark_location optional named list; user-provided landmark targets.
#' @param clean_adt_name logical; whether to sanitise ADT marker names.
#' @param override_landmark optional named list; user overrides for specific landmarks.
#' @param verbose logical; print progress messages.
#' @param is_transformed logical or NULL; if NULL the input mode is auto-detected.
#' @param sqs_threshold numeric; Stain Quality Score threshold (only used when apply_sqs is enabled).
#' @param apply_sqs character; one of "none" (default), "ref_only", "both".
#' @param do_subsample logical; apply hierarchical subsampling on the Ref branch.
#' @param target_ct integer; cap of cells per cell type during Ref subsampling.
#' @param ncores integer; number of cores for parallelisation.
#' @return Normalised proteins x cells matrix.
#' @export
ADTprot_norm <- function(data_counts,
                         metadata, 
                         data.type = c("Ref", "query"), 
                         file.id = "mydata", 
                         bg.prot.pars = NULL, 
                         saveData = TRUE, 
                         return = TRUE, 
                         outdir = ".", 
                         sample_col = "sample", 
                         batch_col = "batch", 
                         marker_to_process = NULL, 
                         exclude_zeroes = FALSE, 
                         bimodal_marker = NULL, 
                         trimodal_marker = NULL, 
                         positive_peak = NULL, 
                         bw_smallest_bi = 1.1, 
                         bw_smallest_tri = 0.8, 
                         bw_smallest_adjustments =  list(CD3 = 0.8, CD4 = 0.8, CD8A = 0.8, CD8B = 0.8), 
                         quantile_clip = 1, 
                         peak_type = "mode", 
                         multi_sample_per_batch = FALSE, 
                         shoulder_valley = TRUE, 
                         shoulder_valley_slope = -0.5, 
                         valley_density_adjust = 3, 
                         landmark_align_type = "negPeak_valley_posPeak", 
                         midpoint_type = "valley", 
                         neg_candidate_thres = NULL, 
                         lower_peak_thres = 0.001, 
                         detect_outlier_valley = FALSE, 
                         target_landmark_location = NULL, 
                         clean_adt_name = FALSE, 
                         override_landmark = NULL, 
                         verbose = TRUE, 
                         is_transformed = NULL,
                         sqs_threshold = 1.0,
                         apply_sqs = "none",
                         # --- Hierarchical subsampling (only applies to the Ref branch) ---
                         # do_subsample must remain TRUE to guarantee that
                         # the background is always built with a balanced number
                         # of cells per (cell_type x tissue x dataset). target_ct
                         # is the cap per cell_type and is the only parameter
                         # exposed to configure the filter size.
                         do_subsample = TRUE,
                         target_ct = 2000L,
                         ncores = 1L){
  data.type <- match.arg(data.type)
  apply_sqs <- match.arg(apply_sqs, choices = c("ref_only", "both", "none"))
  if (data.type == "Ref") {
    res <- .ADTprot_norm_ref(data_counts, metadata, file.id, outdir, sample_col, 
                             batch_col, marker_to_process, exclude_zeroes, 
                             bimodal_marker, trimodal_marker, positive_peak, 
                             bw_smallest_bi, bw_smallest_tri, bw_smallest_adjustments, 
                             quantile_clip, peak_type, multi_sample_per_batch, 
                             shoulder_valley, shoulder_valley_slope, valley_density_adjust, 
                             landmark_align_type, midpoint_type, neg_candidate_thres, 
                             lower_peak_thres, detect_outlier_valley, target_landmark_location, 
                             clean_adt_name, override_landmark, verbose, saveData, is_transformed,
                             sqs_threshold = sqs_threshold, apply_sqs = apply_sqs,
                             do_subsample = do_subsample, target_ct = target_ct,
                             ncores = ncores)
  } else {
    if (!is.null(bg.prot.pars) && !is.null(bg.prot.pars$global_settings)) {
      gs <- bg.prot.pars$global_settings
      if (missing(sqs_threshold) && !is.null(gs$sqs_threshold)) sqs_threshold <- gs$sqs_threshold
      if (missing(apply_sqs)     && !is.null(gs$apply_sqs))     apply_sqs     <- gs$apply_sqs
    }
    res <- .ADTprot_norm_query(data_counts, metadata, bg.prot.pars, file.id, outdir, sample_col, 
                               batch_col, marker_to_process, verbose, saveData, is_transformed,
                               sqs_threshold = sqs_threshold, apply_sqs = apply_sqs,
                               ncores = ncores)
  }
  if (return) return(res)
  invisible(NULL)
}
###############################################################################
# Query prep for scProteinXpress
# ----------------------------------------------------------------------------
# Combines:
#   - Global QC with audit CSV
#   - ncores propagation to ADTprot_norm
#   - verbose = TRUE for ADTprot_norm calls
###############################################################################
#' Prepare a query against an existing reference background
#'
#' Runs the global QC, ADTnorm query alignment to the Ref landmarks, and the
#' compatibility scaleMaxGexp() step.
#'
#' @param query_dat proteins x cells matrix.
#' @param query_metadata data.frame.
#' @param Ref_data Path to the reference directory (containing `_sf.txt`
#'   and `_bg_prot_pars.rds`) or an in-memory list.
#' @param dir Output directory.
#' @param sfback Optional pre-loaded scaling-factor data.frame.
#' @param file.id character(1) used as filename prefix.
#' @param bg.prot.pars Pre-loaded reference normalisation parameters.
#' @param normalise,scale_data,saveData,return.data logical switches.
#' @param is_transformed NULL / TRUE / FALSE; see `run_scProteinXpress()`.
#' @param sqs_threshold,apply_sqs SQS controls; see `run_scProteinXpress()`.
#' @param min.cells.feature,min.coverage.fraction Marker-level coverage filter.
#' @param qc_outlier_cells "no" / "yes"; resolved upstream.
#' @param qc_nmads MAD threshold for scuttle::isOutlier.
#' @param ncores integer >= 1.
#' @return Scaled proteins x cells matrix.
#' @export
PrepQueryData_ADT <- function(query_dat, query_metadata, Ref_data, dir, sfback = NA,
                              file.id = "query", bg.prot.pars, normalise = TRUE,
                              scale_data = TRUE, saveData = TRUE, return.data = TRUE,
                              is_transformed = NULL,
                              ## Aligned with the rest of the chain (none = off).
                              sqs_threshold = 1.0, apply_sqs = "none",
                              # Global coverage filters
                              min.cells.feature     = 10,
                              min.coverage.fraction = 0.10,
                              # Opt-in cell outlier filter, resolved upstream by
                              # run_scProteinXpress -- this function expects the
                              # concrete "yes" / "no" string here.
                              qc_outlier_cells = c("no", "yes"),
                              qc_nmads = 3,
                              # Parallelisation forwarded to ADTprot_norm
                              ncores = 1L){
  qc_outlier_cells <- match.arg(qc_outlier_cells)
  dir_name <- dir
  if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE)

  log.file <- file.path(dir_name, paste0(file.id, "_log.txt"))
  cat("Preparing ADT data: QC + ADTnorm query normalization + compatibility scaling\n",
      file = log.file, append = TRUE)

  # Single global QC with audit CSV (no per-dataset partial filters here).
  qc <- cleanData_ADT(query_dat, nmads = qc_nmads,
                      min.cells.feature     = min.cells.feature,
                      min.coverage.fraction = min.coverage.fraction,
                      audit_csv_path        = file.path(dir_name,
                                                        paste0(file.id, "_excluded_markers.csv")),
                      qc_outlier_cells      = qc_outlier_cells)
  x_counts <- qc$Counts

  query_metadata <- as.data.frame(query_metadata)
  query_metadata <- query_metadata[query_metadata$cell.id %in% colnames(x_counts), , drop = FALSE]
  rownames(query_metadata) <- query_metadata$cell.id
  query_metadata <- query_metadata[colnames(x_counts), , drop = FALSE]

  if (normalise) {
    x_norm <- ADTprot_norm(data_counts    = x_counts,
                           metadata       = query_metadata,
                           data.type      = "query",
                           file.id        = file.id,
                           bg.prot.pars   = bg.prot.pars,
                           saveData       = saveData,
                           return         = TRUE,
                           outdir         = dir_name,
                           sample_col     = "sample",
                           batch_col      = "batch",
                           verbose        = TRUE,
                           is_transformed = is_transformed,
                           sqs_threshold  = sqs_threshold,
                           apply_sqs      = apply_sqs,
                           ncores         = ncores)
  } else {
    x_norm <- x_counts
  }

  if (!is.data.frame(sfback)) {
    ## Defensive read of _sf.txt. The file is a 2-column whitespace-
    ## separated table written by GetBackground(). When absent, fall back to
    ## sf = 1 per feature (scaleMaxGexp is a no-op for sf == 1). When
    ## present but malformed, fail loudly so the caller knows the Ref dir
    ## is corrupted instead of silently scaling by 1.
    scaling.factors.file <- list.files(Ref_data, pattern = "_sf.txt$", full.names = TRUE)
    if (length(scaling.factors.file) > 1L) {
      stop(sprintf("PrepQueryData_ADT: multiple *_sf.txt found in %s: %s",
                   Ref_data, paste(basename(scaling.factors.file), collapse = ", ")))
    }
    if (length(scaling.factors.file) == 1L) {
      sf <- tryCatch(
        read.delim(scaling.factors.file[1], sep = " ", row.names = 1, header = FALSE),
        error = function(e) stop(sprintf(
          "PrepQueryData_ADT: failed to read scaling factors at %s (%s).",
          scaling.factors.file[1], conditionMessage(e))))
      if (!is.data.frame(sf) || nrow(sf) == 0L || ncol(sf) == 0L) {
        stop(sprintf("PrepQueryData_ADT: scaling factors file at %s is empty or malformed.",
                     scaling.factors.file[1]))
      }
    } else {
      sf <- data.frame(V2 = rep(1, nrow(x_norm)), row.names = rownames(x_norm))
    }
  } else {
    sf <- sfback
    if (nrow(sf) == 0L || ncol(sf) == 0L) {
      stop("PrepQueryData_ADT: `sfback` is empty.")
    }
  }

  if (scale_data) {
    feats_common <- intersect(rownames(x_norm), rownames(sf))
    x_scaled <- scaleMaxGexp(x_norm, sf, feats_common)
    if (saveData) saveRDS(x_scaled, file = file.path(dir_name, paste0(file.id, "_ScaledData.rds")))
  } else {
    x_scaled <- x_norm
  }
  if (return.data) return(x_scaled)
  invisible(NULL)
}

# Opt-in cell outlier filter (scuttle::isOutlier on log10 total + detected).
# Invoked from cleanData_ADT() when `qc_outlier_cells = "yes"`. Pre-unification
# this logic was always-on inside the CITE-seq cleanData_ADT body; it is now
# gated by an explicit flag and resolved by modality at the API entry point
# (run_scProteinXpress: "auto" => yes for CITE-seq, no for CyTOF).
.cell_outlier_filter <- function(data_counts, nmads = 3) {
  if (!requireNamespace("scuttle", quietly = TRUE)) {
    stop(".cell_outlier_filter: scuttle is required for qc_outlier_cells != 'no'.")
  }
  qc_input <- as.matrix(data_counts)
  qc_input[is.na(qc_input)] <- 0
  qc <- as.data.frame(scuttle::perCellQCMetrics(qc_input))
  qc$log10_sum <- log10(qc$sum + 1)
  low_sum      <- scuttle::isOutlier(qc$log10_sum, nmads = nmads, type = "lower")
  low_detected <- scuttle::isOutlier(qc$detected,  nmads = nmads, type = "lower")
  bad <- low_sum | low_detected
  list(Counts = data_counts[, !bad, drop = FALSE], QCMetrics = qc[!bad, , drop = FALSE])
}

#' Quality-control filter for ADT count matrices
#'
#' Applies a marker-level filter (both an absolute non-zero-cell threshold
#' and a coverage-fraction threshold), writes a per-marker exclusion audit
#' CSV when one is requested, and optionally applies a cell-level outlier
#' filter (scuttle::isOutlier on log10 total counts and on number of
#' detected markers) before the marker filter.
#'
#' @param data_counts proteins x cells numeric matrix.
#' @param nmads number of MADs for the scuttle cell-outlier filter.
#' @param min.cells.feature integer; minimum number of cells with a non-zero
#'   non-NA value for a marker to be retained.
#' @param min.coverage.fraction numeric in the range 0 to 1; minimum fraction of cells
#'   with any non-NA measurement for a marker to be retained.
#' @param audit_csv_path optional path. When set, exclusions are written here.
#' @param qc_outlier_cells one of "no" (default), "yes". When "yes",
#'   the internal cell-level outlier filter (scuttle::isOutlier on log10 total counts and detected markers) is applied before the marker-level filter. The
#'   value "auto" is resolved upstream in `run_scProteinXpress()`.
#' @return list(Counts, QCMetrics).
#' @export
cleanData_ADT <- function(data_counts, nmads = 3,
                          min.cells.feature = 10,
                          min.coverage.fraction = 0.10,
                          audit_csv_path = NULL,
                          qc_outlier_cells = c("no", "yes")){
  qc_outlier_cells <- match.arg(qc_outlier_cells)
  ## Argument sanity checks. Returning an empty matrix on bad input
  ## would surface as confusing downstream errors (e.g. `nrow == 0` deep
  ## inside ADTprot_norm). Fail at the QC boundary instead.
  if (!is.numeric(min.cells.feature) || min.cells.feature < 1L) {
    stop("cleanData_ADT: min.cells.feature must be >= 1; got ", min.cells.feature)
  }
  if (!is.numeric(min.coverage.fraction) ||
      min.coverage.fraction < 0 || min.coverage.fraction > 1) {
    stop("cleanData_ADT: min.coverage.fraction must be in the range 0 to 1; got ", min.coverage.fraction)
  }
  message("QC (CyTOF/CITE-seq): Filtering markers and empty cells (NAs treated as missing measurements).")
  data_counts <- .as_dense_matrix(data_counts)
  if (nrow(data_counts) == 0L || ncol(data_counts) == 0L) {
    stop(sprintf("cleanData_ADT: input is empty (%d markers x %d cells).",
                 nrow(data_counts), ncol(data_counts)))
  }
  if (is.null(rownames(data_counts))) {
    stop("cleanData_ADT: data_counts must have rownames (marker names).")
  }

  ## Opt-in scuttle cell-outlier filter (qc_outlier_cells = "yes"). Runs
  ## BEFORE marker-level filtering so that the marker coverage statistics
  ## reflect the cell population that will actually be used downstream.
  if (identical(qc_outlier_cells, "yes")) {
    n_before <- ncol(data_counts)
    outlier_res <- .cell_outlier_filter(data_counts, nmads = nmads)
    data_counts <- outlier_res$Counts
    message(sprintf(
      "[cleanData_ADT] scuttle isOutlier dropped %d / %d cells (nmads = %d).",
      n_before - ncol(data_counts), n_before, as.integer(nmads)))
  }

  # Cell-level QC: drop empty cells (no non-zero values across any feature).
  cell_sums   <- colSums(data_counts, na.rm = TRUE)
  valid_cells <- cell_sums > 0
  data_counts <- data_counts[, valid_cells, drop = FALSE]

  n_cells_total <- ncol(data_counts)
  if (n_cells_total == 0L) {
    return(list(Counts = data_counts,
                QCMetrics = data.frame(sum = numeric(0), detected = integer(0))))
  }

  # Feature-level QC: a marker passes iff it satisfies BOTH
  #   (a) absolute count: at least `min.cells.feature` cells with non-zero,
  #       non-NA value (mirrors scGeneXpress' >=10 threshold).
  #   (b) coverage fraction: at least `min.coverage.fraction` of cells have a
  #       non-NA measurement (this is a global filter that catches markers
  #       only present in a tiny fraction of the cohort due to multi-panel
  #       NA padding -- a problem absent in scRNA-seq).
  feat_nonzero_cells <- rowSums(!is.na(data_counts) & data_counts > 0)
  feat_observed_cells <- rowSums(!is.na(data_counts))
  feat_coverage <- feat_observed_cells / n_cells_total

  pass_abs  <- feat_nonzero_cells  >= min.cells.feature
  pass_frac <- feat_coverage       >= min.coverage.fraction
  keep_feat <- pass_abs & pass_frac

  excluded_idx <- which(!keep_feat)
  if (length(excluded_idx) > 0L) {
    reason <- ifelse(!pass_abs[excluded_idx] & !pass_frac[excluded_idx], "abs+frac",
              ifelse(!pass_abs[excluded_idx],                            "abs",
                                                                          "frac"))
    audit <- data.frame(
      marker            = rownames(data_counts)[excluded_idx],
      n_cells_observed  = feat_observed_cells[excluded_idx],
      n_cells_nonzero   = feat_nonzero_cells[excluded_idx],
      coverage_fraction = round(feat_coverage[excluded_idx], 4),
      reason            = reason,
      stringsAsFactors  = FALSE
    )
    message(sprintf(
      "[cleanData_ADT] Excluded %d/%d markers (min.cells.feature=%d, min.coverage.fraction=%.2f). Reasons: %s. Examples: %s",
      length(excluded_idx), nrow(data_counts),
      as.integer(min.cells.feature), min.coverage.fraction,
      paste(sprintf("%s=%d", names(table(reason)), as.integer(table(reason))), collapse = ", "),
      paste(head(audit$marker, 5), collapse = ", ")))

    if (!is.null(audit_csv_path)) {
      tryCatch({
        dir.create(dirname(audit_csv_path), showWarnings = FALSE, recursive = TRUE)
        utils::write.csv(audit, audit_csv_path, row.names = FALSE)
      }, error = function(e) {
        warning(sprintf("Could not write excluded-markers CSV to %s: %s",
                        audit_csv_path, conditionMessage(e)))
      })
    }
  } else {
    message(sprintf(
      "[cleanData_ADT] All %d markers passed filters (min.cells.feature=%d, min.coverage.fraction=%.2f).",
      nrow(data_counts), as.integer(min.cells.feature), min.coverage.fraction))
  }

  data_counts <- data_counts[keep_feat, , drop = FALSE]

  qc_mock <- data.frame(
    sum      = colSums(data_counts, na.rm = TRUE),
    detected = colSums(!is.na(data_counts) & data_counts > 0)
  )
  list(Counts = data_counts, QCMetrics = qc_mock)
}
.scaleMaxGexp_one <- function(xrow, sf_val) {
  if (is.na(sf_val) || sf_val == 0) return(xrow)
  xrow / sf_val
}

#' Scale a normalised matrix by per-marker factors
#'
#' Divides each `genes`-named row of `data_mat` by its scaling factor.
#'
#' @param data_mat proteins x cells matrix.
#' @param scalingfactors named numeric vector or single-column data.frame.
#' @param genes character vector of marker names to scale.
#' @return Scaled matrix.
#' @export
scaleMaxGexp <- function(data_mat, scalingfactors, genes){
  get_sf <- function(g) {
    if (is.numeric(scalingfactors) && !is.null(names(scalingfactors))) return(as.numeric(scalingfactors[g]))
    return(as.numeric(scalingfactors[g, 1]))
  }
  data_s <- foreach::foreach(n = seq_along(genes)) %do% {
    g <- genes[n]
    sfv <- get_sf(g)
    .scaleMaxGexp_one(data_mat[g, ], sfv)
  }
  names(data_s) <- genes
  do.call(rbind, data_s)
}

.optLowFun  <- function(x, empcdf, maxVal) abs(x * (1 - empcdf(x)))
.optHighFun <- function(x, empcdf, maxVal) abs(empcdf(x) * (maxVal - x))

minimizeRectangle <- function(empCDF, maxVal, minVal){
  lowThr  <- optimize(.optLowFun,  interval = c(minVal, maxVal), empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum
  highThr <- optimize(.optHighFun, interval = c(minVal, maxVal), empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum
  c(lowThr, highThr)
}

swap_th <- function(th_dist){
  lapply(th_dist, function(x) if (x[1] > x[2]) c(x[2], x[1]) else c(x[1], x[2]))
}

#' Rectangle-maximisation threshold (original scGeneXpress-derived heuristic)
#'
#' Finds the x that maximises x*(1-F(x)) (low cut) and F(x)*(max-x) (high
#' cut) on the empirical CDF of the pooled marker values. This is a single-
#' curve heuristic: it does not model the negative/positive (/intermediate)
#' populations explicitly, so for two populations that are close together
#' it can place both cuts on the same side of the true crossing point (see
#' `createThresholdDist_mixture()` benchmark in
#' `tests/testthat/test-threshold-mixture.R`), and for a genuinely bimodal
#' (2-population) marker it still always returns two distinct thresholds,
#' manufacturing a "medium" band that may not correspond to any real cell
#' state.
#' @noRd
createThresholdDist_rectangle <- function(gexp){
  gexp <- as.numeric(gexp)
  gexp <- gexp[is.finite(gexp)]
  if (length(gexp) == 0) return(c(0, 0))
  minv <- min(gexp)
  maxv <- max(gexp)
  if (minv >= maxv) return(c(minv, maxv))
  minimizeRectangle(ecdf(gexp), maxv, minv)
}

#' Mixture-model threshold: Bayes decision boundaries between fitted modes
#'
#' Fits a 1-3 component univariate Gaussian mixture (`mclust`, BIC-selected
#' G unless `G` is fixed) to the pooled marker values and returns the
#' crossing point(s) between adjacent weighted component densities -- the
#' actual Bayes-optimal decision boundary under the fitted mixture, which is
#' what "minimizes the theoretical rate of false-positive/false-negative
#' classifications" literally means (the rectangle heuristic approximates
#' this without fitting components at all).
#'
#' Falls back to `createThresholdDist_rectangle()` when: `mclust` is
#' unavailable, there are fewer than `min_n` finite values, the fit errors,
#' or BIC selects a single component (nothing to separate). For a
#' 2-component fit, both returned cuts equal the single crossing point
#' (no natural "medium" band exists for a strictly bimodal marker); for a
#' 3-component fit, the two cuts are the (1|2) and (2|3) crossings.
#'
#' @param gexp Numeric vector of pooled marker values for one bootstrap draw.
#' @param G Integer vector of candidate component counts passed to
#'   `mclust::Mclust` (default `1:3`, matching the low/medium/high design).
#' @param min_n Minimum number of finite values required to attempt a fit;
#'   below this, falls back to the rectangle heuristic (mixture fits are
#'   unstable on very small bootstrap draws).
#' @return `c(lowThr, highThr)`, unsorted (swapped by `swap_th()` downstream
#'   as with the rectangle method).
#' @noRd
createThresholdDist_mixture <- function(gexp, G = 1:3, min_n = 30){
  gexp <- as.numeric(gexp)
  gexp <- gexp[is.finite(gexp)]
  if (length(gexp) < min_n || length(unique(gexp)) < 3 ||
      !requireNamespace("mclust", quietly = TRUE)) {
    return(createThresholdDist_rectangle(gexp))
  }
  fit <- tryCatch(
    mclust::Mclust(gexp, G = G, verbose = FALSE, modelNames = "V"),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$G == 1) {
    return(createThresholdDist_rectangle(gexp))
  }

  pars <- fit$parameters
  mu   <- as.numeric(pars$mean)
  sig2 <- as.numeric(pars$variance$sigmasq)
  if (length(sig2) == 1) sig2 <- rep(sig2, fit$G)
  sig  <- sqrt(sig2)
  pro  <- as.numeric(pars$pro)
  ord  <- order(mu); mu <- mu[ord]; sig <- sig[ord]; pro <- pro[ord]

  .crossing <- function(i, j) {
    f <- function(x) pro[i] * stats::dnorm(x, mu[i], sig[i]) - pro[j] * stats::dnorm(x, mu[j], sig[j])
    tryCatch(stats::uniroot(f, lower = mu[i], upper = mu[j])$root,
             error = function(e) mean(c(mu[i], mu[j])))
  }

  if (fit$G == 2) {
    cut <- .crossing(1, 2)
    return(c(cut, cut))
  }
  c(.crossing(1, 2), .crossing(2, 3))
}

#' Quantile threshold baseline
#'
#' The simplest possible reference-informed baseline: call "low" below the
#' `low_q` quantile and "high" above the `high_q` quantile of the pooled/
#' bootstrapped reference values, with no bootstrap-distribution machinery
#' and no BH correction downstream. Exists so R1-Major4 can be answered
#' directly: rerun the SAME reference/query pipeline with
#' `threshold_method = "quantile"` and compare recall/precision against
#' `"rectangle"`/`"mixture"` on identical data, isolating whether
#' scProtXpress's advantage comes from the reference itself or from the
#' specific bootstrap/discretisation machinery.
#' @param gexp Numeric vector of pooled marker values for one bootstrap draw.
#' @param low_q,high_q Quantiles (default 10th/90th percentile).
#' @return `c(lowThr, highThr)`.
#' @noRd
createThresholdDist_quantile <- function(gexp, low_q = 0.10, high_q = 0.90){
  gexp <- as.numeric(gexp)
  gexp <- gexp[is.finite(gexp)]
  if (length(gexp) == 0) return(c(0, 0))
  if (length(gexp) == 1) return(c(gexp, gexp))
  stats::quantile(gexp, probs = c(low_q, high_q), names = FALSE, type = 7)
}

#' Z-score threshold baseline
#'
#' `c(mean - k*sd, mean + k*sd)` of the pooled/bootstrapped reference
#' values, clamped to the observed data range so it never proposes a
#' threshold outside what was actually measured. Second of the two simple
#' baselines requested for R1-Major4 (see `createThresholdDist_quantile()`).
#' @param gexp Numeric vector of pooled marker values for one bootstrap draw.
#' @param k Number of standard deviations (default 1).
#' @return `c(lowThr, highThr)`.
#' @noRd
createThresholdDist_zscore <- function(gexp, k = 1){
  gexp <- as.numeric(gexp)
  gexp <- gexp[is.finite(gexp)]
  if (length(gexp) == 0) return(c(0, 0))
  if (length(gexp) == 1) return(c(gexp, gexp))
  m <- mean(gexp); s <- stats::sd(gexp)
  if (!is.finite(s) || s == 0) return(c(m, m))
  lo <- max(m - k * s, min(gexp))
  hi <- min(m + k * s, max(gexp))
  c(lo, hi)
}

#' Dispatcher used by GetThresholds(); kept name-stable for back-compat.
#' `method = "rectangle"` reproduces the original behaviour exactly.
#' `method = "mixture"` uses `createThresholdDist_mixture()`.
#' `method = "quantile"`/`"zscore"` use the two simple baselines added for
#' R1-Major4 (see `createThresholdDist_quantile()`/`createThresholdDist_zscore()`).
#' @noRd
createThresholdDist <- function(gexp, method = c("rectangle", "mixture", "quantile", "zscore"), ...){
  method <- match.arg(method)
  switch(method,
    mixture  = createThresholdDist_mixture(gexp, ...),
    quantile = createThresholdDist_quantile(gexp, ...),
    zscore   = createThresholdDist_zscore(gexp, ...),
    createThresholdDist_rectangle(gexp)
  )
}

#' Per-marker threshold estimation over a bootstrap draw
#'
#' Dispatches each marker's pooled/resampled values to `createThresholdDist()`
#' and sorts the resulting (Lower, Upper) pair. Called once per bootstrap
#' iteration by `RunBootstrap()`.
#' @param threshold_method one of "rectangle" (default; original behaviour,
#'   bit-for-bit identical to pre-existing releases) or "mixture" (fits a
#'   1-3 component Gaussian mixture per bootstrap draw and uses the Bayes
#'   decision boundary; see `createThresholdDist_mixture()`). Threaded down
#'   from `run_scProteinXpress(threshold_method = ...)`.
#' @noRd
GetThresholds <- function(processed.back, file.id, threshold_method = c("rectangle", "mixture", "quantile", "zscore")){
  threshold_method <- match.arg(threshold_method)
  check <- names(which(sapply(processed.back, length) < 1))
  if (length(check) > 0) processed.back <- processed.back[setdiff(names(processed.back), check)]
  th_list <- lapply(processed.back, createThresholdDist, method = threshold_method)
  th_list <- swap_th(th_list)
  names(th_list) <- names(processed.back)
  th_list
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}

sample.data <- function(cells_dat, cellid, ncells.sample, nlevels){
  if(nlevels == 1){
    set <- split(cellid, as.factor(cellid))
    
    set_cells <- unlist(sapply(set, function(x){
      pick <- intersect(names(cells_dat), names(x))
      if(length(pick) >= 1 && length(x) >= 50){
        sub_cell <- sample(x = pick, size = ncells.sample, replace = TRUE)
      } else {
        sub_cell <- NA
      }
      return(sub_cell)
    }))
  } else {
    cell.id.split <- sapply(cellid, function(x) { strsplit(x, "[.]")[[1]][1] }) 
    set <- split(cellid, as.factor(cell.id.split))
    
    set_cells <- unlist(sapply(set, function(x){
      cell.id.split.tissue <- sapply(x, function(y) { strsplit(y, "[.]")[[1]][2] })
      set.tissue <- split(x, as.factor(cell.id.split.tissue))
      
      nbc <- unlist(lapply(set.tissue, length))
      n <- round(ncells.sample / length(set.tissue))
      set.tissue <- set.tissue[which(as.numeric(nbc) >= n)]
      
      if(length(set.tissue) == 0) return(NA)
      sn <- round(ncells.sample / length(set.tissue))
      
      set_cells_tis <- unlist(sapply(set.tissue, function(a){
        picktis <- intersect(names(cells_dat), names(a))
        if(length(picktis) >= 1 && length(a) >= sn){
          sub_cell <- sample(x = picktis, size = sn, replace = TRUE)
        } else {
          sub_cell <- NA
        }
        return(sub_cell)
      }))
      
      set_cells_tis <- set_cells_tis[!is.na(set_cells_tis)]
      if(length(set_cells_tis) >= 50){
        if(length(set_cells_tis) != ncells.sample){
          set_cells_tis <- sample(x = set_cells_tis, size = ncells.sample, replace = TRUE)
        }
      } else {
        set_cells_tis <- NA
      }
      return(set_cells_tis)
    }))
  }
  
  if(is.null(length(set_cells))){
    return(NULL)
  } else {
    set_cells <- set_cells[!is.na(set_cells)]
    if(length(set_cells) == 0) return(NULL)
    x <- cells_dat[set_cells]
    names(x) <- set_cells
    return(x)
  }
}

Background_sampling <- function(counts_norm_mat, cell_clusterid_mat, ncells.sample=100, clusterids_vec=c("cell.type"), min.cells.gexp=10){
  nlevels <- length(clusterids_vec)
  cellid <- apply(cell_clusterid_mat, 1, function(x){ paste(x[clusterids_vec], collapse = ".") })
  names(cellid) <- cell_clusterid_mat$cell.id
  
  markers <- rownames(counts_norm_mat)
  
  ids_samples_byGene <- lapply(markers, function(m){
    vals <- counts_norm_mat[m, ]
    vals <- vals[!is.na(vals) & vals > 0]
    if(length(vals) < min.cells.gexp) return(NULL)
    return(sample.data(cells_dat = vals, cellid = cellid, ncells.sample = ncells.sample, nlevels = nlevels))
  })
  
  names(ids_samples_byGene) <- markers
  ids_samples_byGene <- ids_samples_byGene[!sapply(ids_samples_byGene, is.null)]
  return(ids_samples_byGene)
}

RunBackComputations <- function(norm.data, cell_clusterid_mat, ncells.sample, clusterids_vec, file.id){
  sampled.back <- Background_sampling(
    counts_norm_mat = norm.data, 
    cell_clusterid_mat = cell_clusterid_mat, 
    ncells.sample = ncells.sample, 
    clusterids_vec = clusterids_vec, 
    min.cells.gexp = 10
  )
  GetThresholds(processed.back = sampled.back, file.id = file.id)
}
# -----------------------------------------------------------------------------
# .precompute_bootstrap_inputs
# -----------------------------------------------------------------------------

.precompute_bootstrap_inputs <- function(data_norm, metadata, clusterids_vec,
                                         min.cells.gexp = 10) {
  if (!all(clusterids_vec %in% colnames(metadata))) {
    stop(".precompute_bootstrap_inputs: faltan columnas de clusterids_vec en metadata")
  }
  cellid_vec <- do.call(paste, c(metadata[, clusterids_vec, drop = FALSE], sep = "."))
  names(cellid_vec) <- as.character(metadata$cell.id)
  
  nlevels <- length(clusterids_vec)
  cell_ids_all <- as.character(metadata$cell.id)
  
  if (nlevels == 1) {
    # split(cellid, as.factor(cellid)) in the original => groups by the
    # combination (which here is only cell.type). We replicate it.
    by_ct <- split(cell_ids_all, as.factor(cellid_vec))
    # Uniform structure: ct -> "_NONE_" -> cell.ids (to reuse the code
    # from the nlevels>=2 branch). The consumer will know it is nlevels==1 via a flag.
    groups_skeleton <- lapply(by_ct, function(ids) list("_NONE_" = ids))
  } else {
    # In the original, strsplit(cellid, "[.]")[[1]][1] and [2] are used to recover
    # cell.type and tissue. Since we already have the columns in metadata, we use them
    # directly (faster and free of the latent "dot-in-the-name" bug).
    ct_vec <- as.character(metadata[[clusterids_vec[1]]])
    ts_vec <- as.character(metadata[[clusterids_vec[2]]])
    
    # Alphabetical order of cell.types (== levels(as.factor(...)) == split() order)
    ct_levels <- levels(as.factor(ct_vec))
    groups_skeleton <- lapply(ct_levels, function(ct) {
      idx_ct <- which(ct_vec == ct)
      ts_sub <- ts_vec[idx_ct]
      ids_sub <- cell_ids_all[idx_ct]
      ts_levels <- levels(as.factor(ts_sub))
      out <- lapply(ts_levels, function(ts) ids_sub[ts_sub == ts])
      names(out) <- ts_levels
      out
    })
    names(groups_skeleton) <- ct_levels
  }
  markers <- rownames(data_norm)
  vals_by_marker <- vector("list", length(markers))
  groups_by_marker <- vector("list", length(markers))
  names(vals_by_marker) <- markers
  names(groups_by_marker) <- markers
  
  for (m in markers) {
    vals <- data_norm[m, ]
    keep <- !is.na(vals) & vals > 0
    vals <- vals[keep]
    if (length(vals) < min.cells.gexp) {
      # Equivalent to return(NULL) in the original lapply.
      vals_by_marker[[m]] <- NULL
      groups_by_marker[[m]] <- NULL
      next
    }
    vals_by_marker[[m]] <- vals  # named vector
    
    # Replicate EXACTLY the original intersect(names(cells_dat), names(a)):
    # the first argument dictates the order -> use names(vals).
    valid_ids <- names(vals)
    valid_set <- structure(seq_along(valid_ids), names = valid_ids)  # lookup O(1)
    
    groups_by_marker[[m]] <- lapply(groups_skeleton, function(ct_list) {
      lapply(ct_list, function(ids_in_group) {
        # Equivalent to intersect(names(vals), ids_in_group) preserving
        # the order of names(vals). match() returns NA for the missing ones;
        # we sort by position in valid_ids to replicate the order.
        pos <- valid_set[ids_in_group]
        pos <- pos[!is.na(pos)]
        if (length(pos) == 0) return(character(0))
        valid_ids[sort(unname(pos))]
      })
    })
  }
  
  list(
    nlevels          = nlevels,
    cellid_vec       = cellid_vec,
    markers          = markers,
    vals_by_marker   = vals_by_marker,
    groups_by_marker = groups_by_marker
  )
}


# -----------------------------------------------------------------------------
# .sample_data_fast
# Replicates the logic of sample.data() but operates on the structures
# pre-computed in .precompute_bootstrap_inputs(). full_lengths_for_marker
# stores the sizes of each group BEFORE the intersect with valid_ids, just
# like length(x) / length(set.tissue) in the original.
# -----------------------------------------------------------------------------
.sample_data_fast <- function(vals, groups_for_marker, full_lengths_for_marker,
                              ncells.sample, nlevels) {
  if (nlevels == 1) {
    # set <- split(cellid, as.factor(cellid))
    # In our nlevels==1 structure, groups_for_marker[[ct]][["_NONE_"]] are
    # the valid cell.ids; full_lengths_for_marker[[ct]][["_NONE_"]] is
    # length(x) (the full partition before the intersect).
    set_cells <- unlist(lapply(names(groups_for_marker), function(ct) {
      pick   <- groups_for_marker[[ct]][["_NONE_"]]
      x_len  <- full_lengths_for_marker[[ct]][["_NONE_"]]
      if (length(pick) >= 1 && x_len >= 50) {
        sample(x = pick, size = ncells.sample, replace = TRUE)
      } else {
        NA_character_
      }
    }))
  } else {
    # nlevels >= 2: replicates the logic of the original step by step.
    set_cells <- unlist(lapply(names(groups_for_marker), function(ct) {
      ts_list   <- groups_for_marker[[ct]]              # ts -> picktis
      ts_full_l <- full_lengths_for_marker[[ct]]        # ts -> full length
      ts_names  <- names(ts_list)                       # already in alphabetical order

      # n computed with length(set.tissue) BEFORE the filter.
      n_total_tissues <- length(ts_names)
      if (n_total_tissues == 0) return(NA_character_)
      n <- round(ncells.sample / n_total_tissues)

      # Filter tissues with length(a) >= n (length OF THE FULL GROUP, not of picktis)
      keep_ts <- ts_names[vapply(ts_names,
                                 function(ts) ts_full_l[[ts]] >= n,
                                 logical(1))]
      if (length(keep_ts) == 0) return(NA_character_)

      # sn recomputed with set.tissue post-filter (may differ from n).
      sn <- round(ncells.sample / length(keep_ts))

      # Loop over tissues in the SAME ORDER as the original split() (alphabetical).
      set_cells_tis <- unlist(lapply(keep_ts, function(ts) {
        picktis <- ts_list[[ts]]
        a_full_len <- ts_full_l[[ts]]
        if (length(picktis) >= 1 && a_full_len >= sn) {
          sample(x = picktis, size = sn, replace = TRUE)
        } else {
          NA_character_
        }
      }))
      
      set_cells_tis <- set_cells_tis[!is.na(set_cells_tis)]
      if (length(set_cells_tis) >= 50) {
        # Final Sample() adjustment to ncells.sample (identical to the original).
        if (length(set_cells_tis) != ncells.sample) {
          set_cells_tis <- sample(x = set_cells_tis, size = ncells.sample, replace = TRUE)
        }
        set_cells_tis
      } else {
        NA_character_
      }
    }))
  }
  
  if (is.null(length(set_cells))) return(NULL)
  set_cells <- set_cells[!is.na(set_cells)]
  if (length(set_cells) == 0) return(NULL)
  x <- vals[set_cells]
  names(x) <- set_cells
  x
}


# -----------------------------------------------------------------------------
# .Background_sampling_fast
# -----------------------------------------------------------------------------

.Background_sampling_fast <- function(precomp, ncells.sample) {
  markers <- precomp$markers
  nlevels <- precomp$nlevels
  
  ids_samples_byGene <- vector("list", length(markers))
  names(ids_samples_byGene) <- markers
  
  for (m in markers) {
    vals <- precomp$vals_by_marker[[m]]
    if (is.null(vals)) {
      ids_samples_byGene[[m]] <- NULL
      next
    }
    ids_samples_byGene[[m]] <- .sample_data_fast(
      vals                     = vals,
      groups_for_marker        = precomp$groups_by_marker[[m]],
      full_lengths_for_marker  = precomp$full_lengths_by_marker[[m]],
      ncells.sample            = ncells.sample,
      nlevels                  = nlevels
    )
  }
  
  # Same post-processing as the original:
  ids_samples_byGene <- ids_samples_byGene[!sapply(ids_samples_byGene, is.null)]
  ids_samples_byGene
}


# -----------------------------------------------------------------------------
# .RunBackComputations_fast
# -----------------------------------------------------------------------------

.RunBackComputations_fast <- function(precomp, ncells.sample, file.id,
                                      threshold_method = "rectangle") {
  sampled.back <- .Background_sampling_fast(precomp, ncells.sample)
  GetThresholds(processed.back = sampled.back, file.id = file.id,
               threshold_method = threshold_method)
}

# -----------------------------------------------------------------------------
# .SampleBackground_pooled
#
# Unstratified per-marker bootstrap sampler, mirroring the active path of
# the scGeneXpress reference implementation
# (`internalFunctions.R:338-345` of scGeneXpress-main): for each marker,
# draw `ncells.sample` values with replacement from the pool of finite
# values across all cells. `clusterids_vec` is intentionally absent here --
# this branch does not stratify.
#
# Used by RunBootstrap() when bootstrap_stratify == "none".
# -----------------------------------------------------------------------------
.SampleBackground_pooled <- function(norm.data, ncells.sample) {
  mat <- as.matrix(norm.data)
  markers <- rownames(mat)
  out <- lapply(seq_len(nrow(mat)), function(i) {
    vals <- mat[i, ]
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) return(numeric(0))
    sample(vals, ncells.sample, replace = TRUE)
  })
  names(out) <- markers
  out
}

.RunBackComputations_pooled <- function(norm.data, ncells.sample, file.id,
                                        threshold_method = "rectangle") {
  sampled.back <- .SampleBackground_pooled(norm.data, ncells.sample)
  GetThresholds(processed.back = sampled.back, file.id = file.id,
               threshold_method = threshold_method)
}
# -----------------------------------------------------------------------------
# RunBootstrap (optimized version)
# -----------------------------------------------------------------------------
#' Bootstrap per-marker threshold distributions
#'
#' Two-mode bootstrap. `bootstrap_stratify == "none"` (default) samples
#' each marker from the global pool of finite values across cells
#' with replacement, matching the scGeneXpress active path. The
#' `"cell.type"` and `"cell.type_tissue"` modes stratify by the matching
#' metadata columns using a one-shot precompute + per-iteration
#' `.RunBackComputations_fast()`. Either path is run inside
#' `foreach %dopar%` (fork on Unix, PSOCK on Windows) with BLAS/OMP/MKL
#' pinned to 1 thread per worker.
#'
#' @param data_norm proteins x cells normalised matrix.
#' @param dir_name Output directory.
#' @param file.id character(1) used as filename prefix.
#' @param clusterids_vec metadata columns to stratify by. Ignored when
#'   `bootstrap_stratify == "none"`.
#' @param metadata data.frame; must contain `cell.id` and `clusterids_vec`.
#' @param sample_n integer. Cells per bootstrap sample.
#' @param parallelize logical.
#' @param ncores integer. Honours SLURM_CPUS_PER_TASK.
#' @param num.boot.sample integer.
#' @param fixseed integer.
#' @param bootstrap_stratify one of "none" (default), "cell.type",
#'   "cell.type_tissue". Selects the sampling strategy.
#' @return list keyed by marker, each `n_boot x 2 (Lower, Upper)`.
#' @export
RunBootstrap <- function(data_norm, dir_name, file.id, clusterids_vec, metadata,
                         sample_n, parallelize, ncores, num.boot.sample, fixseed,
                         bootstrap_stratify = c("none", "cell.type", "cell.type_tissue"),
                         threshold_method = c("rectangle", "mixture", "quantile", "zscore")){
  bootstrap_stratify <- match.arg(bootstrap_stratify)
  threshold_method   <- match.arg(threshold_method)
  data_norm <- .as_dense_matrix(data_norm)
  set.seed(fixseed)
  cl <- NULL
  
  # --- Robust ncores resolution cross-platform / cross-scheduler ---------
  resolve_ncores <- function(ncores_req, n_tasks) {
    slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                                         unset = NA)))
    base_n <- if (is.numeric(ncores_req) && length(ncores_req) == 1L &&
                  !is.na(ncores_req) && ncores_req >= 1L) {
      as.integer(ncores_req)
    } else if (!is.na(slurm_cpus) && slurm_cpus >= 1L) {
      slurm_cpus
    } else {
      max(1L, parallel::detectCores() - 1L)
    }
    if (!is.na(slurm_cpus) && slurm_cpus >= 1L) {
      base_n <- min(base_n, slurm_cpus)
    }
    min(base_n, as.integer(n_tasks))
  }
  
  effective_ncores <- resolve_ncores(ncores, num.boot.sample)
  use_parallel <- isTRUE(parallelize) && effective_ncores >= 2L && num.boot.sample >= 2L
  
  # --- Pin BLAS/OMP/MKL to 1 thread (parent and workers) -------------------------
  .prev_blas_threads <- NA_integer_
  .prev_omp_threads  <- Sys.getenv("OMP_NUM_THREADS",      unset = NA)
  .prev_mkl_threads  <- Sys.getenv("MKL_NUM_THREADS",      unset = NA)
  .prev_openblas_t   <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  
  if (use_parallel) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      .prev_blas_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
                                     error = function(e) NA_integer_)
      try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(1L),  silent = TRUE)
    }
    Sys.setenv(OMP_NUM_THREADS      = "1")
    Sys.setenv(MKL_NUM_THREADS      = "1")
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  }
  
  # --- PRECOMPUTATION: only in stratified branches -------------------
  # bootstrap_stratify == "none" operates on `data_norm` directly (global
  # pool per marker), so there is nothing to precompute and `data_norm`
  # must remain alive when we enter the foreach.
  precomp <- NULL
  if (bootstrap_stratify != "none") {
  t_pre <- proc.time()
  precomp <- .precompute_bootstrap_inputs(
    data_norm      = data_norm,
    metadata       = metadata,
    clusterids_vec = clusterids_vec,
    min.cells.gexp = 10
  )

  build_full_lengths <- function(metadata, clusterids_vec, markers, nlevels) {
    if (nlevels == 1) {
      cellid_vec <- do.call(paste, c(metadata[, clusterids_vec, drop = FALSE], sep = "."))
      # length(x) in sample.data nlevels==1: the split by cellid groups all
      # cells with that cellid; length(x) is that count.
      tab <- table(cellid_vec)
      ct_levels <- names(tab)
      lens_one <- lapply(ct_levels, function(ct) list("_NONE_" = as.integer(tab[[ct]])))
      names(lens_one) <- ct_levels
      # It is independent of the marker, so we replicate the same list for
      # each marker to keep the interface uniform.
      out <- vector("list", length(markers))
      names(out) <- markers
      for (m in markers) out[[m]] <- lens_one
      return(out)
    }
    ct_vec <- as.character(metadata[[clusterids_vec[1]]])
    ts_vec <- as.character(metadata[[clusterids_vec[2]]])
    ct_levels <- levels(as.factor(ct_vec))
    lens_skeleton <- lapply(ct_levels, function(ct) {
      idx_ct <- which(ct_vec == ct)
      ts_sub <- ts_vec[idx_ct]
      ts_levels <- levels(as.factor(ts_sub))
      tab <- table(factor(ts_sub, levels = ts_levels))
      out <- as.list(as.integer(tab))
      names(out) <- ts_levels
      out
    })
    names(lens_skeleton) <- ct_levels
    # length(set.tissue) in the original is independent of the marker
    # (it is the partition of the full metadata by ct, ts).
    out <- vector("list", length(markers))
    names(out) <- markers
    for (m in markers) out[[m]] <- lens_skeleton
    out
  }
  precomp$full_lengths_by_marker <- build_full_lengths(
    metadata = metadata,
    clusterids_vec = clusterids_vec,
    markers = precomp$markers,
    nlevels = precomp$nlevels
  )
  message(sprintf("[TIMING][bootstrap] precompute elapsed: %.2fs",
                  (proc.time() - t_pre)[["elapsed"]]))
  rm(data_norm); invisible(gc(verbose = FALSE))
  }  # end if (bootstrap_stratify != "none")
  
  # --- Create cluster: fork on Unix, PSOCK on Windows -------------------------
  if (use_parallel) {
    try(unregister_dopar(), silent = TRUE)
    cl <- tryCatch({
      if (.Platform$OS.type == "unix") {
        parallel::makeForkCluster(nnodes = effective_ncores)
      } else {
        parallel::makePSOCKcluster(effective_ncores)
      }
    }, error = function(e) {
      message("[parallel][RunBootstrap] Could not create cluster (",
              conditionMessage(e), "); falling back to sequential.")
      NULL
    })
    
    if (!is.null(cl)) {
      doParallel::registerDoParallel(cl)
      if (.Platform$OS.type != "unix") {
        try(parallel::clusterEvalQ(cl, {
          Sys.setenv(OMP_NUM_THREADS      = "1",
                     MKL_NUM_THREADS      = "1",
                     OPENBLAS_NUM_THREADS = "1")
          if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
            try(RhpcBLASctl::omp_set_num_threads(1L),  silent = TRUE)
          }
        }), silent = TRUE)
      }
    }
  }
  
  n_markers_msg <- if (is.null(precomp)) nrow(data_norm) else length(precomp$markers)
  if (use_parallel && !is.null(cl)) {
    message(sprintf("[PARALLEL][bootstrap] Cluster created: %d workers (%s) for %d iterations (%d markers per iteration, stratify=%s).",
                    effective_ncores,
                    if (.Platform$OS.type == "unix") "fork" else "PSOCK",
                    num.boot.sample,
                    n_markers_msg,
                    bootstrap_stratify))
  } else {
    message(sprintf("[PARALLEL][bootstrap] SEQUENTIAL execution (ncores_requested=%s, effective=%d, parallelize=%s, stratify=%s).",
                    as.character(ncores), effective_ncores,
                    as.character(parallelize), bootstrap_stratify))
  }
  
  on.exit({
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
      try(foreach::registerDoSEQ(), silent = TRUE)
      try(unregister_dopar(), silent = TRUE)
    }
    if (!is.na(.prev_blas_threads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(.prev_blas_threads), silent = TRUE)
    }
    if (is.na(.prev_omp_threads)) Sys.unsetenv("OMP_NUM_THREADS")      else Sys.setenv(OMP_NUM_THREADS      = .prev_omp_threads)
    if (is.na(.prev_mkl_threads)) Sys.unsetenv("MKL_NUM_THREADS")      else Sys.setenv(MKL_NUM_THREADS      = .prev_mkl_threads)
    if (is.na(.prev_openblas_t))  Sys.unsetenv("OPENBLAS_NUM_THREADS") else Sys.setenv(OPENBLAS_NUM_THREADS = .prev_openblas_t)
  }, add = TRUE)
  
  # --- Main loop -------------------------------------------------------
  `%OP%` <- if (use_parallel && !is.null(cl)) foreach::`%dopar%` else foreach::`%do%`
  
  t_boot <- proc.time()
  th_list <- foreach::foreach(
    n = 1:num.boot.sample,
    .packages = c("Matrix", if (identical(threshold_method, "mixture")) "mclust"),
    .export   = c(".sample_data_fast", ".Background_sampling_fast",
                  ".RunBackComputations_fast",
                  ".SampleBackground_pooled", ".RunBackComputations_pooled",
                  ".optLowFun", ".optHighFun",
                  "minimizeRectangle", "createThresholdDist",
                  "createThresholdDist_rectangle", "createThresholdDist_mixture",
                  "createThresholdDist_quantile", "createThresholdDist_zscore",
                  "swap_th", "GetThresholds")
  ) %OP% {
    if (bootstrap_stratify == "none") {
      .RunBackComputations_pooled(data_norm, sample_n, file.id, threshold_method = threshold_method)
    } else {
      .RunBackComputations_fast(precomp, sample_n, file.id, threshold_method = threshold_method)
    }
  }
  message(sprintf("[TIMING][bootstrap] foreach elapsed: %.2fs",
                  (proc.time() - t_boot)[["elapsed"]]))
  
  # --- Post-processing identical to the original -------------------------------------
  ng <- length(th_list[[1]])
  thrs_list <- lapply(seq_len(ng), function(i){
    gthrs <- do.call("rbind", lapply(th_list, function(a) a[[i]]))
    colnames(gthrs) <- c("Lower", "Upper")
    gthrs
  })
  names(thrs_list) <- names(th_list[[1]])
  
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
    try(unregister_dopar(), silent = TRUE)
    cl <- NULL
  }
  
  thrs_list
}

#' Build the reference background threshold distributions
#'
#' Writes `<file.id>_sf.txt` (unit scaling factors) and `<file.id>_th_dist.rds`
#' (per-marker bootstrap thresholds) under `dir_name`.
#'
#' @param data_norm proteins x cells normalised matrix.
#' @param dir_name,file.id Output paths.
#' @param metadata data.frame.
#' @param sample_n,parallelize,ncores,num.boot.sample,fixseed Bootstrap controls.
#' @param bootstrap_stratify one of "none" (default), "cell.type",
#'   "cell.type_tissue". See `run_scProteinXpress()` for semantics.
#'   `"cell.type_tissue"` requires a `tissue` column in `metadata`;
#'   falls back to `"cell.type"` with a warning otherwise.
#' @param threshold_method one of "rectangle" (default, original behaviour)
#'   or "mixture" (Bayes decision boundary from a fitted 1-3 component
#'   Gaussian mixture per bootstrap draw; requires the `mclust` package).
#' @return list of per-marker threshold matrices.
#' @export
GetBackground <- function(data_norm, dir_name, file.id, metadata, sample_n, parallelize, ncores, num.boot.sample, fixseed,
                          bootstrap_stratify = c("none", "cell.type", "cell.type_tissue"),
                          threshold_method   = c("rectangle", "mixture", "quantile", "zscore")){
  bootstrap_stratify <- match.arg(bootstrap_stratify)
  threshold_method   <- match.arg(threshold_method)
  if (identical(threshold_method, "mixture") && !requireNamespace("mclust", quietly = TRUE)) {
    warning("GetBackground: threshold_method = 'mixture' requested but the ",
            "'mclust' package is not installed; falling back to 'rectangle' ",
            "for every marker (createThresholdDist_mixture() degrades ",
            "gracefully per-call, but installing mclust is recommended).",
            call. = FALSE)
  }
  data_norm <- .as_dense_matrix(data_norm)

  ## Explicit resolution of `clusterids_vec` from `bootstrap_stratify`,
  ## driven by the user's choice rather than an implicit heuristic based on
  ## presence/absence of `tissue`.
  if (identical(bootstrap_stratify, "cell.type_tissue") &&
      !all(c("cell.type", "tissue") %in% colnames(metadata))) {
    warning("GetBackground: requested cell.type_tissue stratification but ",
            "tissue column not present; falling back to cell.type.",
            call. = FALSE)
    bootstrap_stratify <- "cell.type"
  }
  clusterids_vec <- switch(bootstrap_stratify,
    "none"             = character(0),
    "cell.type"        = "cell.type",
    "cell.type_tissue" = c("cell.type", "tissue")
  )
  dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)

  feats <- rownames(data_norm)
  sf_vec <- rep(1, length(feats))
  names(sf_vec) <- feats
  utils::write.table(data.frame(feature = feats, sf = sf_vec), file = file.path(dir_name, paste0(file.id, "_sf.txt")), sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

  back.thrs <- RunBootstrap(data_norm, dir_name, file.id, clusterids_vec, metadata, sample_n, parallelize, ncores, num.boot.sample, fixseed,
                            bootstrap_stratify = bootstrap_stratify,
                            threshold_method   = threshold_method)
  saveRDS(back.thrs, file.path(dir_name, paste0(file.id, "_th_dist.rds")))

  ## Provenance sidecar: records the parameters that actually determined
  ## this reference, independent of what any driver script's comments say.
  ## Written every time GetBackground() runs so a stale/mismatched run is
  ## always self-documenting on disk.
  .write_provenance(dir_name, file.id,
                    bootstrap_stratify = bootstrap_stratify,
                    threshold_method   = threshold_method,
                    num.boot.sample    = num.boot.sample,
                    sample_n           = sample_n,
                    fixseed            = fixseed)
  back.thrs
}

#' Write a provenance sidecar (`<file.id>_provenance.json`, falling back to
#' `.txt` if `jsonlite` is unavailable) recording the exact parameters used
#' to build a reference, the package version, and (when available) the git
#' commit of the working directory the call was made from. This exists
#' because `bootstrap_stratify` (and now `threshold_method`) silently
#' resolve to a default when a driver script doesn't pass them explicitly
#' -- this file makes the *effective* choice visible on disk next to the
#' outputs it produced, rather than only in whatever driver script happened
#' to be used.
#' @noRd
.write_provenance <- function(dir_name, file.id, ...){
  params <- list(...)
  pkg_version <- tryCatch(as.character(utils::packageVersion("scProtXpress")),
                          error = function(e) NA_character_)
  git_commit <- tryCatch({
    out <- suppressWarnings(system("git rev-parse HEAD 2>/dev/null", intern = TRUE))
    if (length(out) == 1 && !grepl("fatal", out)) out else NA_character_
  }, error = function(e) NA_character_)

  record <- c(
    list(
      file.id             = file.id,
      timestamp_utc       = format(Sys.time(), tz = "UTC", usetz = TRUE),
      scProtXpress_version = pkg_version,
      git_commit          = git_commit,
      r_version           = R.version.string
    ),
    params
  )

  path_json <- file.path(dir_name, paste0(file.id, "_provenance.json"))
  ok <- FALSE
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    ok <- tryCatch({
      writeLines(jsonlite::toJSON(record, auto_unbox = TRUE, pretty = TRUE), path_json)
      TRUE
    }, error = function(e) FALSE)
  }
  if (!ok) {
    path_txt <- file.path(dir_name, paste0(file.id, "_provenance.txt"))
    lines <- vapply(names(record), function(nm) sprintf("%s: %s", nm, paste(record[[nm]], collapse = ",")), character(1))
    writeLines(lines, path_txt)
  }
  invisible(record)
}

compute_pval <- function(th_list, gexp_mat, genes, p_correct = TRUE){
  p_vals <- foreach::foreach(n = seq_along(genes)) %do% {
    g <- genes[n]
    p_exp  <- 1 - (ecdf(th_list[[g]][,"Upper"])(gexp_mat[g, ]))
    p_nexp <-      ecdf(th_list[[g]][,"Lower"])(gexp_mat[g, ])
    list(p_vals_exp = p_exp, p_vals_nexp = p_nexp)
  }
  
  if (p_correct) {
    sig_vals <- lapply(p_vals, function(x) list(sig_exp = stats::p.adjust(x$p_vals_exp, method = "BH"), sig_nexp = stats::p.adjust(x$p_vals_nexp, method = "BH")))
  } else {
    sig_vals <- lapply(p_vals, function(x) list(sig_exp = x$p_vals_exp, sig_nexp = x$p_vals_nexp))
  }
  
  names(sig_vals) <- genes
  sig_exp  <- do.call("rbind", lapply(sig_vals, `[[`, "sig_exp"))
  sig_nexp <- do.call("rbind", lapply(sig_vals, `[[`, "sig_nexp"))
  colnames(sig_exp)  <- colnames(gexp_mat)
  colnames(sig_nexp) <- colnames(gexp_mat)
  list(sig_exp = sig_exp, sig_nexp = sig_nexp)
}

classify <- function(sig_vals, th, th_i){
  sig_exp <- sig_vals[["sig_exp"]]
  sig_nexp <- sig_vals[["sig_nexp"]]
  discrete_exp <- matrix(0, nrow = nrow(sig_exp), ncol = ncol(sig_exp))
  rownames(discrete_exp) <- rownames(sig_exp)
  colnames(discrete_exp) <- colnames(sig_exp)
  cond_05 <- (sig_nexp > (1 - th_i)) & (sig_exp > th_i)
  discrete_exp[cond_05] <- 0.5
  cond_0 <- (sig_nexp <= th)
  discrete_exp[cond_0] <- 0
  cond_1 <- (sig_exp <= th)
  discrete_exp[cond_1] <- 1
  cond_na <- is.na(sig_exp) | is.na(sig_nexp)
  discrete_exp[cond_na] <- NA_real_
  return(discrete_exp)
}

#' Discretise a scaled query matrix against reference thresholds
#'
#' Computes per-cell p-values against the bootstrap-derived ECDFs and
#' classifies each cell to 0/0.5/1/NA at the configured `th` / `th_i`.
#'
#' @param query_ScaledData proteins x cells matrix.
#' @param Ref_data Path to the reference directory.
#' @param thsback Optional pre-loaded threshold list.
#' @param th Numeric. Significance threshold for 0/1.
#' @param th_i Numeric. Significance threshold for the 0.5 intermediate.
#' @param file.id character(1) used as filename prefix.
#' @param outdir Output directory.
#' @param return.data,saveData logical switches.
#' @param p_correct logical. Whether to BH-adjust the per-cell p-values
#'   computed against the reference ECDFs before classifying. Default
#'   `FALSE`, matching the pre-existing (previously hardcoded, now just
#'   made explicit) behaviour and the discretisation framework of
#'   scGeneXpress: per-cell BH correction over what is typically thousands
#'   to millions of cells is extremely conservative and not where the
#'   multiple-testing burden that matters actually lives. That burden is
#'   at the *population-level* calling step (one test per protein per cell
#'   type per condition) -- see `discrete.expr.levels()`'s new `fdr_method`
#'   argument and `correct_across_scan()` for where to actually control it.
#' @return list(discretised_data, significance).
#' @export
Discretise_data <- function(query_ScaledData, Ref_data, thsback = NA, th = 0.05, th_i = 0.05, file.id = "myquery", outdir = ".", return.data = TRUE, saveData = TRUE, p_correct = FALSE){
  dir_name <- file.path(outdir, paste0("Query_", file.id))
  if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE)
  threshold.dist <- if (length(thsback) == 1) {
    ## Explicit guard around the implicit `_th_dist.rds` resolution.
    f <- list.files(Ref_data, pattern = "_th_dist\\.rds$", full.names = TRUE)
    if (length(f) != 1L) {
      stop(sprintf("Discretise_data: expected exactly one *_th_dist.rds in %s; found %d.",
                   Ref_data, length(f)))
    }
    readRDS(f[1])
  } else thsback
  if (!is.list(threshold.dist) || length(threshold.dist) == 0L) {
    stop("Discretise_data: threshold distribution list is empty or malformed.")
  }
  genes_common <- intersect(rownames(query_ScaledData), names(threshold.dist))
  if (length(genes_common) == 0L) {
    stop(sprintf(
      "Discretise_data: no markers shared between query (n=%d) and threshold distributions (n=%d).",
      nrow(query_ScaledData), length(threshold.dist)))
  }
  query_pvals <- compute_pval(threshold.dist, query_ScaledData, genes_common, p_correct)
  saveRDS(query_pvals, file = file.path(dir_name, paste0(file.id, "_pvals.rds")))
  query_discretised <- classify(query_pvals, th, th_i)
  data <- list(discretised_data = query_discretised, significance = query_pvals)
  if (saveData) saveRDS(data, file = file.path(dir_name, paste0(file.id, "_DiscretisedData.rds")))
  if (return.data) return(data)
  invisible(NULL)
}
#' Compute prob_low/prob_high from an INDEPENDENT population, not the one
#' being tested
#'
#' `Create.sig.frame()` derives `prob_low`/`prob_high` -- the population-
#' level binomial baselines used by `discrete.expr.levels()` -- entirely
#' from the same discretised matrix that is about to be classified against
#' them: it pools the non-zero/medium rate across every protein measured in
#' THIS population, in THIS run. That means "is protein P enriched in
#' population X" is tested against "how leaky/noisy is population X on
#' average across ~150-700 other proteins", not against an independent
#' healthy/reference expectation. Populations with a high baseline rate of
#' non-specific antibody binding (platelets and erythroid cells are the
#' textbook CITE-seq case) get an inflated baseline that makes it *harder*
#' to call true enrichment there, while populations with low background
#' make almost anything look enriched by comparison.
#'
#' This function computes the SAME baseline statistics `Create.sig.frame()`
#' would, but from a caller-supplied independent discretised matrix (e.g.
#' a designated healthy/reference cohort's version of the same population,
#' or the Ref's own self-query) rather than from the query population
#' itself. Pass the result into `discrete.expr.levels(prob_low = ,
#' prob_high = )` (or `discrete.expr.levels_donor()`) for every disease
#' population INSTEAD of re-deriving it per-population from itself.
#'
#' @param reference_discrete.mat proteins x cells discretised matrix from
#'   an INDEPENDENT population/cohort (not the one being tested).
#' @param min.coverage.fraction.sigframe As in `Create.sig.frame()`.
#' @return list(prob_low, prob_high) -- pass straight through to
#'   `discrete.expr.levels()`/`discrete.expr.levels_donor()`.
#' @export
ComputeExternalBaseline <- function(reference_discrete.mat, min.coverage.fraction.sigframe = 0.10){
  sig <- Create.sig.frame(reference_discrete.mat,
                          min.coverage.fraction.sigframe = min.coverage.fraction.sigframe)
  list(prob_low = sig[[3]], prob_high = sig[[4]])
}

#' Isotype/negative-control-based background offset (optional refinement
#' on top of `ComputeExternalBaseline()`)
#'
#' If the panel includes isotype or other negative-control ADTs, their
#' non-zero rate in a given sample/population is a direct estimate of
#' non-specific binding -- a more targeted background estimate than
#' averaging across all real markers (which mixes in genuine biology).
#' When available, blend it into `prob_low` as a floor: the population is
#' never treated as having a lower expected background than its own
#' isotype controls actually show.
#'
#' @param discrete.mat proteins x cells discretised matrix for the
#'   population being tested (the SAME matrix passed to `Create.sig.frame`
#'   for that population, not an external one -- isotype controls
#'   estimate technical background per-sample, so they should come from the
#'   population itself, unlike `ComputeExternalBaseline()`'s biological
#'   baseline which should not).
#' @param isotype_marker_names Character vector of rownames in
#'   `discrete.mat` corresponding to isotype/negative-control ADTs.
#' @param prob_low The externally- or self-derived prob_low to floor.
#' @return Numeric: `max(prob_low, isotype non-zero rate)`.
#' @export
apply_isotype_floor <- function(discrete.mat, isotype_marker_names, prob_low){
  present <- intersect(isotype_marker_names, rownames(discrete.mat))
  if (length(present) == 0) {
    warning("apply_isotype_floor: none of the requested isotype_marker_names ",
            "are present in discrete.mat; returning prob_low unchanged.", call. = FALSE)
    return(prob_low)
  }
  iso <- discrete.mat[present, , drop = FALSE]
  n_nz  <- sum(iso > 0, na.rm = TRUE)
  n_tot <- sum(!is.na(iso))
  if (n_tot == 0) return(prob_low)
  iso_rate <- n_nz / n_tot
  max(prob_low, iso_rate)
}

#' Per-gene signature frame with coverage filter
#'
#' Returns a 5-element list:
#' \enumerate{
#'   \item sig.frame data.frame (gene, n.0, n.0.5, n.1, not.0, not.0.perc,
#'         pct_one, coverage_fraction, n_observed).
#'   \item z2 (legacy z-score retained for back-compat).
#'   \item prob_low (per-pop binomial baseline for the non-zero rate).
#'   \item prob_high (per-pop binomial baseline for 0.5 among non-zero).
#'   \item excluded_genes character vector of genes excluded from the
#'         prob_low / prob_high computation by the coverage filter.
#' }
#'
#' @param discrete.mat proteins x cells discretised matrix.
#' @param min.coverage.fraction.sigframe Numeric in the range 0 to 1.
#' @param audit_csv_path Optional path; writes the excluded-from-baseline
#'   gene table when set.
#' @return list of length 5 (see Description).
#' @noRd
Create.sig.frame <- function(discrete.mat,
                             min.coverage.fraction.sigframe = 0.10,
                             audit_csv_path = NULL){
  ## Argument sanity. `discrete.mat` must be a 0/0.5/1/NA matrix with
  ## rownames (markers). Anything else is a programming error upstream.
  if (!is.matrix(discrete.mat) && !inherits(discrete.mat, "Matrix")) {
    stop("Create.sig.frame: discrete.mat must be a matrix; got ",
         paste(class(discrete.mat), collapse = ","))
  }
  if (is.null(rownames(discrete.mat))) {
    stop("Create.sig.frame: discrete.mat must have rownames (marker names).")
  }
  if (!is.numeric(min.coverage.fraction.sigframe) ||
      min.coverage.fraction.sigframe < 0 ||
      min.coverage.fraction.sigframe > 1) {
    stop("Create.sig.frame: min.coverage.fraction.sigframe must be in the range 0 to 1.")
  }
  # Methodology:
  #   - prob_low / prob_high are computed ONLY over genes whose per-population
  #     coverage (n_observed / n_cells_in_pop) is >= min.coverage.fraction.sigframe.
  #     This protects the population-level baseline from being contaminated by
  #     markers that are technically present in the matrix (they passed the
  #     global QC step) but happen to have very few observations in this
  #     specific population due to multi-panel NA padding.
  #   - Every gene STILL receives its own call in discrete.expr.levels() against
  #     these sanitised probs, including those excluded from the prob calc.
  #     Excluded-from-baseline does NOT mean excluded-from-output.
  #
  # Schema of sig.frame (columns):
  #   gene, n.0, n.0.5, n.1, not.0, not.0.perc, pct_one,
  #   coverage_fraction, n_observed
  # (expr.level is added later by discrete.expr.levels).

  discrete.mat[is.nan(discrete.mat)] <- NA
  n_cells_pop <- ncol(discrete.mat)
  if (n_cells_pop == 0L) {
    return(list(data.frame(), 0, 1e-6, 1e-6, character(0)))
  }

  # Per-gene vectorised counts (single matrix scan per category)
  m0   <- discrete.mat == 0
  m05  <- discrete.mat == 0.5
  m1   <- discrete.mat == 1
  n0   <- rowSums(m0,  na.rm = TRUE)
  n05  <- rowSums(m05, na.rm = TRUE)
  n1   <- rowSums(m1,  na.rm = TRUE)
  n_observed <- n0 + n05 + n1

  keep_gene <- n_observed > 0L
  if (!any(keep_gene)) return(list(data.frame(), 0, 1e-6, 1e-6, character(0)))

  not0  <- (n05 + n1)[keep_gene]
  not0p <- not0 / n_observed[keep_gene]
  pct_one <- n1[keep_gene] / pmax(n_observed[keep_gene], 1L)
  coverage_fraction <- n_observed[keep_gene] / n_cells_pop

  sig.frame <- data.frame(
    gene              = rownames(discrete.mat)[keep_gene],
    n.0               = n0[keep_gene],
    n.0.5             = n05[keep_gene],
    n.1               = n1[keep_gene],
    not.0             = not0,
    not.0.perc        = not0p,
    pct_one           = pct_one,
    coverage_fraction = coverage_fraction,
    n_observed        = n_observed[keep_gene],
    stringsAsFactors  = FALSE
  )
  rownames(sig.frame) <- sig.frame$gene
  sig.frame <- sig.frame[order(sig.frame$not.0.perc, decreasing = TRUE), , drop = FALSE]

  # Legacy z-score statistic (kept for back-compat with previous outputs/plots).
  p_global <- sum(discrete.mat > 0, na.rm = TRUE) /
              (nrow(discrete.mat) * ncol(discrete.mat))
  z2 <- (sqrt(n_cells_pop * p_global * (1 - p_global)) * 2) + (n_cells_pop * p_global)

  # Per-population coverage filter: keep only genes with sufficient coverage
  # to contribute to the binomial baseline.
  # NOTE: sig.frame has been reordered by not.0.perc above, so we must derive
  # passes_cov from the reordered sig.frame$coverage_fraction (NOT from the
  # pre-reorder vector, which would be misaligned with sig.frame$gene).
  passes_cov <- sig.frame$coverage_fraction >= min.coverage.fraction.sigframe
  excluded_genes <- sig.frame$gene[!passes_cov]
  used_genes     <- sig.frame$gene[ passes_cov]

  if (length(used_genes) == 0L) {
    warning(sprintf(
      "[Create.sig.frame] No genes meet min.coverage.fraction.sigframe=%.2f in this population (n_cells_pop=%d). %s",
      min.coverage.fraction.sigframe, n_cells_pop,
      "Falling back to all genes for the binomial baseline."))
    used_idx <- rep(TRUE, nrow(sig.frame))
  } else {
    used_idx <- passes_cov
  }

  # Recompute prob_low / prob_high on the subset of trustworthy genes only.
  # We sum the per-gene counts of those genes (this equals scanning the
  # discrete matrix restricted to those gene rows and all their non-NA cells).
  N0_sub  <- sum(sig.frame$n.0[used_idx])
  N05_sub <- sum(sig.frame$n.0.5[used_idx])
  N1_sub  <- sum(sig.frame$n.1[used_idx])
  Nnz_sub <- N05_sub + N1_sub
  Nall_sub <- N0_sub + Nnz_sub

  prob_low  <- if (Nall_sub > 0) Nnz_sub / Nall_sub else 1e-6
  prob_high <- if (Nnz_sub  > 0) N05_sub / Nnz_sub  else 1e-6
  if (prob_low  >= 1) prob_low  <- 1 - 1e-6
  if (prob_high >= 1) prob_high <- 1 - 1e-6
  if (prob_low  <= 0) prob_low  <- 1e-6
  if (prob_high <= 0) prob_high <- 1e-6

  # Audit CSV: which genes were excluded from the prob_low/high computation
  # in this population, with their coverage and counts.
  if (length(excluded_genes) > 0L && !is.null(audit_csv_path)) {
    audit <- sig.frame[!passes_cov,
                       c("gene", "n_observed", "coverage_fraction",
                         "n.0", "n.0.5", "n.1"),
                       drop = FALSE]
    audit$min.coverage.fraction.sigframe <- min.coverage.fraction.sigframe
    audit$n_cells_in_population <- n_cells_pop
    rownames(audit) <- NULL
    tryCatch({
      dir.create(dirname(audit_csv_path), showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(audit, audit_csv_path, row.names = FALSE)
      message(sprintf(
        "[Create.sig.frame] %d/%d genes excluded from prob_low/high baseline (coverage<%.2f). Audit: %s",
        length(excluded_genes), nrow(sig.frame),
        min.coverage.fraction.sigframe, audit_csv_path))
    }, error = function(e) {
      warning(sprintf("Could not write sig_frame audit CSV: %s", conditionMessage(e)))
    })
  } else if (length(excluded_genes) > 0L) {
    message(sprintf(
      "[Create.sig.frame] %d/%d genes excluded from prob_low/high baseline (coverage<%.2f).",
      length(excluded_genes), nrow(sig.frame), min.coverage.fraction.sigframe))
  }

  list(sig.frame, z2, prob_low, prob_high, excluded_genes)
}
#' Per-donor gene x donor count breakdown, for donor-aware population testing
#'
#' `Create.sig.frame()` pools every cell in a population regardless of which
#' donor/sample it came from, so the downstream `pbinom` test in
#' `discrete.expr.levels()` treats each CELL as an independent trial. Cells
#' from the same donor are correlated (shared genotype, disease state,
#' staining batch), so that test's false-positive rate inflates whenever
#' there is real between-donor heterogeneity -- confirmed by simulation to
#' reach 3-8x the nominal rate at ICC values (0.02-0.20) well within what is
#' typically seen for immune cell-state proportions across donors, even
#' though it is exactly nominal when donors are truly identical
#' (`tests/testthat/test-calibration.R`).
#'
#' This function does NOT replace `Create.sig.frame()` -- it computes the
#' same n.0/n.0.5/n.1 counts *broken out per donor* so that
#' `discrete.expr.levels_donor()` can test across donors instead of across
#' pooled cells. Everything else about the pipeline (thresholds,
#' discretisation) is unchanged.
#'
#' @param discrete.mat proteins x cells discretised matrix (0/0.5/1/NA).
#' @param donor_ids Character/factor vector of length `ncol(discrete.mat)`,
#'   giving the donor/sample each cell belongs to. Must align positionally
#'   with the columns of `discrete.mat`.
#' @return list(n0, n05, n1): each a genes x donors integer matrix.
#' @export
Create.sig.frame_by_donor <- function(discrete.mat, donor_ids){
  if (ncol(discrete.mat) != length(donor_ids)) {
    stop("Create.sig.frame_by_donor: length(donor_ids) must equal ncol(discrete.mat).")
  }
  donor_ids <- as.character(donor_ids)
  donors <- sort(unique(donor_ids))
  genes  <- rownames(discrete.mat)

  n0  <- matrix(0L, nrow = length(genes), ncol = length(donors), dimnames = list(genes, donors))
  n05 <- n0; n1 <- n0

  for (d in donors) {
    sub <- discrete.mat[, donor_ids == d, drop = FALSE]
    n0[,  d] <- rowSums(sub == 0,   na.rm = TRUE)
    n05[, d] <- rowSums(sub == 0.5, na.rm = TRUE)
    n1[,  d] <- rowSums(sub == 1,   na.rm = TRUE)
  }
  list(n0 = n0, n05 = n05, n1 = n1)
}

#' One-way ANOVA (Fleiss-Cuzick-style) estimator of the intra-donor
#' correlation (ICC) for a binary/proportion outcome measured on cells
#' clustered within donors, used by `discrete.expr.levels_donor(test_method
#' = "design_effect")` to deflate the effective sample size before applying
#' the package's existing `pbinom` test. Returns 0 (no correction, falls
#' back to the original pooled behaviour) when there are fewer than 2
#' donors or the between/within decomposition is degenerate.
#' @noRd
.estimate_icc <- function(nnz_per_donor, ntot_per_donor){
  D <- length(ntot_per_donor)
  if (D < 2) return(0)
  n_tot <- sum(ntot_per_donor)
  if (n_tot == 0) return(0)
  p_d   <- nnz_per_donor / pmax(ntot_per_donor, 1)
  p_bar <- sum(nnz_per_donor) / n_tot

  MSB <- sum(ntot_per_donor * (p_d - p_bar)^2) / (D - 1)
  MSW <- sum(ntot_per_donor * p_d * (1 - p_d)) / max(n_tot - D, 1)
  n0  <- (n_tot - sum(ntot_per_donor^2) / n_tot) / (D - 1)

  denom <- MSB + (n0 - 1) * MSW
  icc <- if (is.finite(denom) && denom > 0) (MSB - MSW) / denom else 0
  max(0, min(icc, 0.99))
}

#' Donor-aware alternative to the pooled-cell population test
#'
#' Same output contract and low/not.significant/medium/high cascade as
#' `discrete.expr.levels()`, but the significance test is computed across
#' donors rather than across pooled cells.
#'
#' @param by_donor Output of `Create.sig.frame_by_donor()`.
#' @param prob_low,prob_high Population-level binomial baselines, as from
#'   `Create.sig.frame()` (unchanged -- only the *test*, not the baseline
#'   estimate, changes here; see `ComputeExternalBaseline()` if you also
#'   want the baseline itself to stop being self-referential).
#' @param test_method one of:
#'   \itemize{
#'     \item `"donor_t"` (default) -- one-sample t-test on arcsine-root
#'       transformed per-donor proportions against `asin(sqrt(prob_low))`
#'       (and analogously for the medium/high sub-test). Simplest to
#'       explain ("we test across donors, not cells"); slightly
#'       conservative at high ICC / many donors in simulation.
#'     \item `"donor_wilcox"` -- Wilcoxon signed-rank version of the above;
#'       avoids the normality assumption, useful with few donors or skewed
#'       proportions.
#'     \item `"design_effect"` -- keeps the exact pooled-cell `pbinom` test
#'       from `discrete.expr.levels()`, but deflates `n`/`n_nz` by an
#'       ANOVA-estimated design effect (`1 + (mean donor size - 1) * ICC`)
#'       before applying it. Tracked closest to nominal across the widest
#'       range of simulated ICC values in
#'       `tests/testthat/test-calibration.R` and requires the smallest
#'       conceptual change to the existing statistic.
#'   }
#'   Populations with fewer than 2 donors contributing non-NA cells for a
#'   gene fall back to the original pooled `pbinom` test for that gene
#'   (with a warning tally in the attributes of the returned table).
#' @param min_donors_for_call Integer, default 2. Genes with fewer
#'   contributing donors than this are forced to `"not.significant"`
#'   rather than silently falling back to a pooled-cell test that offers no
#'   real protection against pseudoreplication.
#' @return data.frame: gene, per-donor-summed n.0/n.0.5/n.1 (for
#'   comparability with `discrete.expr.levels()`'s output shape), n_donors,
#'   p_low, p_high, expr.level.
#' @export
discrete.expr.levels_donor <- function(by_donor, prob_low, prob_high,
                                       test_method = c("donor_t", "donor_wilcox", "design_effect"),
                                       min_donors_for_call = 2L){
  test_method <- match.arg(test_method)
  if (prob_low  >= 1) prob_low  <- 1 - 1e-6
  if (prob_high >= 1) prob_high <- 1 - 1e-6
  if (prob_low  <= 0) prob_low  <- 1e-6
  if (prob_high <= 0) prob_high <- 1e-6

  genes <- rownames(by_donor$n0)
  n_fallback <- 0L

  .test_low <- function(nnz_d, ntot_d) {
    switch(test_method,
      donor_t = {
        p_d <- nnz_d / ntot_d
        y <- asin(sqrt(p_d)); mu0 <- asin(sqrt(prob_low))
        if (length(p_d) < 2 || sd(y) == 0) return(NA_real_)
        tryCatch(stats::t.test(y, mu = mu0, alternative = "greater")$p.value, error = function(e) NA_real_)
      },
      donor_wilcox = {
        p_d <- nnz_d / ntot_d
        if (length(p_d) < 2 || all(p_d == p_d[1])) return(NA_real_)
        tryCatch(suppressWarnings(stats::wilcox.test(p_d, mu = prob_low, alternative = "greater")$p.value),
                 error = function(e) NA_real_)
      },
      design_effect = {
        icc <- .estimate_icc(nnz_d, ntot_d)
        n_tot <- sum(ntot_d); m_bar <- mean(ntot_d)
        deff  <- 1 + (m_bar - 1) * icc
        n_eff   <- max(round(n_tot / deff), length(ntot_d))
        nnz_eff <- round(sum(nnz_d) / n_tot * n_eff)
        pbinom(max(nnz_eff - 1L, 0L), size = n_eff, prob = prob_low, lower.tail = FALSE)
      }
    )
  }
  ## Symmetric helper for the medium/high sub-test (tests the 0.5-share
  ## among non-zero cells against prob_high, same three methods).
  .test_high <- function(n05_d, nnz_d) {
    if (sum(nnz_d) == 0) return(NA_real_)
    keep <- nnz_d > 0
    if (sum(keep) < 2) return(NA_real_)
    switch(test_method,
      donor_t = {
        p_d <- n05_d[keep] / nnz_d[keep]
        y <- asin(sqrt(p_d)); mu0 <- asin(sqrt(prob_high))
        if (sd(y) == 0) return(NA_real_)
        tryCatch(stats::t.test(y, mu = mu0, alternative = "greater")$p.value, error = function(e) NA_real_)
      },
      donor_wilcox = {
        p_d <- n05_d[keep] / nnz_d[keep]
        if (all(p_d == p_d[1])) return(NA_real_)
        tryCatch(suppressWarnings(stats::wilcox.test(p_d, mu = prob_high, alternative = "greater")$p.value),
                 error = function(e) NA_real_)
      },
      design_effect = {
        icc <- .estimate_icc(n05_d[keep], nnz_d[keep])
        n_tot <- sum(nnz_d[keep]); m_bar <- mean(nnz_d[keep])
        deff  <- 1 + (m_bar - 1) * icc
        n_eff   <- max(round(n_tot / deff), sum(keep))
        n05_eff <- round(sum(n05_d[keep]) / n_tot * n_eff)
        pbinom(max(n05_eff - 1L, 0L), size = n_eff, prob = prob_high, lower.tail = FALSE)
      }
    )
  }

  out <- lapply(genes, function(g) {
    n0_d  <- by_donor$n0[g, ]; n05_d <- by_donor$n05[g, ]; n1_d <- by_donor$n1[g, ]
    ntot_d <- n0_d + n05_d + n1_d
    nnz_d  <- n05_d + n1_d
    keep_d <- ntot_d > 0   # donors that actually contributed measured cells to this gene
    n_donors <- sum(keep_d)

    N0 <- sum(n0_d); N05 <- sum(n05_d); N1 <- sum(n1_d)
    row <- data.frame(gene = g, n.0 = N0, n.0.5 = N05, n.1 = N1, n_donors = n_donors,
                      stringsAsFactors = FALSE)

    if (n_donors < min_donors_for_call) {
      row$p_low <- NA_real_; row$p_high <- NA_real_; row$expr.level <- "not.significant"
      return(row)
    }

    p_low  <- .test_low(nnz_d[keep_d], ntot_d[keep_d])
    if (is.na(p_low)) { n_fallback <<- n_fallback + 1L; p_low <- 1 }
    row$p_low <- p_low

    if (p_low > 0.05 || (N05 + N1) == 0) {
      row$p_high <- NA_real_; row$expr.level <- "low"; return(row)
    }
    if (N05 == N1) { row$p_high <- NA_real_; row$expr.level <- "not.significant"; return(row) }

    p_high <- .test_high(n05_d[keep_d], nnz_d[keep_d])
    if (is.na(p_high)) { n_fallback <<- n_fallback + 1L; p_high <- 1 }
    row$p_high <- p_high
    row$expr.level <- if (p_high < 0.05) "medium" else "high"
    row
  })
  result <- do.call(rbind, out)
  attr(result, "n_fallback_to_pooled") <- n_fallback
  if (n_fallback > 0) {
    message(sprintf("[discrete.expr.levels_donor] %d/%d genes had a degenerate donor test (e.g. all donors identical) and were scored NA/conservative; check attr(result, 'n_fallback_to_pooled').",
                    n_fallback, length(genes)))
  }
  result
}

#' Assign per-gene expression levels with two binomial tests
#'
#' Vectorised one-sided `pbinom` (alternative = "greater") on each gene's
#' non-zero count versus the per-population baselines `prob_low` and
#' `prob_high` from `Create.sig.frame()`. Returns "low", "medium", "high"
#' or "not.significant" in a new `expr.level` column.
#'
#' Back-compat: `discrete.expr.levels(tbl, zs)` (2 arguments) still works
#' but emits `.Deprecated`. `zs` is no longer used and is kept only for
#' positional API stability; pass `prob_low` and `prob_high` explicitly
#' (preferably from `Create.sig.frame()`) to silence the warning.
#'
#' @param tbl data.frame returned by `Create.sig.frame()` (element 1).
#' @param zs Legacy z-score, ignored.
#' @param prob_low,prob_high Numerics; per-population binomial baselines.
#' @param fdr_method NULL (default; original behaviour, `expr.level` is
#'   thresholded on the raw one-sided p-values exactly as before) or "BH".
#'   When "BH", `p_low`/`p_high` are Benjamini-Hochberg adjusted *across the
#'   genes in `tbl`* before both the reported `q_low`/`q_high` columns and
#'   the `expr.level` cutoff -- this controls the FDR within one
#'   population's whole-panel scan. It does NOT correct across multiple
#'   populations/diseases in a larger discovery run; use
#'   `correct_across_scan()` for that, on the list of tables this function
#'   returns.
#' @param min_n_for_call Integer, default 0 (original behaviour: any n
#'   yields a confident call, `"not.significant"` is only assigned on an
#'   exact `n.0.5 == n.1` tie). When > 0, any gene with `ntot` (or `nnz`,
#'   for the high/medium sub-test) below this value is additionally forced
#'   to `"not.significant"` regardless of how the binomial test comes out,
#'   since a handful of cells should not generate a confident population
#'   call. Recommended: 20-50 depending on how much you trust small
#'   populations in your data.
#' @param effect_size logical, default `TRUE`. Adds `effect_size_low`
#'   (`nnz/ntot - prob_low`) and `effect_size_high` (`n.0.5/nnz -
#'   prob_high`, `NA` when `nnz == 0`) columns -- purely additive, does not
#'   change `expr.level`. Lets you rank/filter calls by strength rather
#'   than only by the categorical label: two "high" calls are not equally
#'   strong evidence.
#' @return The input data.frame with `expr.level` plus (new, additive)
#'   `p_low`, `p_high` (previously computed then discarded), and, per the
#'   arguments above, `q_low`/`q_high` and/or `effect_size_low`/
#'   `effect_size_high`.
#' @noRd
discrete.expr.levels <- function(tbl, zs, prob_low = NULL, prob_high = NULL,
                                 fdr_method = NULL, min_n_for_call = 0L,
                                 effect_size = TRUE){
  # Per-gene call using two binomial tests:
  #   bt_low : P(non-zero >= n_nz | n_total, prob_low),  alternative="greater"
  #            -> rejects "low" if the gene has more non-zero cells than expected.
  #   bt_high: P(0.5 >= n_05 | n_nz, prob_high),         alternative="greater"
  #            -> rejects "high" if among non-zero cells most are intermediate.
  # Uses the population-level prob_low/prob_high supplied by Create.sig.frame
  # (computed only on genes with sufficient per-population coverage). When not
  # supplied, recovers them from the table itself for backward compatibility.

  if (nrow(tbl) == 0) { tbl$expr.level <- character(0); return(tbl) }

  if (!is.null(fdr_method) && !identical(fdr_method, "BH")) {
    stop('discrete.expr.levels: fdr_method must be NULL or "BH".')
  }
  if (!is.numeric(min_n_for_call) || min_n_for_call < 0) {
    stop("discrete.expr.levels: min_n_for_call must be a non-negative integer.")
  }

  if (is.null(prob_low) || is.null(prob_high)){
    N0  <- sum(tbl$n.0)
    N05 <- sum(tbl$n.0.5)
    N1  <- sum(tbl$n.1)
    Nnz <- N05 + N1
    Nall <- N0 + Nnz
    prob_low  <- if (Nall > 0) Nnz / Nall else 1e-6
    prob_high <- if (Nnz  > 0) N05 / Nnz  else 1e-6
  }
  if (prob_low  >= 1) prob_low  <- 1 - 1e-6
  if (prob_high >= 1) prob_high <- 1 - 1e-6
  if (prob_low  <= 0) prob_low  <- 1e-6
  if (prob_high <= 0) prob_high <- 1e-6

  n05  <- as.numeric(tbl$n.0.5)
  n1   <- as.numeric(tbl$n.1)
  nnz  <- n05 + n1
  ntot <- as.numeric(tbl$n.0) + nnz

  # Vectorised one-sided binomial via pbinom:
  # P(X >= k | n, p) = pbinom(k-1, n, p, lower.tail = FALSE).
  p_low  <- pbinom(pmax(nnz - 1L, 0L), size = ntot, prob = prob_low,
                   lower.tail = FALSE)
  p_high <- ifelse(nnz > 0,
                   pbinom(pmax(n05 - 1L, 0L), size = pmax(nnz, 1L),
                          prob = prob_high, lower.tail = FALSE),
                   NA_real_)

  ## --- new, additive: raw p-values were computed above then silently
  ## discarded in every prior release; keep them, and BH-adjust across the
  ## gene panel in `tbl` when requested. `cutoff_low`/`cutoff_high` are what
  ## `is_low`/`med` are actually thresholded against below -- 0.05 on the
  ## raw p when fdr_method is NULL (byte-identical to the original), or the
  ## BH q-value when fdr_method = "BH".
  if (identical(fdr_method, "BH")) {
    q_low  <- stats::p.adjust(p_low,  method = "BH")
    q_high <- stats::p.adjust(p_high, method = "BH")
    cutoff_low  <- q_low
    cutoff_high <- q_high
  } else {
    q_low  <- rep(NA_real_, length(p_low))
    q_high <- rep(NA_real_, length(p_high))
    cutoff_low  <- p_low
    cutoff_high <- p_high
  }

  lvl <- character(length(nnz))
  is_low <- (nnz == 0) | (cutoff_low > 0.05)
  lvl[is_low] <- "low"

  rest <- !is_low
  if (any(rest)){
    tied <- rest & (n05 == n1)
    lvl[tied] <- "not.significant"
    rest2 <- rest & !tied
    if (any(rest2)){
      med <- rest2 & (cutoff_high < 0.05)
      lvl[med]            <- "medium"
      lvl[rest2 & !med]   <- "high"
    }
  }

  ## --- new, additive: force "not.significant" for populations too small
  ## to support a confident call. OFF by default (min_n_for_call = 0L),
  ## so default output is unchanged. Previously the ONLY way to get
  ## "not.significant" was an exact n.0.5 == n.1 tie -- e.g. one 0.5-cell
  ## and one 1-cell out of two total non-zero cells got the same
  ## "not.significant" status as a genuine coin-flip, but two non-zero
  ## cells that happened not to tie got a full-confidence "high"/"medium"
  ## call. This makes the *minimum evidence* required explicit instead.
  if (min_n_for_call > 0L) {
    too_small <- (ntot < min_n_for_call) | (nnz > 0 & nnz < min_n_for_call)
    lvl[too_small] <- "not.significant"
  }

  tbl$expr.level <- lvl
  tbl$p_low  <- p_low
  tbl$p_high <- p_high
  if (identical(fdr_method, "BH")) {
    tbl$q_low  <- q_low
    tbl$q_high <- q_high
  }
  if (isTRUE(effect_size)) {
    tbl$effect_size_low  <- (nnz / ntot) - prob_low
    tbl$effect_size_high <- ifelse(nnz > 0, (n05 / nnz) - prob_high, NA_real_)
  }
  tbl
}

#' Re-adjust p-values for FDR control across an entire discovery scan
#'
#' `discrete.expr.levels()` (optionally with `fdr_method = "BH"`) controls
#' the false discovery rate *within one population's gene panel*. A full
#' biomarker-discovery run tests many (protein x cell type x disease)
#' combinations -- e.g. the manuscript's disease panel applies scProtXpress
#' to 701 proteins across multiple cell types and diseases, which is
#' thousands of individual tests, each with its own within-table BH
#' correction (or none, by default) but nothing correcting *across* tables.
#' This function takes the list of per-population tables produced across a
#' run (e.g. every `sig_tbl.txt` written to `sig_frames/`) and re-adjusts
#' `p_low`/`p_high` jointly across all of them, adding `q_low_global`/
#' `q_high_global` and an `expr.level_global` column that reflects the
#' scan-wide FDR rather than the single-population FDR.
#'
#' @param sig_tbl_list Named list of data.frames as returned by
#'   `discrete.expr.levels()` (each must have `gene`, `p_low`, `p_high`,
#'   `n.0`, `n.0.5`, `n.1`); names are used as a `population` label in the
#'   combined output.
#' @param method `p.adjust` method, default `"BH"`.
#' @return A single combined data.frame (one row per gene x population)
#'   with `population`, the original columns, `q_low_global`,
#'   `q_high_global`, and `expr.level_global`.
#' @export
correct_across_scan <- function(sig_tbl_list, method = "BH"){
  if (length(sig_tbl_list) == 0) stop("correct_across_scan: sig_tbl_list is empty.")
  missing_cols <- vapply(sig_tbl_list, function(x) !all(c("p_low", "p_high") %in% colnames(x)), logical(1))
  if (any(missing_cols)) {
    stop("correct_across_scan: every table must have p_low/p_high columns ",
         "(i.e. must come from a discrete.expr.levels() call, which reports p_low/p_high ",
         "alongside expr.level). ",
         "Offending tables: ", paste(names(sig_tbl_list)[missing_cols], collapse = ", "))
  }
  if (is.null(names(sig_tbl_list)) || any(names(sig_tbl_list) == "")) {
    names(sig_tbl_list) <- paste0("population_", seq_along(sig_tbl_list))
  }
  combined <- do.call(rbind, lapply(names(sig_tbl_list), function(nm) {
    df <- sig_tbl_list[[nm]]
    df$population <- nm
    df
  }))
  combined$q_low_global  <- stats::p.adjust(combined$p_low,  method = method)
  combined$q_high_global <- stats::p.adjust(combined$p_high, method = method)

  ## Recompute expr.level under the global correction, preserving the same
  ## decision cascade as discrete.expr.levels() (low -> tie -> medium/high).
  nnz <- combined$n.0.5 + combined$n.1
  n05 <- combined$n.0.5; n1 <- combined$n.1
  lvl <- character(nrow(combined))
  is_low <- (nnz == 0) | (combined$q_low_global > 0.05)
  lvl[is_low] <- "low"
  rest <- !is_low
  tied <- rest & (n05 == n1)
  lvl[tied] <- "not.significant"
  rest2 <- rest & !tied
  med <- rest2 & (combined$q_high_global < 0.05)
  lvl[med] <- "medium"
  lvl[rest2 & !med] <- "high"
  combined$expr.level_global <- lvl
  combined
}

#' Area under the ROC curve via the Mann-Whitney U statistic
#'
#' Exact, dependency-free AUC (no `pROC`/`ROCR` import needed): the
#' well-known identity AUC = (sum of positive-class ranks - n_pos*(n_pos+1)/2)
#' / (n_pos * n_neg). Ties are handled by `rank()`'s default averaging.
#' @param scores Numeric predicted scores/probabilities.
#' @param labels Logical or 0/1 vector, same length as `scores`; `TRUE`/`1`
#'   is the positive class.
#' @return Numeric AUC in [0, 1], or `NA` if either class is empty.
#' @noRd
.compute_auc <- function(scores, labels){
  labels <- as.logical(labels)
  pos <- scores[labels]; neg <- scores[!labels]
  n_pos <- length(pos); n_neg <- length(neg)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  r_pos <- sum(r[seq_len(n_pos)])
  (r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' Cross-validated patient-level classifier
#'
#' A minimal proof-of-concept classifier for the clinical-translation case:
#' given per-patient summary features (e.g. the proportion of biomarker-
#' positive cells per patient, from `discrete.expr.levels_donor()`'s
#' by-donor breakdown, or any other per-patient numeric summary you
#' compute) and a binary disease label, fits a logistic regression with
#' leave-one-out (default -- appropriate for the small cohorts typical of
#' these studies, e.g. n=28/29) or k-fold cross-validation, and reports
#' out-of-fold AUC, sensitivity, and specificity. This is explicitly a
#' proof-of-concept tool, not a validated diagnostic -- report it as such.
#'
#' @section Small-sample AUC bias -- read before using:
#' Cross-validated AUC computed by pooling out-of-fold predictions (the only
#' option for leave-one-out, where each fold has one held-out point) is
#' NOT centered on 0.5 under a true null for small, weakly- or non-
#' predictive cohorts -- it is biased downward, and the bias gets WORSE as
#' the method gets closer to leave-one-out. For example, at n=40 with a
#' single pure-noise predictor, mean AUC across 60 null draws was 0.33 for
#' LOO, 0.40 for 10-fold, and 0.45 for 3-fold -- all below the naively-
#' expected 0.5. This is a property of small-sample cross-validated AUC
#' generally, not a bug specific to this implementation. **Do not compare
#' an observed AUC to 0.5 as the null expectation for a small cohort.**
#' Use `n_permutations > 0` (below) to get the actual null distribution
#' FOR YOUR data's sample size and CV settings, and report the resulting
#' permutation p-value instead of (or alongside) the raw AUC.
#'
#' @param patient_features data.frame or matrix, one row per patient, one
#'   or more numeric predictor columns (e.g. per-patient proportions of one
#'   or more candidate biomarkers).
#' @param labels Vector aligned to the rows of `patient_features`: logical,
#'   0/1, or a two-level factor/character (the second sorted level, or
#'   `TRUE`/`1`, is treated as the positive/disease class).
#' @param cv_method "loo" (default; leave-one-patient-out) or "kfold".
#' @param k Number of folds when `cv_method = "kfold"` (default 5; ignored
#'   for "loo"). Smaller `k` shows less of the small-sample AUC bias
#'   described above, at the cost of a coarser model per fold.
#' @param threshold Probability cutoff for sensitivity/specificity/the
#'   confusion matrix (default 0.5). AUC does not depend on this.
#' @param fixseed Seed for the k-fold split and for permutations.
#' @param n_permutations Integer, default 0 (off). When > 0, repeats the
#'   ENTIRE cross-validation procedure that many times with `labels`
#'   randomly shuffled each time, building an empirical null distribution
#'   of AUC for your exact sample size, feature set, and CV settings --
#'   the correct way to judge whether the observed AUC is unusual, given
#'   the bias described above. 200 is a reasonable default to request for
#'   a real analysis; runtime scales linearly with this.
#' @return list(predictions, labels, auc, sensitivity, specificity,
#'   confusion (named TP/TN/FP/FN), cv_method, n_patients,
#'   n_excluded_folds, null_aucs, perm_p_value). `n_excluded_folds` counts
#'   folds where the training split had only one class present (degenerate
#'   for logistic regression); those held-out patients get a fallback
#'   prediction equal to the training-fold positive rate rather than a
#'   model-based one, and are flagged rather than silently included as if
#'   they were fully modelled. `null_aucs`/`perm_p_value` are `NULL`/`NA`
#'   unless `n_permutations > 0`; `perm_p_value` is the (1 + #null >=
#'   observed) / (1 + n_permutations) one-sided empirical p-value.
#' @export
fit_patient_classifier <- function(patient_features, labels,
                                   cv_method = c("loo", "kfold"),
                                   k = 5, threshold = 0.5, fixseed = 1234L,
                                   n_permutations = 0L){
  cv_method <- match.arg(cv_method)
  patient_features <- as.data.frame(patient_features)
  n <- nrow(patient_features)
  if (length(labels) != n) {
    stop("fit_patient_classifier: length(labels) must equal nrow(patient_features).")
  }
  if (n < 4) {
    stop("fit_patient_classifier: need at least 4 patients for any form of cross-validation; got ", n, ".")
  }
  if (!is.numeric(n_permutations) || n_permutations < 0) {
    stop("fit_patient_classifier: n_permutations must be a non-negative integer.")
  }
  ## Coerce labels to 0/1 with a clear, reproducible rule: logical as-is;
  ## numeric must already be 0/1; factor/character uses the alphabetically
  ## LATER of the two levels as positive (documented, not left implicit).
  if (is.logical(labels)) {
    y <- as.integer(labels)
  } else if (is.numeric(labels)) {
    if (!all(labels %in% c(0, 1))) stop("fit_patient_classifier: numeric labels must be 0/1.")
    y <- as.integer(labels)
  } else {
    lv <- sort(unique(as.character(labels)))
    if (length(lv) != 2) stop("fit_patient_classifier: labels must have exactly two distinct values; found ",
                              length(lv), ": ", paste(lv, collapse = ", "))
    y <- as.integer(as.character(labels) == lv[2])
    message("fit_patient_classifier: treating '", lv[2], "' as the positive class (vs. '", lv[1], "').")
  }
  if (length(unique(y)) != 2) stop("fit_patient_classifier: labels must have exactly two classes present.")
  if (ncol(patient_features) == 0) stop("fit_patient_classifier: patient_features has no predictor columns.")
  if (identical(cv_method, "kfold") && (k < 2 || k > n)) {
    stop("fit_patient_classifier: k must be between 2 and n_patients.")
  }

  ## Runs the full CV procedure once for a given label vector, returning
  ## out-of-fold predictions. Used for both the real labels and (if
  ## n_permutations > 0) each permutation -- kept as one function so the
  ## real analysis and its null calibration can never silently diverge in
  ## how they're computed.
  .run_cv_once <- function(y_vec, seed_offset){
    fold_id <- if (identical(cv_method, "loo")) {
      seq_len(n)
    } else {
      set.seed(fixseed + seed_offset)
      sample(rep(seq_len(k), length.out = n))
    }
    pred <- rep(NA_real_, n)
    n_excl <- 0L
    for (fold in unique(fold_id)) {
      test_idx  <- which(fold_id == fold)
      train_idx <- setdiff(seq_len(n), test_idx)
      y_train <- y_vec[train_idx]

      if (length(unique(y_train)) < 2) {
        pred[test_idx] <- mean(y_train)
        n_excl <- n_excl + length(test_idx)
        next
      }
      train_df <- patient_features[train_idx, , drop = FALSE]
      train_df$.y <- y_train
      fit <- tryCatch(
        suppressWarnings(stats::glm(.y ~ ., data = train_df, family = stats::binomial())),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        pred[test_idx] <- mean(y_train)
        n_excl <- n_excl + length(test_idx)
        next
      }
      pred[test_idx] <- suppressWarnings(
        stats::predict(fit, newdata = patient_features[test_idx, , drop = FALSE], type = "response")
      )
    }
    list(pred = pred, n_excluded = n_excl)
  }

  real_run <- .run_cv_once(y, seed_offset = 0L)
  pred <- real_run$pred
  n_excluded_folds <- real_run$n_excluded

  auc <- .compute_auc(pred, y)
  pred_class <- as.integer(pred >= threshold)
  TP <- sum(pred_class == 1 & y == 1); TN <- sum(pred_class == 0 & y == 0)
  FP <- sum(pred_class == 1 & y == 0); FN <- sum(pred_class == 0 & y == 1)
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_

  if (n_excluded_folds > 0) {
    warning("fit_patient_classifier: ", n_excluded_folds, "/", n, " held-out patients came from a ",
            "fold whose training split had only one class present, and got a base-rate fallback ",
            "prediction instead of a fitted one. Results are still returned, but with fewer patients ",
            "genuinely cross-validated than n_patients suggests -- consider k-fold with a smaller k ",
            "if this happens often with 'loo'.", call. = FALSE)
  }

  null_aucs <- NULL
  perm_p_value <- NA_real_
  if (n_permutations > 0) {
    set.seed(fixseed)
    perm_labels <- replicate(n_permutations, sample(y), simplify = FALSE)
    null_aucs <- vapply(seq_len(n_permutations), function(i) {
      run <- .run_cv_once(perm_labels[[i]], seed_offset = i)
      .compute_auc(run$pred, perm_labels[[i]])
    }, numeric(1))
    valid_null <- null_aucs[!is.na(null_aucs)]
    perm_p_value <- (1 + sum(valid_null >= auc)) / (1 + length(valid_null))
  }

  list(predictions = pred, labels = y, auc = auc,
      sensitivity = sensitivity, specificity = specificity,
      confusion = c(TP = TP, TN = TN, FP = FP, FN = FN),
      cv_method = cv_method, n_patients = n, n_excluded_folds = n_excluded_folds,
      null_aucs = null_aucs, perm_p_value = perm_p_value)
}


#' Tag markers as "broad" (gating/lineage) or "specific"
#'
#' One-vs-rest differential-abundance tests are structurally disadvantaged
#' for markers that are canonical for MULTIPLE cell types (gating/lineage
#' markers), since the "rest" comparison group is itself partly composed of
#' true-positive cells -- exactly the case scProtXpress's compositional
#' scoring is designed to handle well. This function classifies each marker
#' in a long-format marker-to-cell-type canonical mapping (e.g. from
#' CellMarker DB) as "broad" (canonical for more than one of the profiled
#' cell types) or "specific" (canonical for exactly one), so recall/
#' precision can be reported separately per tier instead of pooled into one
#' number that conflates two different questions.
#'
#' @param marker_celltype_map A long-format data.frame/matrix: one row per
#'   (marker, cell type) canonical pairing, e.g. `data.frame(marker =
#'   c("CD3","CD3","CD19"), cell.type = c("Tcell","NKT","Bcell"))`.
#' @param marker_col,celltype_col Column names (default "marker",
#'   "cell.type").
#' @return data.frame: `marker`, `n_celltypes` (how many distinct cell
#'   types this marker is canonical for in the input), `tier` ("broad" if
#'   `n_celltypes > 1`, else "specific").
#' @export
tag_marker_specificity <- function(marker_celltype_map, marker_col = "marker", celltype_col = "cell.type"){
  if (!all(c(marker_col, celltype_col) %in% colnames(marker_celltype_map))) {
    stop("tag_marker_specificity: marker_celltype_map must have columns '",
         marker_col, "' and '", celltype_col, "'.")
  }
  df <- unique(marker_celltype_map[, c(marker_col, celltype_col)])
  counts <- table(df[[marker_col]])
  out <- data.frame(marker = names(counts), n_celltypes = as.integer(counts),
                    stringsAsFactors = FALSE)
  out$tier <- ifelse(out$n_celltypes > 1, "broad", "specific")
  out[order(out$marker), ]
}

#' Precision/recall/F1/FPR for a set of predicted calls against ground
#' truth, optionally stratified by a grouping column
#'
#' Generic evaluator: takes any table of (marker, cell.type, expr.level)
#' predictions -- e.g. the row-bound `sig_tbl.txt` outputs from a full
#' `run_scProteinXpress()` scan, each tagged with which population it came
#' from -- and a ground-truth set of (marker, cell.type) pairs considered a
#' true positive, and computes the standard classification metrics. Pass
#' `group_col = "tier"` with a column joined in from
#' `tag_marker_specificity()` to compare precision/recall between "broad"
#' and "specific" markers in one call.
#'
#' @param calls data.frame with (at least) `marker`, `cell.type`,
#'   `expr.level` columns.
#' @param ground_truth data.frame with (at least) `marker`, `cell.type`
#'   columns -- pairs where the marker is canonically expected in that cell
#'   type (a true positive).
#' @param positive_levels Character vector of `expr.level` values counted
#'   as a positive call. Default `c("medium", "high")`.
#' @param group_col Optional column name in `calls` to stratify by (e.g.
#'   `"tier"`). When given, one row per group value is returned plus one
#'   `"overall"` row; when `NULL` (default), only the `"overall"` row.
#' @param restrict_to_ground_truth_markers logical, default `TRUE`.
#'   Restricts evaluation to (marker, cell.type) pairs where the MARKER
#'   appears somewhere in `ground_truth` (for any cell type) -- otherwise
#'   every marker scProtXpress called that was simply never assessed for
#'   ground truth (because no compendium is exhaustive) would count as a
#'   false positive purely for not being in an incomplete reference list.
#'   Set `FALSE` only if your `ground_truth` genuinely covers every marker
#'   in `calls`.
#' @return data.frame: `group` ("overall" plus one row per `group_col`
#'   value if given), `n_pairs_evaluated`, `TP`, `FP`, `FN`, `TN`,
#'   `precision`, `recall`, `f1`, `fpr`.
#' @export
evaluate_calls <- function(calls, ground_truth,
                           positive_levels = c("medium", "high"),
                           group_col = NULL,
                           restrict_to_ground_truth_markers = TRUE){
  req_cols <- c("marker", "cell.type", "expr.level")
  if (!all(req_cols %in% colnames(calls))) {
    stop("evaluate_calls: `calls` must have columns: ", paste(req_cols, collapse = ", "))
  }
  if (!all(c("marker", "cell.type") %in% colnames(ground_truth))) {
    stop("evaluate_calls: `ground_truth` must have columns 'marker' and 'cell.type'.")
  }
  if (!is.null(group_col) && !group_col %in% colnames(calls)) {
    stop("evaluate_calls: group_col = '", group_col, "' not found in `calls`.")
  }

  calls <- calls
  if (isTRUE(restrict_to_ground_truth_markers)) {
    calls <- calls[calls$marker %in% unique(ground_truth$marker), , drop = FALSE]
  }
  if (nrow(calls) == 0) {
    stop("evaluate_calls: no rows of `calls` remain after restricting to markers present ",
         "in `ground_truth` -- check that marker names match between the two tables ",
         "(case/whitespace differences are a common cause).")
  }

  gt_pairs <- paste(ground_truth$marker, ground_truth$cell.type, sep = "\u0001")
  calls$is_positive_call <- calls$expr.level %in% positive_levels
  calls$is_true_positive_pair <- paste(calls$marker, calls$cell.type, sep = "\u0001") %in% gt_pairs

  .metrics_for <- function(df, label) {
    TP <- sum(df$is_positive_call & df$is_true_positive_pair)
    FP <- sum(df$is_positive_call & !df$is_true_positive_pair)
    FN <- sum(!df$is_positive_call & df$is_true_positive_pair)
    TN <- sum(!df$is_positive_call & !df$is_true_positive_pair)
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else NA_real_
    fpr <- if ((FP + TN) > 0) FP / (FP + TN) else NA_real_
    data.frame(group = label, n_pairs_evaluated = nrow(df),
              TP = TP, FP = FP, FN = FN, TN = TN,
              precision = precision, recall = recall, f1 = f1, fpr = fpr,
              stringsAsFactors = FALSE)
  }

  out <- .metrics_for(calls, "overall")
  if (!is.null(group_col)) {
    per_group <- do.call(rbind, lapply(sort(unique(calls[[group_col]])), function(g) {
      .metrics_for(calls[calls[[group_col]] == g, , drop = FALSE], g)
    }))
    out <- rbind(per_group, out)
  }
  rownames(out) <- NULL
  out
}


#'
#' Answers "how much reference is enough?" directly: rebuilds thresholds
#' at a series of subsample fractions of the reference (either by CELLS,
#' preserving whichever datasets/cell types happen to be present, or by
#' whole DATASETS, which also tests sensitivity to how many independent
#' studies contribute), repeated `n_replicates` times per fraction to
#' separate sampling noise from a genuine trend, and reports each
#' fraction's threshold estimate against the full reference's.
#'
#' @param data_norm proteins x cells matrix, ADTnorm-normalized.
#' @param metadata data.frame aligned POSITIONALLY to `data_norm`'s
#'   columns; needs a `dataset` column only when `subsample_by =
#'   "datasets"`.
#' @param fractions Numeric vector in (0, 1], default
#'   `c(0.1, 0.25, 0.5, 0.75, 1.0)`. `1.0` is always treated as the full
#'   reference (ground truth to compare every other fraction against) even
#'   if you omit it -- it's added back automatically if missing.
#' @param subsample_by "cells" (default; randomly drop individual cells,
#'   keeping whichever datasets/cell types are already present) or
#'   "datasets" (randomly drop whole datasets -- requires a `dataset`
#'   column in `metadata`, and tests a different question: does the
#'   reference need many independent studies, not just many cells).
#' @param n_replicates Integer, default 3. Independent random subsamples
#'   per fraction (other than 1.0, which has only one "replicate": the
#'   full data).
#' @param markers Optional character vector to restrict to.
#' @param sample_n,num.boot.sample,bootstrap_stratify,threshold_method,ncores,fixseed
#'   Forwarded to every `GetBackground()` call -- keep identical to your
#'   real reference build for an apples-to-apples comparison.
#' @return A data.frame, one row per (marker, fraction, replicate): `marker`,
#'   `fraction`, `replicate`, `n_cells_used`, `Lower`, `Upper`,
#'   `abs_shift_lower`/`abs_shift_upper` (vs. the fraction-1.0 threshold for
#'   that marker). Aggregate with e.g.
#'   `aggregate(abs_shift_upper ~ marker + fraction, data = result, FUN = mean)`
#'   to get a convergence curve per marker, or average across markers for
#'   an overall curve.
#' @export
downsample_reference_stability <- function(
    data_norm, metadata,
    fractions = c(0.1, 0.25, 0.5, 0.75, 1.0),
    subsample_by = c("cells", "datasets"),
    n_replicates = 3,
    markers = NULL,
    sample_n = 100, num.boot.sample = 200,
    bootstrap_stratify = c("none", "cell.type", "cell.type_tissue"),
    threshold_method = c("rectangle", "mixture", "quantile", "zscore"),
    ncores = 1L, fixseed = 1234L
){
  subsample_by       <- match.arg(subsample_by)
  bootstrap_stratify <- match.arg(bootstrap_stratify)
  threshold_method   <- match.arg(threshold_method)

  data_norm <- .as_dense_matrix(data_norm)
  if (ncol(data_norm) != nrow(metadata)) {
    stop("downsample_reference_stability: ncol(data_norm) must equal nrow(metadata) ",
         "(rows of metadata are matched POSITIONALLY to columns of data_norm).")
  }
  if (any(fractions <= 0 | fractions > 1)) {
    stop("downsample_reference_stability: fractions must all be in (0, 1].")
  }
  if (identical(subsample_by, "datasets") && !"dataset" %in% colnames(metadata)) {
    stop("downsample_reference_stability: subsample_by = 'datasets' requires a ",
         "'dataset' column in metadata.")
  }
  if (!1.0 %in% fractions) fractions <- c(fractions, 1.0)
  fractions <- sort(unique(fractions))

  eligible_markers <- if (is.null(markers)) rownames(data_norm) else intersect(markers, rownames(data_norm))
  data_sub <- data_norm[eligible_markers, , drop = FALSE]

  set.seed(fixseed)
  n_cells <- ncol(data_sub)
  dataset_vec <- if ("dataset" %in% colnames(metadata)) as.character(metadata$dataset) else NULL

  .subsample_idx <- function(frac) {
    if (frac >= 1.0) return(seq_len(n_cells))
    if (identical(subsample_by, "cells")) {
      sort(sample(seq_len(n_cells), size = max(1, round(n_cells * frac))))
    } else {
      ds <- unique(dataset_vec)
      keep_ds <- sample(ds, size = max(1, round(length(ds) * frac)))
      which(dataset_vec %in% keep_ds)
    }
  }

  ## Each (fraction, replicate) gets its OWN explicit seed for subsample
  ## selection, rather than relying on the ambient RNG state. This matters
  ## because GetBackground()/RunBootstrap() calls set.seed(fixseed)
  ## internally on every call (so that a single GetBackground() call is
  ## independently reproducible) -- which has the side effect that the RNG
  ## state immediately AFTER any .run_gb() call is fully determined by
  ## fixseed and the call's data size, not by what happened before it. Left
  ## unaddressed, this would silently collapse different "replicates" of
  ## the same fraction into duplicates of each other.
  ## sample.int(.Machine$integer.max, 1) after seeding on a per-cell
  ## deterministic combination of (fraction, replicate) gives a distinct,
  ## reproducible seed for every combination without depending on call order.
  .cell_seed <- function(frac, rep_i, tag) {
    set.seed(fixseed)  # start from a known state so the derived seed is reproducible
    key <- as.integer(round(frac * 1e4)) * 100003L + as.integer(rep_i) * 7L + (if (identical(tag, "subsample")) 1L else 2L)
    set.seed(key)
    sample.int(.Machine$integer.max, 1)
  }

  .run_gb <- function(idx, label, seed_for_gb) {
    GetBackground(data_sub[, idx, drop = FALSE], tempfile(paste0("downsample_", label, "_")),
                  label, metadata[idx, , drop = FALSE],
                  sample_n = sample_n, parallelize = (ncores > 1), ncores = ncores,
                  num.boot.sample = num.boot.sample, fixseed = seed_for_gb,
                  bootstrap_stratify = bootstrap_stratify, threshold_method = threshold_method)
  }

  full_th <- .run_gb(seq_len(n_cells), "full", fixseed)

  rows <- list()
  for (frac in fractions) {
    n_rep_this_frac <- if (frac >= 1.0) 1L else n_replicates
    for (rep_i in seq_len(n_rep_this_frac)) {
      if (frac >= 1.0) {
        idx <- seq_len(n_cells)
        th  <- full_th
      } else {
        set.seed(.cell_seed(frac, rep_i, "subsample"))
        idx <- .subsample_idx(frac)
        th  <- .run_gb(idx, paste0("f", round(frac * 100), "_r", rep_i),
                       .cell_seed(frac, rep_i, "bootstrap"))
      }
      n_used <- length(idx)
      for (m in eligible_markers) {
        if (is.null(th[[m]])) next
        lo <- mean(th[[m]][, "Lower"]); hi <- mean(th[[m]][, "Upper"])
        full_lo <- mean(full_th[[m]][, "Lower"]); full_hi <- mean(full_th[[m]][, "Upper"])
        rows[[length(rows) + 1]] <- data.frame(
          marker = m, fraction = frac, replicate = rep_i, n_cells_used = n_used,
          Lower = lo, Upper = hi,
          abs_shift_lower = abs(lo - full_lo), abs_shift_upper = abs(hi - full_hi),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

#' Per-population sig-table computation, dispatching to the pooled-cell test
#' or a donor-aware one
#'
#' Used by `run_scProteinXpress()`'s `sig.frames` block. Kept as a
#' top-level function (not nested inside `run_scProteinXpress()`) so this
#' dispatch logic -- pooled vs. donor-aware, and the fallback when donor
#' metadata isn't available -- can be unit-tested directly
#' (`tests/testthat/test-sig-test-dispatch.R`) without needing a full
#' pipeline run.
#'
#' @param discretised_data proteins x cells discretised matrix for one
#'   population.
#' @param pop_meta data.frame with (at least) `cell.id` and `sample`
#'   columns for the same population, or `NULL` if unavailable (e.g. a
#'   `Query_<pop>` directory built before donor metadata was persisted
#'   alongside it).
#' @param min.coverage.fraction.sigframe,audit_path As in
#'   `Create.sig.frame()`.
#' @param sig_test_method "pooled" (default; original behaviour, byte-for-
#'   byte identical output) or one of `discrete.expr.levels_donor()`'s
#'   `test_method` options ("donor_t", "donor_wilcox", "design_effect").
#'   Falls back to "pooled" with a warning if `pop_meta` is `NULL` or has
#'   no `sample` column, rather than silently doing nothing.
#' @param min_donors_for_call As in `discrete.expr.levels_donor()`.
#' @return A sig-table data.frame (same shape either way: `discrete.expr.levels()`'s
#'   or `discrete.expr.levels_donor()`'s output).
#' @noRd
.compute_sig_table_for_population <- function(discretised_data, pop_meta,
                                              min.coverage.fraction.sigframe,
                                              audit_path,
                                              sig_test_method = "pooled",
                                              min_donors_for_call = 2L){
  ## Prefer the canonical 'donor' column (present whenever pop_meta was
  ## persisted by a version of run_scProteinXpress() that resolves
  ## donor_col). Fall back to 'sample' for Query_<pop>/*_Metadata.rds files
  ## written before donor_col existed, where donor == sample was implicit
  ## rather than explicit -- same effective behaviour, just not
  ## configurable for those older files without regenerating them.
  donor_col_present <- !is.null(pop_meta) && "cell.id" %in% colnames(pop_meta) &&
    ("donor" %in% colnames(pop_meta) || "sample" %in% colnames(pop_meta))
  donor_col_name <- if (!is.null(pop_meta) && "donor" %in% colnames(pop_meta)) "donor" else "sample"

  if (identical(sig_test_method, "pooled") || !donor_col_present) {
    if (!identical(sig_test_method, "pooled")) {
      warning("sig_test_method = '", sig_test_method, "' requested but no per-cell ",
              "donor/sample metadata is available for this population (a Query_<pop> ",
              "directory built before donor-aware metadata was persisted, or metadata ",
              "with no 'donor'/'sample' column) -- falling back to the pooled-cell test ",
              "for THIS population only. This test treats every cell as an independent ",
              "trial and is anti-conservative whenever cells are correlated within donor; ",
              "treat any resulting call for this population with that in mind.",
              call. = FALSE)
    }
    sig.res <- Create.sig.frame(discretised_data,
                                min.coverage.fraction.sigframe = min.coverage.fraction.sigframe,
                                audit_csv_path = audit_path)
    return(discrete.expr.levels(sig.res[[1]], sig.res[[2]],
                                prob_low = sig.res[[3]], prob_high = sig.res[[4]]))
  }

  ## Donor-aware path. Note: this still derives prob_low/prob_high from the
  ## SAME pooled population (Create.sig.frame()'s self-referential
  ## baseline) -- the donor fix addresses pseudoreplication in the TEST
  ## across cells vs. donors; it does not by itself fix the separate
  ## self-referential-baseline issue (see ComputeExternalBaseline() for
  ## that, which composes with either test).
  sig.res <- Create.sig.frame(discretised_data,
                              min.coverage.fraction.sigframe = min.coverage.fraction.sigframe,
                              audit_csv_path = audit_path)
  donor_ids <- pop_meta[[donor_col_name]][match(colnames(discretised_data), pop_meta$cell.id)]
  if (anyNA(donor_ids)) {
    warning(".compute_sig_table_for_population: could not align donor IDs (column '",
            donor_col_name, "') to every discretised cell for this population (",
            sum(is.na(donor_ids)), "/", length(donor_ids), " unmatched) -- falling ",
            "back to the pooled-cell test for this population.", call. = FALSE)
    return(discrete.expr.levels(sig.res[[1]], sig.res[[2]],
                                prob_low = sig.res[[3]], prob_high = sig.res[[4]]))
  }
  by_donor <- Create.sig.frame_by_donor(discretised_data, donor_ids)
  discrete.expr.levels_donor(by_donor, sig.res[[3]], sig.res[[4]],
                             test_method = sig_test_method,
                             min_donors_for_call = min_donors_for_call)
}
