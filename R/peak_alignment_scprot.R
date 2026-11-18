###############################################################################
# peak_alignment_scprot.R
#
# Local copy of ADTnorm::peak_alignment() with fixed-mode alignment.
#
# (a) Why
#     Query alignment needs to warp each query sample onto the landmarks and
#     the warping grid saved in the reference (<file.name>_bg_prot_pars.rds:
#     `fixed_target_landmark`, `warp_parameters`). The published
#     ADTnorm::peak_alignment() (https://github.com/yezhengSTAT/ADTnorm, v1.0)
#     has neither the `alignment_mode` / `warp_parameters_fixed` arguments
#     nor the `target_landmark_used` / `warp_parameters` outputs, so with it
#     the reference is saved without landmarks and every query marker is left
#     unaligned (arcsinh values only) without any error. The pipeline
#     therefore calls this internal copy, `.peak_alignment_scprot()`, instead.
#
# (b) Changes relative to ADTnorm::peak_alignment() (v1.0)
#     - New arguments `alignment_mode` (stored in the output) and
#       `warp_parameters_fixed`: when given, the B-spline range (`from`,
#       `to`), the clipping bounds (`lower_bound`, `upper_bound`), `lambda`,
#       `nbreaks` and `nb` are taken from the reference instead of being
#       recomputed from the query.
#     - When the last of >= 3 landmarks is above `upper_bound` it is
#       clipped to `upper_bound`, as already done for the valley.
#     - Returns `target_landmark_used`, `alignment_mode`, `marker_min_raw`,
#       `marker_max_raw` and `warp_parameters` next to `cell_x_adt_norm`.
#     Apart from that clip, without `warp_parameters_fixed` (reference
#     building) the alignment is the same as in ADTnorm::peak_alignment().
#
# License: GPL-3 (inherited from ADTnorm; see ADTnorm LICENSE.md).
# Original author: Ye Zheng <yzheng23@fredhutch.org>.
###############################################################################

#' @noRd
.peak_alignment_scprot = function(
    cell_x_adt,
    cell_x_feature = NULL,
    landmark_matrix = NULL,
    target_landmark = NULL,
    alignment_mode = NA_character_,
    warp_parameters_fixed = NULL
) {
  ## get parameters
  grouping = NULL
  monwrd = TRUE
  subsample = NULL
  peakNr = NULL
  clipRange = 0.01
  nbreaks = 11
  bwFac = 2
  warpFuns = FALSE
  chunksinze = 10
  newNcFile = NULL
  z = NULL
  nb = 1001

  ## smoothing / warping penalty
  lambda = 1e-5

  exp_data = cell_x_adt
  cell_x_adt_norm = cell_x_adt
  samples = levels(cell_x_feature$sample)

  ## raw marker range of the CURRENT input
  marker_min_raw = min(cell_x_adt, na.rm = TRUE)
  marker_max_raw = max(cell_x_adt, na.rm = TRUE)

  ## mean landmark positions across samples
  landmark_means_raw = colMeans(landmark_matrix, na.rm = TRUE)

  ## actual target used for alignment
  if (!is.null(target_landmark)) {
    target_landmark_used = target_landmark
  } else {
    target_landmark_used = landmark_means_raw
  }

  ## set up fda parameters
  extend = 0.15

  ## IMPORTANT:
  ## - if warp_parameters_fixed is provided -> use saved REFERENCE grid/range
  ## - otherwise -> keep original behavior (REFERENCE BUILDING unchanged)
  if (!is.null(warp_parameters_fixed)) {
    if (!is.list(warp_parameters_fixed)) {
      stop("warp_parameters_fixed must be a list.")
    }

    needed <- c("lambda", "nbreaks", "nb", "from", "to", "lower_bound", "upper_bound")
    miss <- setdiff(needed, names(warp_parameters_fixed))
    if (length(miss) > 0) {
      stop("warp_parameters_fixed is missing: ", paste(miss, collapse = ", "))
    }

    lambda      = warp_parameters_fixed$lambda
    nbreaks     = warp_parameters_fixed$nbreaks
    nb          = warp_parameters_fixed$nb
    from        = warp_parameters_fixed$from
    to          = warp_parameters_fixed$to
    lower_bound = warp_parameters_fixed$lower_bound
    upper_bound = warp_parameters_fixed$upper_bound

  } else {
    from = min(c(marker_min_raw, target_landmark_used[1], min(landmark_matrix, na.rm = TRUE))) -
      diff(range(cell_x_adt, na.rm = TRUE)) * extend
    to = max(c(marker_max_raw, target_landmark_used[length(target_landmark_used)], max(landmark_matrix, na.rm = TRUE))) +
      diff(range(cell_x_adt, na.rm = TRUE)) * extend

    lower_bound = marker_min_raw
    upper_bound = marker_max_raw
  }

  wbasis = fda::create.bspline.basis(
    rangeval = c(from, to),
    norder = 4,
    breaks = seq(from, to, len = nbreaks)
  )
  Wfd0 = fda::fd(matrix(0, wbasis$nbasis, 1), wbasis)
  WfdPar = fda::fdPar(Wfd0, 1, lambda)

  density_y = c()
  for (sample in samples) {
    cell_ind_tmp = which(cell_x_feature$sample == sample)
    cell_ind = cell_ind_tmp[which(!is.na(cell_x_adt[cell_ind_tmp]))]

    if (length(cell_ind) > 0) {
      density_y = cbind(
        density_y,
        stats::density(cell_x_adt[cell_ind], from = from, to = to, n = nb, na.rm = TRUE)$y
      )
    } else {
      density_y = cbind(density_y, rep(NA, nb))
    }
  }
  colnames(density_y) = samples

  arg_vals = seq(from, to, len = nb)
  fdobj = fda::Data2fd(arg_vals, density_y, wbasis)

  if (ncol(landmark_matrix) == 1) { ## only one landmark: offset
    offsets = landmark_matrix - stats::median(landmark_matrix, na.rm = TRUE)
    names(offsets) = samples

    funs = funsBack = vector("list", nrow(landmark_matrix))
    names(funs) = samples
    names(funsBack) = samples

    for (j in seq_along(funs)) {
      funs[[samples[j]]] = function(x) x - z
      e1 = new.env(hash = TRUE)
      e1$z = offsets[samples[j]]
      environment(funs[[samples[j]]]) = e1

      funsBack[[samples[j]]] = function(x) x + z
      e2 = new.env(hash = TRUE)
      e2$z = offsets[samples[j]]
      environment(funsBack[[samples[j]]]) = e2
    }

  } else { ## more than one landmark: warping

    ## if any valley is beyond the upper bound of range, replace by upper bound
    if (ncol(landmark_matrix) >= 2 && any(landmark_matrix[, 2] > upper_bound, na.rm = TRUE)) {
      landmark_matrix[which(landmark_matrix[, 2] > upper_bound), 2] = upper_bound
      print(paste0(
        "Warning: some valley landmarks are larger than the upper bound of the range. ",
        "They are replaced by the upper bound used for alignment."
      ))
    }

    ## if first landmark goes below lower bound, replace by lower bound
    if (any(landmark_matrix[, 1] < lower_bound, na.rm = TRUE)) {
      landmark_matrix[which(landmark_matrix[, 1] < lower_bound), 1] = lower_bound
      print("Warning: some first landmarks are smaller than the lower bound of the range. They are replaced by the lower bound used for alignment.")
    }

    ## if last landmark goes beyond upper bound, replace by upper bound
    if (ncol(landmark_matrix) >= 3 && any(landmark_matrix[, ncol(landmark_matrix)] > upper_bound, na.rm = TRUE)) {
      landmark_matrix[which(landmark_matrix[, ncol(landmark_matrix)] > upper_bound), ncol(landmark_matrix)] = upper_bound
    }

    args = list(
      "unregfd" = fdobj,
      "fdobj" = fdobj,
      "ximarks" = landmark_matrix,
      "WfdPar" = WfdPar,
      "monwrd" = monwrd,
      "x0marks" = target_landmark_used
    )

    args_run = args[intersect(names(formals(fda::landmarkreg)), names(args))]
    regDens = do.call(fda::landmarkreg, args_run, quote = TRUE)

    warpfdobj = regDens$warpfd
    warpedX = fda::eval.fd(warpfdobj, arg_vals)
    warpedX[1, ] = utils::head(arg_vals, 1)
    warpedX[nrow(warpedX), ] = utils::tail(arg_vals, 1)

    funs = apply(warpedX, 2, stats::approxfun, arg_vals)
    funsBack = apply(warpedX, 2, function(a, b) stats::approxfun(b, a), arg_vals)
  }

  names(funs) = names(funsBack) = samples

  warped_landmark_matrix = landmark_matrix

  ## transform the raw data using the warping functions
  for (i in seq_along(funs)) {
    cell_ind_tmp = which(cell_x_feature$sample == samples[i])
    cell_index = cell_ind_tmp[which(!is.na(cell_x_adt[cell_ind_tmp]))]

    if (length(cell_index) > 0) {
      thisDat = t(t(cell_x_adt[cell_index]))
      newDat = as.matrix(funs[[i]](thisDat))
      newDat[is.na(newDat)] = thisDat[is.na(newDat)]
      cell_x_adt_norm[cell_index] = newDat
      warped_landmark_matrix[i, ] = funs[[i]](landmark_matrix[i, ])
    } else {
      warped_landmark_matrix[i, ] = NA
    }
  }

  return(list(
    # normalized values for this marker
    cell_x_adt_norm = cell_x_adt_norm,

    # raw range of this CURRENT marker input
    marker_min_raw = marker_min_raw,
    marker_max_raw = marker_max_raw,

    # target landmarks actually used
    target_landmark_used = target_landmark_used,

    # alignment mode used for this marker
    alignment_mode = alignment_mode,

    # compact warp settings needed for reproducibility / later fixed-mode processing
    warp_parameters = list(
      monwrd = monwrd,
      lambda = lambda,
      nbreaks = nbreaks,
      nb = nb,
      from = from,
      to = to,
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
  ))
}

###############################################################################
# .peak_alignment_by_sample(): warp of one marker that does not lose the whole
# marker when fda::landmarkreg fails for some sample (bug B).
#
# The whole call is tried first, as before. If it fails, every sample is
# warped on its own with the settings the whole call would have used
# (B-spline range, clipping bounds, lambda, grid; or `warp_parameters_fixed`
# when given). The warp of a sample does not depend on the other samples, so
# the samples that do not fail get the same values as in a successful whole
# call. Samples that still fail are set to NA and listed in `warp_failed`; if
# every sample fails, the call fails (the marker is then not aligned).
###############################################################################

#' @noRd
.peak_alignment_by_sample <- function(cell_x_adt,
                                      cell_x_feature = NULL,
                                      landmark_matrix = NULL,
                                      target_landmark = NULL,
                                      alignment_mode = NA_character_,
                                      warp_parameters_fixed = NULL) {
  res <- tryCatch(.peak_alignment_scprot(cell_x_adt, cell_x_feature, landmark_matrix, target_landmark,
                                         alignment_mode, warp_parameters_fixed),
                  error = function(e) e)
  if (!inherits(res, "error")) {
    res$warp_failed <- character(0)
    return(res)
  }

  x  <- cell_x_adt
  tl <- if (!is.null(target_landmark)) target_landmark else colMeans(landmark_matrix, na.rm = TRUE)
  wp <- warp_parameters_fixed
  if (is.null(wp)) {
    ## as computed by .peak_alignment_scprot() for the whole call
    rng <- diff(range(x, na.rm = TRUE))
    wp <- list(monwrd = TRUE, lambda = 1e-5, nbreaks = 11, nb = 1001,
               from = min(c(min(x, na.rm = TRUE), tl[1], min(landmark_matrix, na.rm = TRUE))) - rng * 0.15,
               to   = max(c(max(x, na.rm = TRUE), tl[length(tl)], max(landmark_matrix, na.rm = TRUE))) + rng * 0.15,
               lower_bound = min(x, na.rm = TRUE), upper_bound = max(x, na.rm = TRUE))
  }
  samples <- levels(cell_x_feature$sample)
  lm <- as.matrix(landmark_matrix)
  if (is.null(rownames(lm))) rownames(lm) <- samples
  out <- x
  failed <- character(0)
  n_tried <- 0L
  for (s in samples) {
    i <- which(cell_x_feature$sample == s)
    if (!length(i)) next
    n_tried <- n_tried + 1L
    r <- tryCatch(suppressWarnings(.peak_alignment_scprot(
      x[i], data.frame(sample = factor(rep(s, length(i)), levels = s)),
      lm[s, , drop = FALSE], tl, alignment_mode, wp)), error = function(e) NULL)
    if (is.null(r)) {
      out[i] <- NA_real_
      failed <- c(failed, s)
    } else {
      out[i] <- r$cell_x_adt_norm
    }
  }
  if (length(failed) == n_tried) {
    stop("the warp failed for every sample (", conditionMessage(res), ")")
  }
  list(cell_x_adt_norm = out, marker_min_raw = min(x, na.rm = TRUE), marker_max_raw = max(x, na.rm = TRUE),
       target_landmark_used = tl, alignment_mode = alignment_mode, warp_parameters = wp,
       warp_failed = failed)
}
