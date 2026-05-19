#' Run the unified scProtXpress pipeline (CyTOF or CITE-seq)
#'
#' Builds a reference background (peak alignment + bootstrap thresholds),
#' optionally discretises a query against it (per-cell-type and global),
#' and emits per-population signature frames with a binomial expression-
#' level call. The same entry point handles both modalities; only the
#' ingest step branches on `modality`.
#'
#' @section Modality dispatch:
#' `modality` controls three things:
#' \enumerate{
#'   \item Ingest helper: `.ingest_cytof()` returns a dense base matrix;
#'         `.ingest_citeseq()` returns a sparse dgCMatrix.
#'   \item `qc_outlier_cells = "auto"` activates `scuttle::isOutlier`
#'         only when `modality == "CITE-seq"`; "yes" / "no" force the
#'         choice regardless of modality.
#'   \item `num.cells.sample = NULL` resolves to 1000L for CyTOF or
#'         100L for CITE-seq, matching the historical defaults.
#' }
#'
#' @section Metadata column resolution (A1):
#' Each `*_col` argument is resolved against `metadata` by trying the
#' explicit value first then a list of fallbacks:
#' \describe{
#'   \item{cell_id_col}{user, "cell_id", "cell.id".}
#'   \item{cell_type_col}{user, "cell_type", "cell.type",
#'         "cell.type_harmonised".}
#'   \item{dataset_col}{user, "dataset_id", "dataset".}
#'   \item{sample_col}{user, "sample_id", "sample".}
#'   \item{batch_col}{user, "batch_id", "batch", falls back to `sample`.}
#'   \item{tissue_col}{user, "tissue", "tissue_id", falls back to
#'         `tissue_default`.}
#' }
#' The first column found is renamed in place to the internal canonical
#' name used by everything downstream (`cell.id`, `cell.type`, `dataset`,
#' `sample`, `batch`, `tissue`).
#'
#' @section Input transformation state:
#' `is_transformed = NULL` (default) delegates detection to
#' `.detect_input_mode()` (autodetect by integer-coercibility).
#' Explicit TRUE / FALSE overrides detection (TRUE = "already
#' transformed, skip arcsinh"; FALSE = "force apply arcsinh").
#'
#' @param data Proteins x cells matrix / Matrix / data.frame. Driver helpers
#'   `unpack_sce_cytof()` and `unpack_list_citeseq()` extract this from
#'   modality-specific containers.
#' @param metadata data.frame with rows = cells; column names match the
#'   user-supplied `*_col` arguments (or one of the documented fallbacks).
#' @param org reserved; not used by the current pipeline.
#' @param dir Output root.
#' @param file.name Run identifier (used as `Ref_<file.name>` folder name).
#' @param precBack If a valid directory, load Ref outputs from there
#'   instead of building.
#' @param discretize logical. Run `Discretise_data()` after Ref/Query.
#' @param sig.frames logical. Run `Create.sig.frame()` + `discrete.expr.levels()`.
#' @param num.boot.sample integer. Bootstrap iterations.
#' @param num.cells.sample integer or NULL. Cells per bootstrap sample.
#'   NULL resolves to 1000L (CyTOF) or 100L (CITE-seq).
#' @param ncores integer >= 1. Workers for the inner `foreach %dopar%` loops.
#' @param fixseed integer. Seed for reproducibility.
#' @param tissue_default character(1). Used when metadata has no `tissue`.
#' @param modality "CyTOF" or "CITE-seq". Selects the ingest helper and
#'   resolves auto-defaults for `qc_outlier_cells` and `num.cells.sample`.
#' @param cell_id_col,cell_type_col,dataset_col,sample_col,batch_col,tissue_col
#'   Column names in `metadata` (see Metadata column resolution).
#' @param is_transformed NULL, TRUE, FALSE. See Input transformation state.
#' @param sqs_threshold Numeric. SQS cutoff when `apply_sqs != "none"`.
#' @param apply_sqs "none" (default), "ref_only", "both".
#' @param min.cells.feature integer. Minimum non-zero cells per marker.
#' @param min.coverage.fraction Numeric in [0, 1]. Minimum coverage fraction
#'   per marker.
#' @param min.coverage.fraction.sigframe Numeric in [0, 1]. Per-population
#'   coverage filter for `Create.sig.frame()` baseline.
#' @param qc_outlier_cells "auto" (default), "yes", "no".
#' @param qc_nmads integer. Used by the scuttle filter when active.
#' @param do_subsample logical. Apply hierarchical subsampling of the Ref.
#' @param target_ct integer. Per-cell-type cap for hierarchical subsampling.
#' @param bootstrap_stratify one of "none" (default), "cell.type",
#'   "cell.type_tissue". Controls how cells are pooled when generating each
#'   per-marker bootstrap iteration:
#'   \itemize{
#'     \item `"none"` reproduces the original `scGeneXpress` bootstrap:
#'       per-marker sampling from the global pool of non-zero values across
#'       all cells, with replacement. Default; recommended for methodological
#'       fidelity with the scGeneXpress family.
#'     \item `"cell.type"` stratifies the bootstrap by cell type. Each cell
#'       type contributes an equal share of `num.cells.sample` cells per
#'       iteration. Advanced option.
#'     \item `"cell.type_tissue"` stratifies by `cell.type` x `tissue` and
#'       reproduces the behaviour shipped in Hito 2 T-3, kept for
#'       back-compatibility with the CyTOF reference run. Falls back to
#'       `"cell.type"` (with a warning) when no `tissue` column is present.
#'   }
#' @param ... Forwarded to internal normalisation helpers.
#' @return `invisible(NULL)`. Side effects: files under `dir/Ref_<file.name>/`.
#' @export
run_scProteinXpress <- function(
    data,
    metadata,
    org = "human",
    dir,
    file.name,
    precBack = NA,
    discretize = TRUE,
    sig.frames = TRUE,
    num.boot.sample = 1000,
    num.cells.sample = NULL,        # NULL => resolves by modality (see header).
    ncores = 1L,
    fixseed = 1234,
    tissue_default = "PBMC",
    ## ---- modality dispatch ----
    modality = c("CyTOF", "CITE-seq"),
    ## ---- metadata column mapping (A1 fallbacks applied internally) ----
    cell_id_col   = "cell_id",
    cell_type_col = "cell_type",
    dataset_col   = "dataset_id",
    sample_col    = "sample_id",
    batch_col     = "batch",
    tissue_col    = "tissue",
    ## ---- input transformation state ----
    is_transformed = NULL,          # NULL => autodetect by .detect_input_mode().
    ## ---- SQS (off by default) ----
    sqs_threshold = 1.0,
    apply_sqs = c("none", "ref_only", "both"),
    ## ---- P2: COVERAGE FILTERS (Mode C: single global filter after union) ----
    min.cells.feature     = 10,     # absolute: at least N cells with non-zero, non-NA value
    min.coverage.fraction = 0.10,   # fraction:  at least 10% of cells have a non-NA measurement
    ## ---- P3: SIG_FRAME PER-POPULATION COVERAGE FILTER (binomial baseline) ----
    min.coverage.fraction.sigframe = 0.10,
    ## ---- QC: cell-level outlier filter (scuttle, opt-in unified) ----
    qc_outlier_cells = c("auto", "yes", "no"),
    qc_nmads = 3,
    ## ---- Hierarchical subsampling of the reference background ----
    do_subsample = TRUE,
    target_ct    = 2000L,
    ## ---- Bootstrap stratification (T-7 revert; default = scGeneXpress) ----
    bootstrap_stratify = c("none", "cell.type", "cell.type_tissue"),
    ...){

  modality           <- match.arg(modality)
  apply_sqs          <- match.arg(apply_sqs)
  qc_outlier_cells   <- match.arg(qc_outlier_cells)
  bootstrap_stratify <- match.arg(bootstrap_stratify)

  ## Resolve modality-dependent default for num.cells.sample (A4).
  if (is.null(num.cells.sample)) {
    num.cells.sample <- if (identical(modality, "CyTOF")) 1000L else 100L
  }
  ## is_transformed is forwarded as-is (NULL / TRUE / FALSE) to
  ## .detect_input_mode(), which interprets:
  ##   NULL  => autodetect by integer-coercibility heuristic;
  ##   TRUE  => force-skip arcsinh;
  ##   FALSE => force-apply arcsinh.
  ## See the .detect_input_mode() docstring in internalFunctions.R.

  ## Resolve `qc_outlier_cells = "auto"` to a concrete "yes" / "no". This is
  ## the single source of truth for the scuttle outlier filter activation;
  ## downstream code (cleanData_ADT, PrepQueryData_ADT) only ever sees the
  ## resolved value (Hito 3, Task #14 follow-up to Hito 2 D2).
  qc_outlier_cells_resolved <- if (identical(qc_outlier_cells, "auto")) {
    if (identical(modality, "CITE-seq")) "yes" else "no"
  } else {
    qc_outlier_cells
  }

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  if (base::endsWith(dir, "/")) dir <- base::substr(dir, 1, base::nchar(dir) - 1)
  use_precBack <- is.character(precBack) && length(precBack) == 1 && nzchar(precBack) && dir.exists(precBack)

  ## -------------------------------------------------------------------------
  ## Entry-point validations (V1-V3, V14-V25). Cell-level validations that
  ## need post-translation metadata (V4-V13, V28) happen after the metadata
  ## translation block below.
  ## -------------------------------------------------------------------------
  if (!(is.matrix(data) || inherits(data, "Matrix") || is.data.frame(data))) {
    stop("data must be matrix, Matrix, or data.frame; got ", paste(class(data), collapse = ","))
  }
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame; got ", paste(class(metadata), collapse = ","))
  if (!is.null(is_transformed) && !(isTRUE(is_transformed) || isFALSE(is_transformed))) {
    stop("is_transformed must be NULL, TRUE, or FALSE")
  }
  if (!is.numeric(ncores) || length(ncores) != 1L || ncores < 1L) stop("ncores must be a single integer >= 1")
  if (!is.numeric(target_ct) || length(target_ct) != 1L || target_ct < 1L) stop("target_ct must be a single integer >= 1")
  if (!is.numeric(min.coverage.fraction) || min.coverage.fraction < 0 || min.coverage.fraction > 1) {
    stop("min.coverage.fraction must be in [0, 1]")
  }
  if (!is.numeric(min.coverage.fraction.sigframe) || min.coverage.fraction.sigframe < 0 || min.coverage.fraction.sigframe > 1) {
    stop("min.coverage.fraction.sigframe must be in [0, 1]")
  }
  if (!is.numeric(min.cells.feature) || min.cells.feature < 1L) stop("min.cells.feature must be >= 1")
  if (!nzchar(file.name)) stop("file.name must be a non-empty string")
  if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("cannot create output dir: ", dir)
  }
  if (qc_outlier_cells != "no") {
    ## scuttle is only an Imports for CITE-seq usage; defer the requirement
    ## to the point where it would actually be invoked.
    will_filter <- identical(qc_outlier_cells, "yes") ||
                   (identical(qc_outlier_cells, "auto") && identical(modality, "CITE-seq"))
    if (will_filter && !requireNamespace("scuttle", quietly = TRUE)) {
      stop("qc_outlier_cells would activate scuttle::isOutlier; install scuttle or set qc_outlier_cells = 'no'.")
    }
  }
  if (ncores >= 2L) {
    if (!requireNamespace("foreach", quietly = TRUE) || !requireNamespace("doParallel", quietly = TRUE)) {
      stop("parallel execution requires foreach + doParallel")
    }
  }
  if (use_precBack) {
    pb_pars <- list.files(precBack, pattern = "_bg_prot_pars\\.rds$", full.names = TRUE)
    pb_th   <- list.files(precBack, pattern = "_th_dist\\.rds$",      full.names = TRUE)
    if (length(pb_pars) != 1L || length(pb_th) != 1L) {
      stop(sprintf("precBack does not contain a single *_bg_prot_pars.rds and *_th_dist.rds: %s", precBack))
    }
  }

  ## -------------------------------------------------------------------------
  ## Modality dispatch: ingest stage.
  ## -------------------------------------------------------------------------
  ingest <- switch(modality,
    "CyTOF"    = .ingest_cytof(data, metadata),
    "CITE-seq" = .ingest_citeseq(data, metadata)
  )
  data     <- ingest$data
  metadata <- ingest$metadata

  ## -------------------------------------------------------------------------
  ## Metadata translation with A1 fallbacks: each *_col is resolved against
  ## metadata trying the user value first then a list of fallbacks. The
  ## first match wins and is renamed in place to the canonical internal
  ## name. Internal canonical names are: cell.id, cell.type, dataset,
  ## sample, batch, tissue.
  ## -------------------------------------------------------------------------
  metadata <- .resolve_meta_col(metadata, cell_id_col,
                                fallbacks = c("cell_id", "cell.id"),
                                internal_name = "cell.id")
  metadata <- .resolve_meta_col(metadata, cell_type_col,
                                fallbacks = c("cell_type", "cell.type", "cell.type_harmonised"),
                                internal_name = "cell.type")
  metadata <- .resolve_meta_col(metadata, dataset_col,
                                fallbacks = c("dataset_id", "dataset"),
                                internal_name = "dataset")
  metadata <- .resolve_meta_col(metadata, sample_col,
                                fallbacks = c("sample_id", "sample"),
                                internal_name = "sample")
  ## `batch` is optional in the published ADTnorm: if absent, falls back
  ## to `sample` (matches the codebase/run_scProteinXpress.R behaviour).
  ## Driver-level callers (CyTOF) always pass an explicit batch, so this
  ## fallback only fires for callers that genuinely lack it.
  metadata <- .resolve_meta_col(metadata, batch_col,
                                fallbacks = c("batch_id", "batch"),
                                internal_name = "batch",
                                fallback_value = NULL,
                                fallback_from  = "sample")
  ## tissue is optional: missing column or empty cells fall back to tissue_default.
  metadata <- .resolve_meta_col(metadata, tissue_col,
                                fallbacks = c("tissue", "tissue_id"),
                                internal_name = "tissue",
                                fallback_value = tissue_default)

  metadata$cell.id   <- as.character(metadata$cell.id)
  metadata$cell.type <- as.character(metadata$cell.type)
  metadata$dataset   <- as.character(metadata$dataset)
  metadata$sample    <- as.character(metadata$sample)
  metadata$batch     <- as.character(metadata$batch)
  metadata$tissue    <- as.character(metadata$tissue)
  metadata$batch[is.na(metadata$batch) | metadata$batch == ""] <-
    metadata$sample[is.na(metadata$batch) | metadata$batch == ""]
  metadata$tissue[is.na(metadata$tissue) | metadata$tissue == ""] <- tissue_default

  rownames(metadata) <- metadata$cell.id

  ## -------------------------------------------------------------------------
  ## Cell-level alignment validations (A5 strengthening of V13).
  ## -------------------------------------------------------------------------
  if (ncol(data) < 50L) stop(sprintf("data has %d cells; at least 50 are required (published cell-type threshold).", ncol(data)))
  if (any(duplicated(colnames(data)))) stop("duplicate cell ids in colnames(data)")
  if (any(duplicated(metadata$cell.id))) stop("duplicate cell ids in metadata$cell.id")
  if (any(is.na(metadata$cell.id) | metadata$cell.id == "")) stop("metadata$cell.id has empty/NA values")
  if (any(duplicated(rownames(data)))) stop("duplicate marker names in rownames(data)")

  meta_only <- setdiff(metadata$cell.id, colnames(data))
  data_only <- setdiff(colnames(data), metadata$cell.id)
  if (length(meta_only) > 0L || length(data_only) > 0L) {
    stop(sprintf(
      "data and metadata are not aligned: %d cell.id only in metadata (e.g. %s), %d cell ids only in data (e.g. %s). Align inputs before calling run_scProteinXpress().",
      length(meta_only), paste(utils::head(meta_only, 5), collapse = ", "),
      length(data_only), paste(utils::head(data_only, 5), collapse = ", ")
    ))
  }
  metadata <- metadata[match(colnames(data), metadata$cell.id), , drop = FALSE]
  ## Defensive contract: post-conditions of the ingest+translation block.
  stopifnot(all(colnames(data) == metadata$cell.id))
  
  ref_dir <- file.path(dir, paste0("Ref_", file.name))
  
  #############################################################################
  ################### BACKGROUND CONSTRUCTION / SELECTION ######################
  #############################################################################
  
  if (!use_precBack) {
    cell.type.freq <- as.data.frame(base::table(metadata$cell.type))
    cell.type.to.keep <- base::as.character(
      cell.type.freq[which(cell.type.freq$Freq >= 50), "Var1"]
    )
    
    metadataBack <- metadata[base::which(metadata$cell.type %in% cell.type.to.keep), , drop = FALSE]
    dataBack <- data[, base::as.character(metadataBack$cell.id), drop = FALSE]
    
    if (base::length(cell.type.to.keep) < 2) {
      stop("Data contains less than 2 celltypes with at least 50 cells. Background reference cannot be created.")
    }
    
    if (!dir.exists(ref_dir)) dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
    
    bg_norm_path <- file.path(ref_dir, paste0(file.name, "_NormData.rds"))
    bg_pars_path <- file.path(ref_dir, paste0(file.name, "_bg_prot_pars.rds"))
    th_path      <- file.path(ref_dir, paste0(file.name, "_th_dist.rds"))
    # Sanitised raw counts + metadata used to build the background. Persisted so
    # that the Ref self-query (when precBack=NA) reuses *exactly* the same cells
    # and feature set, including NA-padding for markers absent in some datasets.
    bg_counts_path <- file.path(ref_dir, paste0(file.name, "_BackCounts.rds"))
    bg_meta_path   <- file.path(ref_dir, paste0(file.name, "_BackMetadata.rds"))
    
    if (!(file.exists(bg_norm_path) && file.exists(bg_pars_path))) {
      message("Building background...")
      
      # P2 (Mode C): build a single NA-padded union matrix without per-dataset
      # filtering, then apply ONE global cleanData_ADT with the configured
      # coverage thresholds. This avoids the per-dataset bias where a marker
      # measured in only one tiny panel could pass a per-dataset filter and
      # later contaminate global statistics with a 5%-coverage padded matrix.
      cells_by_ds <- split(metadataBack$cell.id, metadataBack$dataset)
      
      ds_list <- lapply(names(cells_by_ds), function(ds) {
        cells <- intersect(cells_by_ds[[ds]], colnames(dataBack))
        if (length(cells) < 50L) {
          message(sprintf("[Ref build] Dataset '%s' has <50 cells; skipped.", ds))
          return(NULL)
        }
        list(ds = ds, mat = dataBack[, cells, drop = FALSE], cells = cells)
      })
      ds_list <- Filter(Negate(is.null), ds_list)
      
      if (length(ds_list) == 0L) {
        stop("No dataset has >=50 cells after subsetting; cannot build background.")
      }
      
      # Union of features observed in ANY dataset, padded with NA where absent.
      union_feats <- Reduce(union, lapply(ds_list, function(z) rownames(z$mat)))
      dataBack_union <- do.call(cbind, lapply(ds_list, function(z) {
        mat  <- as.matrix(z$mat)
        miss <- setdiff(union_feats, rownames(mat))
        if (length(miss) > 0L) {
          na_mat <- matrix(NA_real_, nrow = length(miss), ncol = ncol(mat),
                           dimnames = list(miss, colnames(mat)))
          mat <- rbind(mat, na_mat)
        }
        mat[union_feats, , drop = FALSE]
      }))
      
      message(sprintf(
        "[Ref build] Union matrix before global filter: %d markers x %d cells (from %d datasets).",
        nrow(dataBack_union), ncol(dataBack_union), length(ds_list)))
      
      # Single global QC pass (Mode C). Audit CSV lands in the Ref dir so it
      # can be inspected and tied to this exact run.
      audit_csv <- file.path(ref_dir, paste0(file.name, "_excluded_markers.csv"))
      qc_global <- suppressWarnings(cleanData_ADT(
        data_counts           = dataBack_union,
        nmads                 = qc_nmads,
        min.cells.feature     = min.cells.feature,
        min.coverage.fraction = min.coverage.fraction,
        audit_csv_path        = audit_csv,
        qc_outlier_cells      = qc_outlier_cells_resolved
      ))
      dataBack_qc <- qc_global$Counts
      
      if (nrow(dataBack_qc) == 0L) {
        stop(sprintf(
          "All markers were filtered out by the global QC (min.cells.feature=%d, min.coverage.fraction=%.2f). %s",
          as.integer(min.cells.feature), min.coverage.fraction,
          "Lower the thresholds or check the input data."))
      }
      
      message(sprintf(
        "[Ref build] After global filter: %d markers x %d cells.",
        nrow(dataBack_qc), ncol(dataBack_qc)))
      
      kept_all     <- colnames(dataBack_qc)
      metadataBack <- metadataBack[metadataBack$cell.id %in% kept_all, , drop = FALSE]
      metadataBack <- metadataBack[match(colnames(dataBack_qc), metadataBack$cell.id), , drop = FALSE]
      
      dataBack <- dataBack_qc
      
      # Persist the sanitised counts + metadata BEFORE normalisation so that the
      # Ref self-query (precBack=NA branch) can reuse the exact same matrix and
      # cell set on subsequent reruns, even when the rest of the Ref is cached.
      saveRDS(dataBack,     bg_counts_path)
      saveRDS(metadataBack, bg_meta_path)
      
      # BACKGROUND NORMALIZATION
      t_ref <- proc.time()
      norm.data <- suppressWarnings(ADTprot_norm(
        data_counts    = dataBack,
        metadata       = metadataBack,
        data.type      = "Ref",
        file.id        = file.name,
        bg.prot.pars   = NULL,
        saveData       = TRUE,
        return         = TRUE,
        outdir         = ref_dir,
        sample_col     = "sample",
        batch_col      = "batch",
        verbose        = TRUE,
        is_transformed = is_transformed,
        sqs_threshold  = sqs_threshold,
        apply_sqs      = apply_sqs,
        ncores         = ncores,
        ...))
      message(sprintf("[TIMING] ADTprot_norm (Ref) elapsed: %.1fs",
                      (proc.time() - t_ref)[["elapsed"]]))
      bg.prot.pars <- readRDS(bg_pars_path)
      
    } else {
      norm.data    <- readRDS(bg_norm_path)
      bg.prot.pars <- readRDS(bg_pars_path)
      # When the Ref is cached, load the persisted sanitised counts/metadata too.
      # If they don't exist (older Ref built before this fix), fall back to
      # rebuilding the sanitised data using the SAME Mode C logic as the main
      # build path (single global filter), so self-query stays consistent.
      if (file.exists(bg_counts_path) && file.exists(bg_meta_path)) {
        dataBack     <- readRDS(bg_counts_path)
        metadataBack <- readRDS(bg_meta_path)
      } else {
        warning(sprintf(
          "Cached Ref at %s has no persisted BackCounts/BackMetadata (older run). %s",
          ref_dir,
          "Self-query will rebuild the sanitised data using Mode C global filter."
        ))
        cells_by_ds <- split(metadataBack$cell.id, metadataBack$dataset)
        ds_list <- lapply(names(cells_by_ds), function(ds) {
          cells <- intersect(cells_by_ds[[ds]], colnames(dataBack))
          if (length(cells) < 50L) return(NULL)
          list(ds = ds, mat = dataBack[, cells, drop = FALSE], cells = cells)
        })
        ds_list <- Filter(Negate(is.null), ds_list)
        union_feats <- Reduce(union, lapply(ds_list, function(z) rownames(z$mat)))
        dataBack_union <- do.call(cbind, lapply(ds_list, function(z) {
          mat  <- as.matrix(z$mat)
          miss <- setdiff(union_feats, rownames(mat))
          if (length(miss) > 0L) {
            na_mat <- matrix(NA_real_, nrow = length(miss), ncol = ncol(mat),
                             dimnames = list(miss, colnames(mat)))
            mat <- rbind(mat, na_mat)
          }
          mat[union_feats, , drop = FALSE]
        }))
        # Single global filter (same thresholds as the main build).
        qc_global <- suppressWarnings(cleanData_ADT(
          data_counts           = dataBack_union,
          nmads                 = qc_nmads,
          min.cells.feature     = min.cells.feature,
          min.coverage.fraction = min.coverage.fraction,
          audit_csv_path        = file.path(ref_dir, paste0(file.name, "_excluded_markers.csv")),
          qc_outlier_cells      = qc_outlier_cells_resolved
        ))
        dataBack <- qc_global$Counts
        if (nrow(dataBack) == 0L) {
          stop("Fallback rebuild produced an empty matrix after global filter.")
        }
        kept_all     <- colnames(dataBack)
        metadataBack <- metadataBack[metadataBack$cell.id %in% kept_all, , drop = FALSE]
        metadataBack <- metadataBack[match(colnames(dataBack), metadataBack$cell.id), , drop = FALSE]
        # Persist for next time so we don't pay this cost again.
        saveRDS(dataBack,     bg_counts_path)
        saveRDS(metadataBack, bg_meta_path)
      }
    }
    
    if (!file.exists(th_path)) {
      t_bs <- proc.time()
      thrs.data <- GetBackground(norm.data, ref_dir, file.name, metadataBack,
                                 sample_n           = num.cells.sample,
                                 parallelize        = TRUE,
                                 ncores             = ncores,
                                 num.boot.sample    = num.boot.sample,
                                 fixseed            = fixseed,
                                 bootstrap_stratify = bootstrap_stratify)
      message(sprintf("[TIMING] GetBackground/RunBootstrap elapsed: %.1fs",
                      (proc.time() - t_bs)[["elapsed"]]))
    } else {
      thrs.data <- readRDS(th_path)
    }
    
  } else {
    th_file   <- list.files(path = precBack, pattern = "_th_dist\\.rds$", full.names = TRUE)
    pars_file <- list.files(path = precBack, pattern = "_bg_prot_pars\\.rds$", full.names = TRUE)
    if (length(th_file) != 1 || length(pars_file) != 1) stop("Files not found in precBack")
    thrs.data    <- readRDS(th_file[1])
    bg.prot.pars <- readRDS(pars_file[1])
  }
  
  #############################################################################
  ################### PROTEIN QUANTIFICATION (QUERY) ##########################
  #############################################################################
  
  # Branch by mode:
  #   * Ref + self-query (precBack=NA): the "query" is the SAME dataset that
  #     built the Ref. To keep the discretisation methodologically coherent,
  #     it must use the *sanitised* dataBack/metadataBack (per-dataset QC,
  #     min.cells.feature filter, NA-padded union of features) -- not the raw
  #     `data`/`metadata`. Cell-type filter must also match the Ref (>=50).
  #   * External query (precBack supplied): use raw `data`/`metadata` as-is and
  #     apply the looser >=10 cells per cell-type filter, since the query is
  #     independent of the Ref's QC choices.
  if (!use_precBack) {
    query.cell.type.freq    <- as.data.frame(table(metadataBack$cell.type))
    query.cell.type.to.keep <- as.character(query.cell.type.freq[query.cell.type.freq$Freq >= 50, "Var1"])
    metadataQuery <- metadataBack[metadataBack$cell.type %in% query.cell.type.to.keep, , drop = FALSE]
    dataQuery     <- dataBack[, as.character(metadataQuery$cell.id), drop = FALSE]
    message(sprintf(
      "[Ref+self-query] Using sanitised background as query: %d cells x %d features, %d cell types (>=50 cells).",
      ncol(dataQuery), nrow(dataQuery), length(query.cell.type.to.keep)))
  } else {
    query.cell.type.freq    <- as.data.frame(table(metadata$cell.type))
    query.cell.type.to.keep <- as.character(query.cell.type.freq[query.cell.type.freq$Freq >= 10, "Var1"])
    metadataQuery <- metadata[metadata$cell.type %in% query.cell.type.to.keep, , drop = FALSE]
    dataQuery     <- data[, as.character(metadataQuery$cell.id), drop = FALSE]
    message(sprintf(
      "[External query] %d cells x %d features, %d cell types (>=10 cells).",
      ncol(dataQuery), nrow(dataQuery), length(query.cell.type.to.keep)))
  }
  
  if (discretize) {
    set.seed(fixseed)
    base_outdir <- if (!use_precBack) ref_dir else dir
    query_all_dir <- file.path(base_outdir, "Query_ALL")
    if (!dir.exists(query_all_dir)) dir.create(query_all_dir, recursive = TRUE)
    
    full_query_norm_path <- file.path(query_all_dir, "ALL_ScaledData.rds")
    full_query_meta_path <- file.path(query_all_dir, "ALL_metadata_afterQC.rds")
    
    if (file.exists(full_query_norm_path) && file.exists(full_query_meta_path)) {
      query_all_norm <- readRDS(full_query_norm_path)
      metadataQuery2 <- readRDS(full_query_meta_path)
    } else {
      # QUERY PREPARATION
      t_q <- proc.time()
      query_all_norm <- suppressWarnings(PrepQueryData_ADT(
        query_dat             = dataQuery,
        query_metadata        = metadataQuery,
        Ref_data              = if (!use_precBack) ref_dir else precBack,
        dir                   = query_all_dir,
        file.id               = "ALL",
        bg.prot.pars          = bg.prot.pars,
        normalise             = TRUE,
        scale_data            = TRUE,
        saveData              = TRUE,
        return.data           = TRUE,
        is_transformed        = is_transformed,
        sqs_threshold         = sqs_threshold,
        apply_sqs             = apply_sqs,
        min.cells.feature     = min.cells.feature,
        min.coverage.fraction = min.coverage.fraction,
        qc_outlier_cells      = qc_outlier_cells_resolved,
        qc_nmads              = qc_nmads,
        ncores                = ncores
      ))
      message(sprintf("[TIMING] PrepQueryData_ADT (Query ALL) elapsed: %.1fs",
                      (proc.time() - t_q)[["elapsed"]]))
      metadataQuery2 <- metadataQuery[metadataQuery$cell.id %in% colnames(query_all_norm), , drop = FALSE]
      metadataQuery2 <- metadataQuery2[colnames(query_all_norm), , drop = FALSE]
      saveRDS(metadataQuery2, full_query_meta_path)
    }
    
    lapply(query.cell.type.to.keep, function(pop) {
      pop_outdir <- base_outdir
      disc_path <- file.path(pop_outdir, paste0("Query_", pop), paste0(pop, "_DiscretisedData.rds"))
      if (file.exists(disc_path)) return(NULL)
      
      pop.meta <- metadataQuery2[metadataQuery2$cell.type == pop, , drop = FALSE]
      if (nrow(pop.meta) == 0) return(NULL)
      pop.norm <- query_all_norm[, pop.meta$cell.id, drop = FALSE]
      
      pop_dir <- file.path(pop_outdir, paste0("Query_", pop))
      if (!dir.exists(pop_dir)) dir.create(pop_dir, recursive = TRUE)
      saveRDS(pop.norm, file = file.path(pop_dir, paste0(pop, "_NormData.rds")))
      
      Discretise_data(pop.norm, if (!use_precBack) ref_dir else precBack, thrs.data, outdir = pop_outdir, file.id = pop)
    })
  }
  
  if (sig.frames) {
    base_outdir <- if (!use_precBack) ref_dir else dir
    sig_dir <- file.path(base_outdir, "sig_frames")
    if (!dir.exists(sig_dir)) dir.create(sig_dir, recursive = TRUE)
    sig_audit_dir <- file.path(sig_dir, "audit_excluded_genes")
    if (!dir.exists(sig_audit_dir)) dir.create(sig_audit_dir, recursive = TRUE)
    
    lapply(query.cell.type.to.keep, function(pop) {
      dfile <- file.path(base_outdir, paste0("Query_", pop), paste0(pop, "_DiscretisedData.rds"))
      if (!file.exists(dfile)) return(NULL)
      
      discrete.data <- readRDS(dfile)
      audit_path <- file.path(sig_audit_dir, paste0(pop, "_excluded_from_baseline.csv"))
      sig.res <- Create.sig.frame(
        discrete.data$discretised_data,
        min.coverage.fraction.sigframe = min.coverage.fraction.sigframe,
        audit_csv_path                 = audit_path
      )
      # Use the per-population sanitised probs for the binomial calls.
      sig.tbl <- discrete.expr.levels(sig.res[[1]], sig.res[[2]],
                                      prob_low  = sig.res[[3]],
                                      prob_high = sig.res[[4]])
      
      write.table(sig.tbl, file.path(sig_dir, paste0(pop, "_sig_tbl.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    })
  }
  return(NULL)
}
