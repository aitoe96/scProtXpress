#' Unpack a SingleCellExperiment for the CyTOF path
#'
#' Extracts a dense `proteins x cells` matrix from `assay(sce, assay_name)`
#' and a base `data.frame` from `colData(sce)`. Also fills `metadata$batch`
#' from `batch_default_col` when present so the user can call
#' `run_scProteinXpress()` without extra plumbing.
#'
#' @param sce A `SingleCellExperiment`.
#' @param assay_name character(1). Defaults to "exprs" (arcsinh-transformed).
#' @param batch_default_col Column in `colData(sce)` used as batch when
#'   metadata has no `batch` (defaults to "dataset_id").
#' @return A list with `data` (matrix), `metadata` (data.frame), and
#'   `is_transformed` (TRUE if `assay_name == "exprs"`, else NULL for
#'   autodetect).
#' @export
unpack_sce_cytof <- function(sce, assay_name = "exprs",
                             batch_default_col = "dataset_id") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("unpack_sce_cytof requires SingleCellExperiment + SummarizedExperiment. ",
         "Install them or pass `data` + `metadata` directly.")
  }
  data <- as.matrix(SummarizedExperiment::assay(sce, assay_name))
  meta <- as.data.frame(SingleCellExperiment::colData(sce))
  if (!"batch" %in% colnames(meta) && batch_default_col %in% colnames(meta)) {
    meta$batch <- meta[[batch_default_col]]
  }
  list(
    data           = data,
    metadata       = meta,
    is_transformed = if (identical(assay_name, "exprs")) TRUE else NULL
  )
}

#' Unpack the CITE-seq query list `{$counts, $cell_x_feature}`
#'
#' Pulls `$counts` (dgCMatrix or coercible) and `$cell_x_feature`
#' (data.frame) out of the conventional CITE-seq query object and returns
#' them in the canonical `list(data, metadata, is_transformed)` shape
#' expected by `run_scProteinXpress()`.
#'
#' @param query_obj A list with `counts` and `cell_x_feature`.
#' @return list(data, metadata, is_transformed = NULL).
#' @export
unpack_list_citeseq <- function(query_obj) {
  if (!is.list(query_obj) || is.null(query_obj$counts) ||
      is.null(query_obj$cell_x_feature)) {
    stop("unpack_list_citeseq: expected a list with `counts` and `cell_x_feature`.")
  }
  list(
    data           = query_obj$counts,
    metadata       = as.data.frame(query_obj$cell_x_feature,
                                   stringsAsFactors = FALSE),
    is_transformed = NULL
  )
}
