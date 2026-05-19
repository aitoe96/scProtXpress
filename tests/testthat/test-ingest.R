test_that("unpack_sce_cytof requires SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")

  toy <- make_toy_cytof()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(exprs = toy$data),
    colData = S4Vectors::DataFrame(toy$metadata)
  )

  unpacked <- unpack_sce_cytof(sce, assay_name = "exprs",
                               batch_default_col = "dataset_id")
  expect_true(is.matrix(unpacked$data))
  expect_equal(dim(unpacked$data), dim(toy$data))
  expect_true("batch" %in% colnames(unpacked$metadata))
  expect_true(isTRUE(unpacked$is_transformed))
})

test_that("unpack_list_citeseq returns canonical shape", {
  qobj <- make_toy_citeseq()
  out  <- unpack_list_citeseq(qobj)

  expect_true(inherits(out$data, "dgCMatrix"))
  expect_true(is.data.frame(out$metadata))
  expect_null(out$is_transformed)
  expect_true("cell.type_harmonised" %in% colnames(out$metadata))
})

test_that("unpack_list_citeseq fails on a malformed input", {
  expect_error(unpack_list_citeseq(list(counts = NULL,
                                        cell_x_feature = data.frame())),
               "counts.*cell_x_feature")
  expect_error(unpack_list_citeseq(list(counts = matrix(0, 2, 2))),
               "counts.*cell_x_feature")
})
