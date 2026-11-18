## Thresholds per group of cell types (background_groups / query_groups):
## validation of the tables, outputs per group, groups not built, populations
## without a group, and that the global outputs do not change (also when the
## groups are added to a background built before without them).

test_that("group tables are validated", {
  chk <- scProtXpress:::.check_group_table
  expect_null(chk(NULL, "x"))
  expect_error(chk(data.frame(g = "a", l = "b"), "x"), "columns 'group' and 'label'")
  expect_error(chk(data.frame(group = c("a", NA), label = c("b", "c")), "x"), "empty or missing")
  expect_equal(nrow(chk(data.frame(group = c("a", "a"), label = c("b", "b")), "x")), 1L)
})

.run_toy <- function(toy, outdir, ...) {
  suppressWarnings(suppressMessages(run_scProteinXpress(
    data = toy$data, metadata = toy$metadata, dir = outdir, file.name = "toy", modality = "CyTOF",
    cell_id_col = "cell_id", cell_type_col = "cell_type", dataset_col = "dataset_id",
    sample_col = "sample_id", batch_col = "batch", tissue_default = "PBMC", is_transformed = TRUE,
    num.boot.sample = 5L, num.cells.sample = 20L, min.cells.feature = 1L, min.coverage.fraction = 0,
    ncores = 1L, do_subsample = FALSE, bootstrap_stratify = "cell.type", ...)))
}

test_that("background and query thresholds per group", {
  skip_if_not_installed("ADTnorm")
  skip_if_not_installed("flowCore")
  skip_on_cran()
  toy <- make_toy_cytof(n_cells_per_ct = 120L, n_markers = 6L)
  ## G12 = ct1 + ct2 (240 cells, 3 donors): built; G3 = ct3 (120 cells): not built
  bg_groups <- data.frame(group = c("G12", "G12", "G3"), label = c("ct1", "ct2", "ct3"))

  d_plain <- tempfile("grp_plain_"); d_grp <- tempfile("grp_"); d_late <- tempfile("grp_late_")
  for (d in c(d_plain, d_grp, d_late)) dir.create(d)
  .run_toy(toy, d_plain)
  .run_toy(toy, d_grp, background_groups = bg_groups)
  ref_plain <- file.path(d_plain, "Ref_toy"); ref_grp <- file.path(d_grp, "Ref_toy")

  ## groups: saved, G3 not built, self-query discretised per group
  thg <- readRDS(file.path(ref_grp, "toy_th_dist_groups.rds"))
  expect_identical(names(thg), "G12")
  log <- utils::read.csv(file.path(ref_grp, "toy_groups.csv"))
  expect_identical(log$built, c(TRUE, FALSE))
  expect_identical(log$n_cells, c(240L, 120L))
  expect_true(file.exists(file.path(ref_grp, "groups", "G12", "Query_ct1", "ct1_DiscretisedData.rds")))
  expect_true(file.exists(file.path(ref_grp, "groups", "G12", "sig_frames", "ct2_sig_tbl.txt")))
  expect_false(dir.exists(file.path(ref_grp, "groups", "G12", "Query_ct3")))
  um <- utils::read.csv(file.path(ref_grp, "groups_unmatched.csv"))
  expect_identical(um$population, "ct3")
  expect_identical(um$reason, "groups not built")

  ## the global outputs do not change with the groups
  for (f in c("toy_NormData.rds", "toy_th_dist.rds", "Query_ALL/ALL_ScaledData.rds",
              "Query_ct1/ct1_DiscretisedData.rds")) {
    expect_identical(readRDS(file.path(ref_grp, f)), readRDS(file.path(ref_plain, f)), info = f)
  }
  expect_identical(readLines(file.path(ref_grp, "sig_frames", "ct3_sig_tbl.txt")),
                   readLines(file.path(ref_plain, "sig_frames", "ct3_sig_tbl.txt")))

  ## groups added later to a background built without them: the cached
  ## background is reused and the result equals building it with groups
  .run_toy(toy, d_late)
  th_before <- readRDS(file.path(d_late, "Ref_toy", "toy_th_dist.rds"))
  .run_toy(toy, d_late, background_groups = bg_groups)
  expect_identical(readRDS(file.path(d_late, "Ref_toy", "toy_th_dist.rds")), th_before)
  expect_identical(readRDS(file.path(d_late, "Ref_toy", "toy_th_dist_groups.rds")), thg)
  expect_identical(readRDS(file.path(d_late, "Ref_toy", "groups", "G12", "Query_ct1", "ct1_DiscretisedData.rds")),
                   readRDS(file.path(ref_grp, "groups", "G12", "Query_ct1", "ct1_DiscretisedData.rds")))
  ## other groups on the same background are refused
  expect_error(.run_toy(toy, d_late, background_groups = bg_groups[1:2, ]), "built with other groups")

  ## query: ct1, ct2 -> G12; ct3 has no group
  dq <- tempfile("grp_query_"); dir.create(dq)
  .run_toy(toy, dq, precBack = ref_grp,
           query_groups = data.frame(group = c("G12", "G12"), label = c("ct1", "ct2")))
  expect_true(file.exists(file.path(dq, "groups", "G12", "Query_ct2", "ct2_DiscretisedData.rds")))
  expect_true(file.exists(file.path(dq, "groups", "G12", "sig_frames", "ct1_sig_tbl.txt")))
  expect_true(file.exists(file.path(dq, "Query_ct3", "ct3_DiscretisedData.rds")))
  um_q <- utils::read.csv(file.path(dq, "groups_unmatched.csv"))
  expect_identical(um_q$reason, "no group")
  ## query groups need a background with groups
  dq2 <- tempfile("grp_query2_"); dir.create(dq2)
  expect_error(.run_toy(toy, dq2, precBack = ref_plain, query_groups = data.frame(group = "G12", label = "ct1")),
               "no thresholds per group")
})
