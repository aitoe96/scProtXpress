## Landmark detection on cell-type balanced copies and the alignment fixes of
## 0.3.0: per-iteration bootstrap seeds (bug C), per-sample warp fallback
## (bug B), missing valley imputed, per-marker seeds, minimum cells per sample,
## and queries that follow the detection settings of their background; the
## positive peak of the background, single-peak query samples and the debug
## files.

## Toy background: 3 cell types with different marker levels, 4 samples; the
## last sample has only 12 cells (below min_cells_align = 25).
.toy_balance_data <- function(seed = 7L) {
  set.seed(seed)
  n_by_sample <- c(s1 = 150L, s2 = 120L, s3 = 90L, s4 = 12L)
  sample <- rep(names(n_by_sample), n_by_sample)
  ct <- unlist(lapply(n_by_sample, function(n) sample(c("B", "T", "mono"), n, replace = TRUE,
                                                     prob = c(0.6, 0.3, 0.1))), use.names = FALSE)
  level <- c(B = 0.5, T = 3, mono = 1.8)
  counts <- t(vapply(seq_len(4L), function(k) {
    stats::rnorm(length(ct), mean = level[ct] * (1 + 0.2 * (k - 1)), sd = 0.35)
  }, numeric(length(ct))))
  counts[counts < 0] <- 0
  cell_id <- sprintf("c%04d", seq_along(ct))
  dimnames(counts) <- list(paste0("M", 1:4), cell_id)
  meta <- data.frame(cell.id = cell_id, cell.type = ct, sample = sample, batch = sample,
                     tissue = "PBMC", dataset = "DS1", row.names = cell_id, stringsAsFactors = FALSE)
  list(counts = counts, meta = meta)
}

.toy_ref <- function(toy, outdir, ...) {
  suppressMessages(ADTprot_norm(toy$counts, toy$meta, data.type = "Ref", file.id = "toy",
                                outdir = outdir, is_transformed = TRUE, do_subsample = FALSE,
                                verbose = FALSE, ...))
}

test_that(".balance_index gives the sizes of the balanced copy", {
  bi <- scProtXpress:::.balance_index
  ct <- c(rep("a", 150), rep("b", 5), rep("c", 1))
  set.seed(1)
  i <- bi(ct, 300, 1L, FALSE)                         # K = 3: 100 cells per type
  expect_equal(as.vector(table(ct[i])), c(100L, 100L, 100L))
  expect_true(all((150L + seq_len(5)) %in% i))        # a smaller type keeps all its cells
  expect_false(any(duplicated(i[ct[i] == "a"])))      # a larger one is down-sampled without replacement
  i <- bi(ct, 300, 2L, FALSE)                         # "c" (1 cell) left out: K = 2
  expect_equal(as.vector(table(factor(ct[i], levels = c("a", "b", "c")))), c(150L, 150L, 0L))
  expect_length(bi(ct, "sample", 1L, FALSE), 3L * round(156 / 3))
  expect_equal(as.vector(table(ct[bi(ct, "maxK", 1L, FALSE)])), c(150L, 150L, 150L))
  expect_identical(bi(rep("a", 10), 300, 1L, FALSE), seq_len(10))   # K < 2: sample as it is
  expect_identical(bi(ct, 300, 10L, FALSE), seq_along(ct))          # one eligible type
  ## types in C-locale order ("A" < "a" < "b") whatever the locale of the session
  lab <- c("b", "a", "A", "b", "A")
  expect_identical(unique(lab[bi(lab, 6, 1L, FALSE)]), c("A", "a", "b"))
})

test_that("RunBootstrap draws distinct replicates and does not depend on the number of cores", {
  skip_on_cran()
  set.seed(11)
  m <- matrix(abs(stats::rnorm(3 * 400, 2, 1)), nrow = 3,
              dimnames = list(paste0("M", 1:3), paste0("c", 1:400)))
  meta <- data.frame(cell.id = colnames(m), cell.type = rep(c("x", "y"), 200), row.names = colnames(m))
  run <- function(nc) suppressMessages(RunBootstrap(m, tempdir(), "bt", "cell.type", meta, sample_n = 50,
                                                    parallelize = nc > 1, ncores = nc, num.boot.sample = 20,
                                                    fixseed = 1234, bootstrap_stratify = "cell.type"))
  b1 <- run(1L)
  expect_equal(nrow(unique(b1[["M1"]])), 20L)
  expect_identical(run(2L), b1)
})

test_that("bootstrap and subsampling do not depend on the collation, which is restored", {
  skip_on_cran()
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
  skip_if(identical(suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8")), ""),
          "en_US.UTF-8 collation not available")
  en <- Sys.getlocale("LC_COLLATE")
  ## "B" sorts before "a" in C and after it in en_US; same for the tissues
  set.seed(11)
  m <- matrix(abs(stats::rnorm(2 * 400, 2, 1)), nrow = 2,
              dimnames = list(paste0("M", 1:2), paste0("c", 1:400)))
  meta <- data.frame(cell.id = colnames(m), cell.type = rep(c("a", "B"), 200),
                     tissue = rep(c("PBMC", "bone marrow"), each = 200), row.names = colnames(m))
  boot <- function() suppressMessages(RunBootstrap(m, tempdir(), "bt", c("cell.type", "tissue"), meta,
                                                   sample_n = 50, parallelize = FALSE, ncores = 1L,
                                                   num.boot.sample = 10, fixseed = 1234,
                                                   bootstrap_stratify = "cell.type_tissue"))
  toy <- .toy_balance_data()
  toy$meta$cell.type <- c(B = "B", T = "t", mono = "mono")[toy$meta$cell.type]
  toy$meta$tissue <- rep(c("PBMC", "bone marrow"), length.out = nrow(toy$meta))
  toy$meta$dataset <- rep(c("DS1", "ds2"), length.out = nrow(toy$meta))
  subs <- function() {
    d <- tempfile("collate_"); dir.create(d)
    colnames(suppressMessages(ADTprot_norm(toy$counts, toy$meta, data.type = "Ref", file.id = "toy",
                                           outdir = d, is_transformed = TRUE, do_subsample = TRUE,
                                           target_ct = 20L, verbose = FALSE)))
  }
  b_en <- boot(); s_en <- subs()
  expect_identical(Sys.getlocale("LC_COLLATE"), en)
  Sys.setlocale("LC_COLLATE", "C")
  expect_identical(boot(), b_en)
  expect_identical(subs(), s_en)
  expect_identical(Sys.getlocale("LC_COLLATE"), "C")
})

test_that(".peak_alignment_by_sample keeps the samples that do not fail", {
  set.seed(5)
  x <- c(stats::rnorm(300, 1, 0.3), stats::rnorm(300, 1.2, 0.3), stats::rnorm(300, 0.9, 0.3))
  x <- x + rep(c(0, 3), length.out = 900) * (stats::runif(900) < 0.3)
  fea <- data.frame(sample = factor(rep(c("a", "b", "c"), each = 300)))
  lm <- rbind(a = c(1, 2.5, 4), b = c(1.2, 2.7, 4.2), c = c(0.9, 2.4, 3.9))
  ok <- scProtXpress:::.peak_alignment_by_sample(x, fea, lm, NULL)
  expect_identical(ok$warp_failed, character(0))
  expect_equal(ok$cell_x_adt_norm, scProtXpress:::.peak_alignment_scprot(x, fea, lm, NULL)$cell_x_adt_norm)

  ## fda::landmarkreg failing for sample "b" (the failure is numerically
  ## chaotic, so it is simulated): the whole call and "b" alone fail
  real <- scProtXpress:::.peak_alignment_scprot
  local_mocked_bindings(.peak_alignment_scprot = function(cell_x_adt, cell_x_feature, ...) {
    if ("b" %in% levels(cell_x_feature$sample)) stop("leading minor of order 2 is not positive")
    real(cell_x_adt, cell_x_feature, ...)
  })
  res <- scProtXpress:::.peak_alignment_by_sample(x, fea, lm, NULL)
  expect_identical(res$warp_failed, "b")
  expect_true(all(is.na(res$cell_x_adt_norm[fea$sample == "b"])))
  ## the other samples get the values of the whole call that did not fail
  expect_equal(res$cell_x_adt_norm[fea$sample != "b"], ok$cell_x_adt_norm[fea$sample != "b"])
  expect_equal(res$warp_parameters, ok$warp_parameters)
  expect_error(scProtXpress:::.peak_alignment_by_sample(x[fea$sample == "b"],
                                                        data.frame(sample = factor(rep("b", 300))),
                                                        lm["b", , drop = FALSE], NULL),
               "every sample")
})

test_that("a single peak without valley or shoulder gets a missing valley, not an error", {
  set.seed(9)
  x <- c(stats::rnorm(200, 1, 0.3), stats::rnorm(200, 1.1, 0.3))
  m <- matrix(x, ncol = 1, dimnames = list(NULL, "tmpName1"))
  fea <- data.frame(sample = factor(rep(c("a", "b"), each = 200)), batch = rep(c("a", "b"), each = 200))
  ## sample "b": its only peak is to the right of every value -> no valley, no shoulder
  pk <- matrix(c(1, 100), ncol = 1, dimnames = list(c("a", "b"), NULL))
  args <- list(m, fea, "tmpName1", pk, TRUE, NULL, FALSE, adjust = 3, min_fc = 20,
               shoulder_valley_slope = -0.5, neg_candidate_thres = asinh(1/5 + 1),
               arcsinh_transform_flag = FALSE)
  expect_error(suppressWarnings(do.call(ADTnorm::get_valley_location, args)), "No valley is detected")
  res <- suppressWarnings(do.call(scProtXpress:::.get_valley_location_patched, args))
  expect_identical(res$no_valley_samples, "b")
  expect_true(is.na(res$valley_location_list["b", 1]))
  expect_false(is.na(res$valley_location_list["a", 1]))
})

test_that("Ref with balancing: settings saved, small samples NA and logged, same result with 1 or 2 cores", {
  skip_on_cran()
  toy <- .toy_balance_data()
  out1 <- tempfile("bal1-"); dir.create(out1)
  out2 <- tempfile("bal2-"); dir.create(out2)
  n1 <- .toy_ref(toy, out1, ncores = 1L)
  n2 <- .toy_ref(toy, out2, ncores = 2L)
  expect_identical(n1, n2)

  gs <- readRDS(file.path(out1, "toy_bg_prot_pars.rds"))$global_settings
  expect_true(gs$balance_celltypes)
  expect_equal(gs$balance_n, 2000L)
  expect_equal(gs$min_cells_align, 25L)
  expect_setequal(gs$bimodal_marker, paste0("M", 1:4))    # "auto" with balancing -> every marker
  expect_identical(gs$alignment_version, "0.3.0")

  s4 <- toy$meta$sample == "s4"
  expect_true(all(is.na(n1[, s4])))
  expect_false(anyNA(n1[, !s4]))
  log <- utils::read.csv(file.path(out1, "toy_alignment_log.csv"))
  expect_true(all(c("stage", "file.id", "marker", "sample", "reason") %in% names(log)))
  expect_setequal(log$marker[log$reason == "too_small"], paste0("M", 1:4))
  expect_true(all(log$sample[log$reason == "too_small"] == "s4"))

  ## balancing changes the detection input, hence the result
  out3 <- tempfile("bal3-"); dir.create(out3)
  n3 <- .toy_ref(toy, out3, balance_celltypes = FALSE)
  expect_false(identical(n1, n3))
  expect_null(readRDS(file.path(out3, "toy_bg_prot_pars.rds"))$global_settings$bimodal_marker)

  ## same seed -> same result; balancing without a cell type fails clearly
  out4 <- tempfile("bal4-"); dir.create(out4)
  expect_identical(.toy_ref(toy, out4, ncores = 1L), n1)
  meta_noct <- toy$meta; meta_noct$cell.type <- NULL
  expect_error(suppressMessages(ADTprot_norm(toy$counts, meta_noct, data.type = "Ref", file.id = "x",
                                             outdir = tempdir(), is_transformed = TRUE,
                                             do_subsample = FALSE, verbose = FALSE)),
               "balance_celltypes = FALSE")
})

test_that("queries follow the bimodal markers of their background; balancing only with balance_celltypes_query", {
  skip_on_cran()
  toy <- .toy_balance_data()
  out <- tempfile("balq-"); dir.create(out)
  .toy_ref(toy, out, bimodal_marker = "M2")
  bg <- readRDS(file.path(out, "toy_bg_prot_pars.rds"))
  expect_identical(bg$global_settings$bimodal_marker, "M2")

  seen <- new.env(); seen$bimodal <- list(); seen$balanced <- 0L
  real_peak <- scProtXpress:::.get_peak_mode_patched
  real_bal  <- scProtXpress:::.balanced_detection_input
  local_mocked_bindings(
    .get_peak_mode_patched = function(cell_x_adt, cell_x_feature, adt_marker_select, adt_marker_index,
                                      bwFac_smallest, bimodal_marker_index, ...) {
      seen$bimodal[[adt_marker_select]] <- adt_marker_index %in% bimodal_marker_index
      real_peak(cell_x_adt, cell_x_feature, adt_marker_select, adt_marker_index, bwFac_smallest,
                bimodal_marker_index, ...)
    },
    .balanced_detection_input = function(...) { seen$balanced <- seen$balanced + 1L; real_bal(...) })
  qdir <- tempfile("q-"); dir.create(qdir)
  q <- suppressMessages(ADTprot_norm(toy$counts, toy$meta, data.type = "query", file.id = "q",
                                     bg.prot.pars = bg, outdir = qdir, is_transformed = TRUE,
                                     verbose = FALSE))
  expect_identical(unlist(seen$bimodal)[order(names(seen$bimodal))],
                   c(tmpName1 = FALSE, tmpName2 = TRUE, tmpName3 = FALSE, tmpName4 = FALSE))
  expect_equal(seen$balanced, 0L)
  expect_true(all(is.na(q[, toy$meta$sample == "s4"])))
  expect_true(file.exists(file.path(qdir, "q_alignment_log.csv")))

  suppressMessages(ADTprot_norm(toy$counts, toy$meta, data.type = "query", file.id = "q_bal",
                                bg.prot.pars = bg, outdir = qdir, is_transformed = TRUE,
                                verbose = FALSE, balance_celltypes_query = TRUE))
  expect_equal(seen$balanced, 4L)

  ## background of scProtXpress <= 0.2.0: no balancing, with a message
  old <- bg; old$global_settings$alignment_version <- NULL
  seen$balanced <- 0L
  expect_message(ADTprot_norm(toy$counts, toy$meta, data.type = "query", file.id = "q_old",
                              bg.prot.pars = old, outdir = qdir, is_transformed = TRUE, verbose = FALSE,
                              balance_celltypes_query = TRUE),
                 "scProtXpress <= 0.2.0")
  expect_equal(seen$balanced, 0L)
})

test_that("Ref positive peak: median of the per-sample maxima if no sample detects it, sample maximum if at or below the valley", {
  smp <- factor(rep(c("a", "b", "c"), each = 2))
  x <- c(1, 3, 1, 5, 1, 6)                                   # maxima 3, 5, 6
  lm <- cbind(c(a = 1, b = 1.2, c = 0.9), c(2, 2.2, 5.5))    # negative peak, valley
  r <- scProtXpress:::.impute_positive_peak(lm, x, smp, c("a", "b", "c"))
  expect_equal(unname(r$landmark_matrix[, 3]), c(5, 5, 6))   # c: median 5 <= valley 5.5
  expect_identical(unname(r$method), c("median_of_max", "median_of_max", "sample_max"))
  r <- scProtXpress:::.impute_positive_peak(cbind(lm, c(4, 1.5, 3)), x, smp, c("a", "b", "c"))
  expect_equal(unname(r$landmark_matrix[, 3]), c(4, 5, 6))
  expect_identical(unname(r$method), c(NA, "sample_max", "sample_max"))
})

test_that(".align_by_landmark_set warps every set of landmarks on its own", {
  set.seed(5)
  x <- c(stats::rnorm(300, 1, 0.3), stats::rnorm(300, 1.2, 0.3), stats::rnorm(300, 0.9, 0.3))
  x <- x + rep(c(0, 3), length.out = 900) * (stats::runif(900) < 0.3)
  smp <- factor(rep(c("a", "b", "c"), each = 300))
  lm <- rbind(a = c(1, 2.5, 4), b = c(1.2, 2.7, 4.2), c = c(0.9, 2.4, 3.9))
  tl <- c(1, 2.5, 4)
  wp <- scProtXpress:::.peak_alignment_scprot(x, data.frame(sample = smp), lm, tl)$warp_parameters
  all3 <- scProtXpress:::.peak_alignment_by_sample(x, data.frame(sample = smp), lm, tl, "fixed", wp)
  used <- matrix(TRUE, 3, 3, dimnames = list(rownames(lm), NULL))
  expect_equal(scProtXpress:::.align_by_landmark_set(x, smp, lm, used, tl, wp)$cell_x_adt_norm,
               all3$cell_x_adt_norm)
  used["b", 3] <- FALSE
  r <- scProtXpress:::.align_by_landmark_set(x, smp, lm, used, tl, wp)
  b <- smp == "b"
  expect_equal(r$cell_x_adt_norm[!b], all3$cell_x_adt_norm[!b])
  expect_equal(r$cell_x_adt_norm[b],
               scProtXpress:::.peak_alignment_by_sample(x[b], data.frame(sample = factor(rep("b", sum(b)))),
                                                        lm["b", 1:2, drop = FALSE], tl[1:2], "fixed",
                                                        wp)$cell_x_adt_norm)
})

test_that("query single peaks: nearest background target or always the negative one, no positive peak imputed", {
  skip_on_cran()
  toy <- .toy_balance_data()
  out <- tempfile("sp-"); dir.create(out)
  .toy_ref(toy, out, bimodal_marker = NULL)
  bg <- readRDS(file.path(out, "toy_bg_prot_pars.rds"))
  tl <- lapply(bg$markers, `[[`, "fixed_target_landmark")
  expect_true(all(lengths(tl) == 3L))

  ## q1: two peaks; q2: one peak at the negative target; q3: one peak at the positive target
  set.seed(3)
  n <- 300L
  counts <- t(vapply(names(tl), function(m) {
    t3 <- tl[[m]]
    pmax(0, c(ifelse(stats::runif(n) < 0.6, stats::rnorm(n, t3[1], 0.2), stats::rnorm(n, t3[3], 0.2)),
              stats::rnorm(n, t3[1], 0.2), stats::rnorm(n, t3[3], 0.2)))
  }, numeric(3L * n)))
  cell_id <- sprintf("q%04d", seq_len(3L * n))
  colnames(counts) <- cell_id
  meta <- data.frame(cell.id = cell_id, sample = rep(c("q1", "q2", "q3"), each = n), row.names = cell_id)
  meta$batch <- meta$sample
  run <- function(single_peak) {
    qdir <- tempfile("spq-"); dir.create(qdir)
    q <- suppressMessages(ADTprot_norm(counts, meta, data.type = "query", file.id = "q", bg.prot.pars = bg,
                                       outdir = qdir, is_transformed = TRUE, verbose = FALSE,
                                       single_peak = single_peak, debug = TRUE))
    list(q = q, dbg = utils::read.csv(file.path(qdir, "q_debug_landmarks.csv"), stringsAsFactors = FALSE))
  }
  near <- run("nearest")
  neg  <- run("negative")
  for (m in names(tl)) {
    d <- near$dbg[near$dbg$marker == m, ]; rownames(d) <- d$sample
    expect_identical(d[c("q2", "q3"), "single_peak"], c("negative", "positive"))
    expect_identical(d[c("q2", "q3"), "landmarks_used"], c("neg,valley", "valley,pos"))
    d <- neg$dbg[neg$dbg$marker == m, ]; rownames(d) <- d$sample
    expect_identical(d["q3", "landmarks_used"], "neg,valley")
  }
  ## the aligned single peak of q3 ends nearer to the positive target (nearest)
  ## or to the negative one (negative); markers whose positive target lies in
  ## the warp range of the background, to which the landmarks are clipped
  q3 <- meta$sample == "q3"
  ub <- vapply(bg$markers, function(p) p$warp_parameters$upper_bound, numeric(1))
  in_range <- names(tl)[vapply(tl, `[`, numeric(1), 3) <= ub[names(tl)]]
  expect_gt(length(in_range), 0L)
  for (m in in_range) {
    dist <- function(q, j) abs(stats::median(q[m, q3]) - tl[[m]][j])
    expect_lt(dist(near$q, 3), dist(near$q, 1))
    expect_lt(dist(neg$q, 1), dist(neg$q, 3))
  }
})

test_that("debug = TRUE writes the landmarks of every marker and sample, the balanced cells and the settings", {
  skip_on_cran()
  toy <- .toy_balance_data()
  out0 <- tempfile("dbg0-"); dir.create(out0)
  out1 <- tempfile("dbg1-"); dir.create(out1)
  n0 <- .toy_ref(toy, out0)
  n1 <- .toy_ref(toy, out1, debug = TRUE)
  expect_identical(n1, n0)
  tab <- utils::read.csv(file.path(out1, "toy_debug_landmarks.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(tab), 16L)                                # 4 markers x 4 samples
  expect_true(all(tab$status[tab$sample == "s4"] == "too_small"))
  expect_true(all(tab$status[tab$sample != "s4"] == "aligned"))
  expect_true(all(tab$landmarks_used[tab$sample != "s4"] == "neg,valley,pos"))
  r <- readRDS(file.path(out1, "toy_debug", "M1.rds"))
  expect_true(all(c("peaks_raw", "peaks", "valleys", "lm", "targets", "balanced_cells", "marker_pars") %in% names(r)))
  expect_setequal(unique(r$balanced_cells$sample), c("s1", "s2", "s3"))
  expect_true(all(r$balanced_cells$cell_id %in% toy$meta$cell.id))
  expect_gt(sum(r$balanced_cells$n), nrow(r$balanced_cells))  # cells drawn more than once
  expect_true(file.exists(file.path(out1, "toy_debug_settings.csv")))
})

test_that("0.3.1 FDR: p_low adjusted over every gene, p_high only over the genes tested for medium/high", {
  ## gene g1: clearly expressed and intermediate; g2-g5: not expressed; g6: tie
  tbl <- data.frame(gene = paste0("g", 1:6),
                    n.0   = c(10, 95, 96, 97, 98, 50),
                    n.0.5 = c(80,  3,  2,  2,  1, 25),
                    n.1   = c(10,  2,  2,  1,  1, 25))
  r <- scProtXpress:::discrete.expr.levels(tbl, NULL, prob_low = 0.3, prob_high = 0.5, fdr_method = "BH")
  expect_equal(r$q_low, stats::p.adjust(r$p_low, "BH"))
  tested <- r$expr.level %in% c("medium", "high")
  expect_true(all(is.na(r$q_high[!tested])))
  expect_equal(r$q_high[tested], stats::p.adjust(r$p_high[tested], "BH"))
  expect_true(all(c("expr.level_raw", "q_low", "q_high") %in% names(r)))
  r0 <- scProtXpress:::discrete.expr.levels(tbl, NULL, prob_low = 0.3, prob_high = 0.5)
  expect_identical(r0$expr.level, r$expr.level_raw)
})
