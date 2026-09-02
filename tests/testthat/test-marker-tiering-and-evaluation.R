## Tests for tag_marker_specificity() and evaluate_calls().

test_that("tag_marker_specificity correctly classifies broad vs specific markers", {
  map <- data.frame(
    marker    = c("CD3", "CD3", "CD3", "CD19", "CD14", "CD56", "CD56"),
    cell.type = c("Tcell", "NKT", "Treg", "Bcell", "Monocyte", "NK", "NKT"),
    stringsAsFactors = FALSE
  )
  tiers <- tag_marker_specificity(map)

  expect_equal(tiers$tier[tiers$marker == "CD3"],  "broad")
  expect_equal(tiers$n_celltypes[tiers$marker == "CD3"], 3L)
  expect_equal(tiers$tier[tiers$marker == "CD19"], "specific")
  expect_equal(tiers$n_celltypes[tiers$marker == "CD19"], 1L)
  expect_equal(tiers$tier[tiers$marker == "CD56"], "broad")
  expect_equal(tiers$n_celltypes[tiers$marker == "CD56"], 2L)
})

test_that("tag_marker_specificity deduplicates repeated (marker, cell type) rows before counting", {
  map <- data.frame(
    marker    = c("CD3", "CD3", "CD3"),        # same pair listed twice
    cell.type = c("Tcell", "Tcell", "NKT"),
    stringsAsFactors = FALSE
  )
  tiers <- tag_marker_specificity(map)
  expect_equal(tiers$n_celltypes[tiers$marker == "CD3"], 2L)  # not 3
})

test_that("tag_marker_specificity errors clearly on missing columns", {
  expect_error(
    tag_marker_specificity(data.frame(gene = "CD3", type = "Tcell")),
    "marker.*cell.type"
  )
})

test_that("evaluate_calls computes the correct confusion matrix and metrics by hand", {
  gt <- data.frame(marker = c("CD3", "CD19", "CD56"), cell.type = c("Tcell", "Bcell", "NK"),
                   stringsAsFactors = FALSE)
  calls <- data.frame(
    marker     = c("CD3", "CD19", "CD56", "CD14"),
    cell.type  = c("Tcell", "Bcell", "NK", "Bcell"),
    expr.level = c("high", "low", "medium", "high"),
    tier       = c("broad", "specific", "broad", "specific"),
    stringsAsFactors = FALSE
  )
  res <- evaluate_calls(calls, gt, group_col = "tier")
  overall <- res[res$group == "overall", ]

  ## CD3-Tcell (high) and CD56-NK (medium) are true positives correctly
  ## called; CD19-Bcell (low) is a true positive MISSED (false negative);
  ## CD14 is excluded entirely (not in ground_truth's marker universe) by
  ## the default restrict_to_ground_truth_markers = TRUE.
  expect_equal(overall$TP, 2L)
  expect_equal(overall$FN, 1L)
  expect_equal(overall$FP, 0L)
  expect_equal(overall$precision, 1.0)
  expect_equal(overall$recall, 2/3)
  expect_equal(overall$f1, 0.8)

  broad <- res[res$group == "broad", ]
  expect_equal(broad$TP, 2L); expect_equal(broad$FN, 0L)
  specific <- res[res$group == "specific", ]
  expect_equal(specific$TP, 0L); expect_equal(specific$FN, 1L)
})

test_that("restrict_to_ground_truth_markers = FALSE counts an off-ground-truth marker as a false positive instead of excluding it", {
  gt <- data.frame(marker = "CD3", cell.type = "Tcell", stringsAsFactors = FALSE)
  calls <- data.frame(marker = c("CD3", "CD14"), cell.type = c("Tcell", "Bcell"),
                      expr.level = c("high", "high"), stringsAsFactors = FALSE)

  restricted <- evaluate_calls(calls, gt)
  expect_equal(restricted$n_pairs_evaluated, 1L)   # CD14 excluded

  unrestricted <- evaluate_calls(calls, gt, restrict_to_ground_truth_markers = FALSE)
  expect_equal(unrestricted$n_pairs_evaluated, 2L)
  expect_equal(unrestricted$FP, 1L)                # CD14-Bcell now counted as a false positive
})

test_that("evaluate_calls errors clearly on missing required columns", {
  gt <- data.frame(marker = "CD3", cell.type = "Tcell", stringsAsFactors = FALSE)
  expect_error(
    evaluate_calls(data.frame(marker = "CD3", cell.type = "Tcell"), gt),
    "expr.level"
  )
  expect_error(
    evaluate_calls(data.frame(marker = "CD3", cell.type = "Tcell", expr.level = "high"),
                  data.frame(gene = "CD3")),
    "marker.*cell.type"
  )
})
