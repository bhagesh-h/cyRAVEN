# The contract fig_marker_umaps_by_group() has to keep.
#
# 1. The POOLED panel is written for every marker, always. Each facet of a split
#    figure holds only a subset of the cells, so without it nothing shows a
#    marker over the whole cohort at full size.
# 2. EVERY category present is split by, not only the one named as
#    --group-column. Naming a statistical group says nothing about which
#    categories are worth seeing, and a cohort usually carries several.
# 3. Numeric columns and identifiers are never faceted. Faceting an age gives one
#    panel per distinct age; faceting patient_id multiplies the folder by the
#    donor count.
#
# All three have failed at some point. An earlier form wrote the pooled panel OR
# the split one, never both. A later one split only by --group-column, so a run
# grouped by timepoint drew a single merged panel for infection focus even though
# the column was in the sheet and on the UMAP overview.
#
# Note a full-pipeline test cannot cover the no-category path: the pipeline
# defaults its group column to "cohort", and the demonstration sheet has one.

marker_cells <- function(n = 300L, cats = TRUE) {
  set.seed(11)
  d <- data.frame(
    umap_1 = rnorm(n), umap_2 = rnorm(n),
    CD3 = rnorm(n), CD4 = rnorm(n),
    sample_id = rep(paste0("s", 1:6), length.out = n),
    patient_id = rep(paste0("p", 1:3), length.out = n),
    age_years = rep(c(61, 72, 55), length.out = n),
    stringsAsFactors = FALSE)
  if (cats) {
    d$timepoint <- rep(c("d0", "d3", "d7"), length.out = n)
    d$focus     <- rep(c("Pulmonary", "Cutaneous"), length.out = n)
    d$one_level <- rep("only", n)
  }
  d
}

test_that("every category is split by, not only the named group column", {
  out <- withr::local_tempdir()
  f <- basename(cyRAVEN::fig_marker_umaps_by_group(
    marker_cells(), c("CD3", "CD4"), out, group_col = "timepoint"))

  # 2 markers x (1 pooled + timepoint + focus). one_level is single-valued.
  expect_setequal(f, c(
    "umap_CD3.png", "umap_CD3_by_timepoint.png", "umap_CD3_by_focus.png",
    "umap_CD4.png", "umap_CD4_by_timepoint.png", "umap_CD4_by_focus.png"))
  expect_true(all(file.size(file.path(out, f)) > 1000))
})

test_that("a category is drawn even when no group column is named", {
  # The regression this exists for: infection focus went undrawn on a run
  # grouped by timepoint, and on ungrouped runs nothing was split at all.
  out <- withr::local_tempdir()
  f <- basename(cyRAVEN::fig_marker_umaps_by_group(
    marker_cells(), "CD3", out, group_col = NULL))

  expect_true("umap_CD3.png" %in% f)
  expect_true("umap_CD3_by_timepoint.png" %in% f)
  expect_true("umap_CD3_by_focus.png" %in% f)
})

test_that("numeric columns, identifiers and single-level columns are not faceted", {
  out <- withr::local_tempdir()
  f <- basename(cyRAVEN::fig_marker_umaps_by_group(
    marker_cells(), "CD3", out, group_col = "timepoint"))

  expect_false(any(grepl("age_years|patient_id|sample_id|one_level", f)))
})

test_that("with no category at all, only the pooled panel is written", {
  cells <- marker_cells(cats = FALSE)
  for (gc in list(NULL, "nosuchcolumn")) {
    out <- withr::local_tempdir()
    f <- basename(cyRAVEN::fig_marker_umaps_by_group(
      cells, c("CD3", "CD4"), out, group_col = gc))
    expect_setequal(f, c("umap_CD3.png", "umap_CD4.png"))
    expect_false(any(grepl("_by_", f)))
  }
})

test_that("marker_facet_cols keeps categories and drops the rest", {
  cols <- cyRAVEN:::marker_facet_cols(marker_cells(),
                                      feature_cols = c("CD3", "CD4"))
  expect_true(all(c("timepoint", "focus") %in% cols))
  expect_false(any(c("age_years", "patient_id", "sample_id", "one_level") %in% cols))
})

test_that("it returns early rather than writing an empty folder", {
  out <- withr::local_tempdir()
  expect_length(cyRAVEN::fig_marker_umaps_by_group(
    marker_cells(), "NOT_A_MARKER", out, group_col = "timepoint"), 0L)
  expect_length(cyRAVEN::fig_marker_umaps_by_group(
    marker_cells()[0, ], "CD3", out, group_col = "timepoint"), 0L)
})


test_that("the split figure is wider than the pooled one, one panel per level", {
  # NOT a test of the panel border. The border is set inside a local closure that
  # is not reachable from here, and detecting a drawn rectangle in a PNG is not
  # worth the machinery; it was checked by eye when it was added. What this pins
  # is the thing a border change could plausibly break by accident -- that the
  # split figure still gets a canvas sized to its number of levels, rather than
  # the pooled figure's width with the panels squeezed into it.
  out <- withr::local_tempdir()
  d <- marker_cells()
  f <- basename(suppressMessages(cyRAVEN::fig_marker_umaps_by_group(
    d, "CD3", out, group_col = "timepoint")))
  expect_true("umap_CD3.png" %in% f)
  expect_true("umap_CD3_by_timepoint.png" %in% f)

  png_width <- function(p) {
    con <- file(p, "rb"); on.exit(close(con), add = TRUE)
    readBin(con, "raw", 16L)
    readBin(con, "integer", 1L, size = 4L, endian = "big")
  }
  w_pooled <- png_width(file.path(out, "umap_CD3.png"))
  w_split  <- png_width(file.path(out, "umap_CD3_by_timepoint.png"))
  expect_gt(w_split, w_pooled)
})
