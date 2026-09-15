# =============================================================================
# Splitting the figure set by timepoint
#
# The claims: a level is only split out when there is enough of it to draw; the
# levels come back in collection order rather than alphabetical; subsetting is
# by the samples drawn at that visit; and a cross-sectional sheet is left alone
# rather than being split into one-sample figures.
# =============================================================================

sheet <- function(tp, ids = paste0("S", seq_along(tp)))
  data.frame(sample_id = ids, timepoint = tp, stringsAsFactors = FALSE)

test_that("timepoints come back in collection order, not alphabetical", {
  # d10 after d7. Sorted as text it lands between d0 and d3, which puts the
  # sections of the report in an order that misreads as the study design.
  s <- sheet(c("d0", "d0", "d3", "d3", "d7", "d7", "d10", "d10"))
  expect_identical(timepoint_levels(s), c("d0", "d3", "d7", "d10"))
})

test_that("a sheet with no timepoint, or only one level, is not split", {
  expect_null(timepoint_levels(NULL))
  expect_null(timepoint_levels(
    data.frame(sample_id = c("S1", "S2"), stringsAsFactors = FALSE)))
  expect_null(timepoint_levels(sheet(c("d0", "d0", "d0"))))
  # Blank and NA are not levels.
  expect_null(timepoint_levels(sheet(c("d0", "d0", "", NA))))
})

test_that("a level with a single sample is dropped rather than drawn", {
  # One acquisition gives one-column boxplots and a one-sample UMAP: a figure
  # set that looks like a result and carries none.
  s <- sheet(c("d0", "d0", "d3", "d3", "d7"))
  expect_identical(timepoint_levels(s), c("d0", "d3"))
  # And if dropping it leaves fewer than two levels, nothing is split at all.
  expect_null(timepoint_levels(sheet(c("d0", "d0", "d3"))))
})

test_that("samples are selected by the visit they were drawn at", {
  s <- sheet(c("d0", "d3", "d0", "d7"), ids = c("A", "B", "C", "D"))
  expect_identical(samples_at_timepoint(s, "d0"), c("A", "C"))
  expect_identical(samples_at_timepoint(s, "d7"), "D")
  expect_identical(samples_at_timepoint(s, "d99"), character(0))
})

test_that("subsetting a sample-keyed table keeps only those samples", {
  tab <- data.frame(sample_id = c("A", "B", "C"), value = 1:3,
                    stringsAsFactors = FALSE)
  out <- subset_by_sample(tab, c("A", "C"))
  expect_identical(out$sample_id, c("A", "C"))
  expect_identical(out$value, c(1L, 3L))
})

test_that("a subset that empties a table yields NULL, not an empty frame", {
  # Every figure function here treats NULL as "nothing to draw" and an empty
  # frame as "draw an empty figure", so the distinction decides whether a
  # blank PNG appears in the output.
  tab <- data.frame(sample_id = c("A", "B"), value = 1:2,
                    stringsAsFactors = FALSE)
  expect_null(subset_by_sample(tab, "Z"))
  expect_null(subset_by_sample(NULL, "A"))
  expect_null(subset_by_sample(data.frame(x = 1), "A"))   # no sample_id column
})

test_that("a level name becomes a safe directory name", {
  expect_identical(split_dir_name("d0"), "d0")
  expect_identical(split_dir_name("day 0 / pre"), "day_0_pre")
})

test_that("the split writes one directory per timepoint, with figures in each", {
  smap <- sheet(c("d0", "d0", "d0", "d3", "d3", "d7", "d7"),
                ids = paste0("S", 1:7))
  set.seed(4)
  freq <- do.call(rbind, lapply(smap$sample_id, function(s)
    data.frame(sample_id = s,
               population = c("CD4 T cells", "CD8 T cells", "B cells"),
               pct_of_cd45_pos = runif(3, 5, 40),
               count = sample(500:2000, 3),
               stringsAsFactors = FALSE)))
  out <- withr::local_tempdir()
  written <- split_figures_by_timepoint(
    list(smap = smap, freq = freq, colors = fcs_colors()), out)

  root <- file.path(out, "by_timepoint")
  expect_setequal(list.dirs(root, full.names = FALSE, recursive = FALSE),
                  c("d0", "d3", "d7"))
  for (lv in c("d0", "d3", "d7"))
    expect_true(file.exists(file.path(root, lv, "population_frequencies.png")))
  expect_length(written, 3L)
})

test_that("the split draws each visit from its own samples only", {
  # The figure for d7 must not be the pooled figure under a d7 heading. Two
  # visits with deliberately different values give different file contents.
  smap <- sheet(c("d0", "d0", "d7", "d7"), ids = paste0("S", 1:4))
  freq <- data.frame(
    sample_id = rep(paste0("S", 1:4), each = 2),
    population = rep(c("CD4 T cells", "B cells"), 4),
    pct_of_cd45_pos = c(10, 5, 11, 6,    60, 30, 62, 31),
    count = rep(1000L, 8), stringsAsFactors = FALSE)
  out <- withr::local_tempdir()
  split_figures_by_timepoint(
    list(smap = smap, freq = freq, colors = fcs_colors()), out)
  f0 <- file.path(out, "by_timepoint", "d0", "population_frequencies.png")
  f7 <- file.path(out, "by_timepoint", "d7", "population_frequencies.png")
  expect_true(file.exists(f0) && file.exists(f7))
  expect_false(identical(readBin(f0, "raw", file.size(f0)),
                         readBin(f7, "raw", file.size(f7))))
})

test_that("a sheet with nothing to split by writes nothing and does not fail", {
  out <- withr::local_tempdir()
  w <- split_figures_by_timepoint(
    list(smap = sheet(c("d0", "d0", "d0")), freq = NULL), out)
  expect_length(w, 0L)
  expect_false(dir.exists(file.path(out, "by_timepoint")))
})
