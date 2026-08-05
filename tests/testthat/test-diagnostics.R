# =============================================================================
# Batch diagnostics, LISI, figures, and package state
# =============================================================================

test_that("LISI stays within its definitional bounds", {
  cells <- synth_cells()
  co <- as.matrix(cells[, c("umap_1", "umap_2")])
  s <- lisi_score(co, cells$batch, k = 20L)
  # 1 = every neighbour shares the cell's batch; the number of batches = fully
  # mixed. Anything outside that range is a bug, not a finding.
  expect_true(all(s >= 1 - 1e-9))
  expect_true(all(s <= length(unique(cells$batch)) + 1e-9))
})

test_that("partial neighbour selection picks the same set as a full sort", {
  # Optimisation guard: lisi_score() replaced order(dr)[seq_len(k)] with a
  # partial select. Only the SET matters to the Simpson index, but the sets must
  # match exactly, ties included.
  withr::local_seed(8)
  for (trial in 1:20) {
    dr <- c(stats::rnorm(200), rep(0.5, 5))          # deliberate ties
    k <- 10L
    slow <- sort(order(dr)[seq_len(k)])
    kth <- sort.int(dr, partial = k)[k]
    nb <- which(dr <= kth)
    if (length(nb) > k) nb <- nb[order(dr[nb])[seq_len(k)]]
    expect_identical(sort(nb), slow)
  }
})

test_that("batch report detects total batch/cohort confounding", {
  cells <- synth_cells()          # batch is perfectly nested in cohort
  br <- suppressMessages(
    batch_mixing_report(cells, "batch", group_col = "cohort", k = 20L, n_perm = 8L))
  expect_false(is.null(br))
  expect_identical(nrow(br$summary), 1L)
  expect_gte(br$confounding$cramers_v, 0.8)
  expect_match(br$confounding$verdict, "^SEVERE")
})

test_that("batch report restores the RNG stream", {
  cells <- synth_cells()
  set.seed(31); before <- stats::runif(3)
  set.seed(31)
  invisible(suppressMessages(
    batch_mixing_report(cells, "batch", k = 10L, n_perm = 3L)))
  expect_identical(before, stats::runif(3))
})

test_that("a single batch level yields nothing rather than a vacuous figure", {
  cells <- synth_cells(); cells$batch <- "only"
  expect_null(suppressMessages(batch_mixing_report(cells, "batch")))
})

test_that("plot text wraps to the figure width and preserves hard breaks", {
  long <- paste(rep("word", 60), collapse = " ")
  w <- wrap_plot_text(long, width_in = 6)
  expect_gt(length(strsplit(w, "\n")[[1]]), 1L)
  # Explicit newlines in captions are intentional and must survive wrapping.
  two <- wrap_plot_text("short one\nshort two", width_in = 12)
  expect_identical(two, "short one\nshort two")
})

test_that("tile text colour is chosen by contrast, not by a fixed threshold", {
  # Dark fills must take white ink, light fills dark ink. Comparing WCAG contrast
  # ratios keeps this correct if the palette is re-themed.
  cols <- contrast_text_colour(c(0, 50, 100), limits = c(0, 100), option = "D")
  expect_length(cols, 3L)
  expect_identical(cols[1], "white")     # viridis "D" is near-black at 0
  expect_true(cols[3] != "white")        # and bright yellow at 100
})

test_that("expected_positive_markers reads the spec it audits", {
  e <- expected_positive_markers(default_population_spec())
  expect_setequal(e[["CD3 pos CD4 pos"]], c("CD3", "CD4"))
  # "intermediate" is not "above": a population defined by an intermediate band
  # has no required-positive marker and must not be outlined as if it did.
  expect_false("CD14 int CD16 int" %in% names(e))
})

test_that("figures are written and are not blank", {
  od <- withr::local_tempdir()
  mfi <- synth_mfi(); cells <- synth_cells()

  f1 <- file.path(od, "heat.png")
  suppressMessages(fig_population_marker_heatmap(
    mfi, f1, annotate_expected = expected_positive_markers(default_population_spec())))
  expect_true(file.exists(f1)); expect_gt(file.size(f1), 5000)

  f2 <- file.path(od, "conf.png")
  d <- suppressMessages(fig_cohort_confusion(cells, f2, group_col = "cohort"))
  expect_true(file.exists(f2)); expect_gt(file.size(f2), 5000)
  # Each population's cohort shares are a composition and must close to 100%.
  expect_true(all(abs(tapply(d$pct, d$population, sum) - 100) < 1e-6))
})

test_that("the run manifest records what produced a result set", {
  p <- file.path(withr::local_tempdir(), "manifest.txt")
  write_run_manifest(p, opt = list(outdir = ".", seed = 42),
                     files = character(0), status = "completed")
  x <- readLines(p)
  expect_true(any(grepl("^status: completed", x)))
  expect_true(any(grepl("^r_version:", x)))
  expect_true(any(grepl("ggplot2:", x)))
  expect_true(any(grepl("sessionInfo", x)))
})

test_that("package state is settable and restorable", {
  old <- set_fcs_colors(list(gate_highlight = "#0072F0"))
  expect_identical(fcs_colors()$gate_highlight, "#0072F0")
  # A partial override must not blank the keys it omits.
  expect_false(is.null(fcs_colors()$population_palette))
  set_fcs_colors(old)

  oldg <- set_reference_group("HC")
  expect_identical(reference_group(), "HC")
  set_reference_group(oldg)
})

test_that("verbosity is controlled by an option, not by an argument", {
  withr::local_options(cyRAVEN.verbose = "none")
  expect_silent(log_msg("this must not appear"))
  withr::local_options(cyRAVEN.verbose = "inform")
  expect_message(log_msg("this must appear"), "must appear")
})
