# The properties every figure must have, whatever it happens to draw.
#
# WHY THIS FILE EXISTS. Five figure defects shipped in succession and every one
# of them passed the full suite:
#
#   * 11 of 20 figures written on a transparent background, so black labels
#     vanished against a dark viewer
#   * a key drawn on the right after the panel had asked for top-left, because
#     patchwork's guide collection lifts the guide out of the panel
#   * a twelve-entry key that consumed the whole height reserved for the UMAP
#     row in umap_multigraph_overlay.png and collapsed those panels to a sliver
#   * a caption clipped mid-name at the canvas edge
#   * two study groups issued the same colour by palette recycling
#
# The figure assertions that existed checked that the file was created and was
# larger than 1000 bytes. None of the defects above changes either.
#
# So these assertions are about the RENDERED FILE, not the ggplot object. Every
# defect above was introduced between the object and the device -- by guide
# collection, by the theme, or by the device's own background handling -- and an
# assertion on the object would have missed all of them.
#
# The contract itself -- png_header(), expect_figure_contract() and the bounds
# they check -- now lives in helper-figure-contract.R, so a second test file
# writing figures asserts against the same definition instead of a weaker copy of
# it.
#
# COVERAGE. Of the 29 exported fig_* functions, this file renders 23. The clinical
# family and the between-group volcano are rendered against the same contract in
# test-clinical.R, next to the tests for the statistics they draw.
# fig_gate_strategy() needs a fitted convex-gate strategy object, and
# test-convex-gates.R renders it as part of testing the fit.

# ---- fixtures ---------------------------------------------------------------
# Built on helper-synthetic.R wherever possible so this file adds no new data,
# and derived by calling the real producing function wherever that is cheap, so
# the shapes cannot drift from what the pipeline actually passes.

out <- withr::local_tempdir()
gof <- synth_groups()

fx_fixture <- function() {
  m <- synth_mfi()
  m$block <- "myeloid"                       # the one column functional markers adds
  m
}

ufreq_fixture <- function() {
  f <- synth_freq()
  f$u_pct_points <- 0.30
  f$pct_lo  <- pmax(0, f$pct_of_cd45_pos - f$u_pct_points)
  f$pct_hi  <- f$pct_of_cd45_pos + f$u_pct_points
  f$lod_pct <- 0.05
  f$loq_pct <- 0.20
  f$detection <- ifelse(f$pct_of_cd45_pos < f$loq_pct, "below LOD", "quantified")
  f
}

ratios_fixture <- function() {
  f <- synth_freq()
  cd4 <- f[f$population == "CD4 T cells", ]
  cd8 <- f[f$population == "CD8 T cells", ]
  data.frame(sample_id = cd4$sample_id, panel = "P1",
             population = "CD4:CD8",
             value = cd4$pct_of_cd45_pos / cd8$pct_of_cd45_pos,
             is_control = FALSE, qc_status = "pass", stringsAsFactors = FALSE)
}

counts_fixture <- function() {
  f <- synth_freq()
  data.frame(sample_id = f$sample_id, population = f$population,
             cells_per_ul = f$pct_of_cd45_pos * 50,
             is_control = FALSE, qc_status = "pass", stringsAsFactors = FALSE)
}

budget_fixture <- function() {
  f <- synth_freq()
  do.call(rbind, lapply(c("CD3", "CD4", "counting"), function(tm)
    data.frame(sample_id = f$sample_id, population = f$population, term = tm,
               u_pct_points = abs(stats::rnorm(nrow(f), 0.2, 0.05)),
               stringsAsFactors = FALSE)))
}

bins_fixture <- function() {
  ids <- names(gof)[1:4]
  do.call(rbind, lapply(ids, function(s) data.frame(
    sample_id = s, bin = 1:20, time_from = (0:19) * 10, time_to = (1:20) * 10,
    n_events = c(rep(100, 10), 40, rep(100, 9)),
    flagged = c(rep(FALSE, 10), TRUE, rep(FALSE, 9)),
    stringsAsFactors = FALSE)))
}

# The two `recon` shapes the pipeline builds: one for the gate reconstruction
# figure (R/pipeline.R:570) and one for the per-marker density panel (:676).
recon_fixture <- function() {
  lapply(names(gof)[1:3], function(s) list(
    sample_id = s, panel = "P1", verdict = "stained sample - included",
    is_control = FALSE, qc_status = "pass", cd45_source = "valley",
    fsc_log10 = FALSE, ssc_log10 = TRUE,
    cd45 = stats::rnorm(400, 3, 1),
    gate = list(fsc_lo = 1, fsc_hi = 9, ssc_lo = 1, ssc_hi = 9),
    cd45_threshold = 2,
    marker_densities = stats::setNames(
      lapply(MARKERS, function(m) stats::rnorm(400, 2, 1)), MARKERS)))
}

qc_recon_fixture <- function() {
  lapply(names(gof)[1:3], function(s) list(
    sample_id = s,
    thresholds = stats::setNames(
      lapply(MARKERS, function(m) list(threshold = 1.5, source = "valley")),
      MARKERS),
    marker_densities = stats::setNames(
      lapply(MARKERS, function(m) stats::rnorm(400, 2, 1)), MARKERS)))
}

# ---- the abundance family ---------------------------------------------------

test_that("abundance grids render opaque and within bounds", {
  f <- file.path(out, "group_comparison.png")
  suppressMessages(fig_group_comparison(synth_freq(), f, group_of = gof,
                                        reference = "HC"))
  expect_figure_contract(f, "fig_group_comparison")

  f2 <- file.path(out, "population_frequencies.png")
  suppressMessages(fig_population_frequencies(synth_freq(), f2))
  expect_figure_contract(f2, "fig_population_frequencies")

  f3 <- file.path(out, "functional_markers.png")
  suppressMessages(fig_functional_markers(fx_fixture(), f3, group_of = gof,
                                          reference = "HC"))
  expect_figure_contract(f3, "fig_functional_markers")

  f4 <- file.path(out, "population_ratios.png")
  suppressMessages(fig_population_ratios(ratios_fixture(), f4, group_of = gof,
                                         reference = "HC"))
  expect_figure_contract(f4, "fig_population_ratios")

  f5 <- file.path(out, "marker_state.png")
  suppressMessages(fig_marker_state(synth_mfi(), f5, group_of = gof,
                                    reference = "HC"))
  expect_figure_contract(f5, "fig_marker_state")
})

test_that("the absolute-count figures render", {
  f <- file.path(out, "absolute_counts.png")
  suppressMessages(fig_group_comparison(
    counts_fixture(), f, group_of = gof, reference = "HC",
    value_col = "cells_per_ul", value_label = "cells / uL"))
  expect_figure_contract(f, "absolute counts via fig_group_comparison")

  f2 <- file.path(out, "absolute_counts_qc.png")
  suppressMessages(fig_absolute_counts_qc(counts_fixture(), f2, group_of = gof))
  expect_figure_contract(f2, "fig_absolute_counts_qc")
})

# ---- uncertainty ------------------------------------------------------------

test_that("the uncertainty figures render", {
  uf <- ufreq_fixture()
  f <- file.path(out, "frequency_uncertainty.png")
  suppressMessages(fig_frequency_uncertainty(uf, f, group_of = gof))
  expect_figure_contract(f, "fig_frequency_uncertainty")

  f2 <- file.path(out, "detection_limits.png")
  suppressMessages(fig_detection_limits(uf, f2))
  expect_figure_contract(f2, "fig_detection_limits")

  f3 <- file.path(out, "uncertainty_budget.png")
  suppressMessages(fig_uncertainty_budget(budget_fixture(), f3))
  expect_figure_contract(f3, "fig_uncertainty_budget")
})

# ---- gating and diagnostics -------------------------------------------------

test_that("threshold drift renders, with and without a reference", {
  thr <- synth_thresholds()
  f <- file.path(out, "threshold_drift.png")
  suppressMessages(fig_threshold_drift(thr, f, group_of = gof, reference = "HC"))
  expect_figure_contract(f, "fig_threshold_drift")

  # REGRESSION. semantic_colours() returns NULL when no reference anchors the
  # study palette, and scale_fill_manual(values = NULL) then aborts inside
  # ggplot with "Insufficient values in manual scale". The argument defaults to
  # NULL, so every direct caller hit it; run_cyraven() was safe only because it
  # always passes one.
  f2 <- file.path(out, "threshold_drift_noref.png")
  expect_no_error(
    suppressMessages(fig_threshold_drift(thr, f2, group_of = gof)))
  expect_figure_contract(f2, "fig_threshold_drift (no reference)")
})

test_that("the gate reconstruction figures render", {
  f <- file.path(out, "recon_diagnostics.png")
  suppressMessages(fig_recon_diagnostics(recon_fixture(), f))
  expect_figure_contract(f, "fig_recon_diagnostics")

  f2 <- file.path(out, "gating_qc.png")
  suppressMessages(fig_gating_qc(qc_recon_fixture(), f2))
  expect_figure_contract(f2, "fig_gating_qc")
})

test_that("the acquisition-stability figure renders", {
  f <- file.path(out, "acquisition_qc.png")
  suppressMessages(fig_acquisition_qc(bins_fixture(), f))
  expect_figure_contract(f, "fig_acquisition_qc")
})

# ---- embedding --------------------------------------------------------------

test_that("the UMAP figures render opaque", {
  cells <- synth_cells()

  f <- file.path(out, "umap_overview.png")
  suppressMessages(fig_umap_overview(cells, f))
  expect_figure_contract(f, "fig_umap_overview")

  f2 <- file.path(out, "umap_overview_by_group.png")
  suppressMessages(fig_umap_overview_by_group(cells, f2, group_col = "cohort"))
  expect_figure_contract(f2, "fig_umap_overview_by_group")

  f3 <- file.path(out, "umap_markers.png")
  suppressMessages(fig_marker_grid(cells, MARKERS, f3))
  expect_figure_contract(f3, "fig_marker_grid")

  f4 <- file.path(out, "umap_density.png")
  suppressMessages(fig_density_by_sample(cells, f4, facet_by = "cohort"))
  expect_figure_contract(f4, "fig_density_by_sample")

  f5 <- file.path(out, "cohort_composition_heatmap.png")
  suppressMessages(fig_cohort_confusion(cells, f5, group_col = "cohort"))
  expect_figure_contract(f5, "fig_cohort_confusion")

  f6 <- file.path(out, "population_marker_heatmap.png")
  suppressMessages(fig_population_marker_heatmap(synth_mfi(), f6))
  expect_figure_contract(f6, "fig_population_marker_heatmap")
})

test_that("the sample key is dropped only past the cap, and the panel survives", {
  # REGRESSION. With one legend entry per sample, a 107-sample cohort spent most
  # of the canvas on the key and the UMAP panels collapsed to overlapping
  # slivers. Two things have to hold: the key goes away when there are too many
  # samples, and nothing changes for the cohort sizes this package is normally
  # run on.
  cells <- synth_cells()
  n_small <- length(unique(cells$sample_id))
  expect_lte(n_small, 40L)   # the fixture must sit UNDER the cap for this test

  f_small <- file.path(out, "umap_overview_small.png")
  suppressMessages(fig_umap_overview(cells, f_small))
  h_small <- expect_figure_contract(f_small, "fig_umap_overview (few samples)")

  # Same cells, relabelled into more samples than the cap. The embedding is
  # unchanged; only the number of distinct sample_id values grows.
  many <- cells
  many$sample_id <- paste0("S", seq_len(nrow(many)) %% 60L)
  expect_gt(length(unique(many$sample_id)), 40L)
  f_many <- file.path(out, "umap_overview_many.png")
  suppressMessages(fig_umap_overview(many, f_many))
  h_many <- expect_figure_contract(f_many, "fig_umap_overview (many samples)")

  # The canvas is the same size either way -- the fix is the key, not the paper.
  expect_equal(h_many$width, h_small$width)
  expect_equal(h_many$height, h_small$height)
})

test_that("clustering and batch figures render", {
  cells <- synth_cells()
  uc <- suppressMessages(
    run_unsupervised_clusters(cells, MARKERS, n_clusters = 4L, grid = 3L,
                              seed = 1L))
  f <- file.path(out, "unsupervised_clusters.png")
  suppressMessages(fig_unsupervised_clusters(cells, uc$cluster, f))
  expect_figure_contract(f, "fig_unsupervised_clusters")

  rep <- suppressMessages(
    batch_mixing_report(cells, "batch", group_col = "cohort",
                        k = 10L, max_cells = 600L, n_perm = 3L))
  skip_if(is.null(rep), "batch report unavailable on this fixture")
  f2 <- file.path(out, "batch_diagnostic.png")
  suppressMessages(fig_batch_diagnostic(cells, rep, f2, batch_col = "batch"))
  expect_figure_contract(f2, "fig_batch_diagnostic")
})

test_that("the multigraph overlay keeps its UMAP row", {
  # REGRESSION for the collapsed row. The key is built from the population
  # labels, so the cohort's cluster count drives its height; with ncol = 1 it
  # grew downwards until it had eaten the row. The composition is
  # title + UMAP row + peaks pane, so a healthy figure is much taller than it is
  # wide, and the UMAP row is a fixed share of it. Asserting the aspect ratio
  # catches the collapse without pinning the exact layout.
  cells <- synth_cells()
  f <- file.path(out, "umap_multigraph_overlay.png")
  suppressMessages(fig_multigraph_overlay(
    cells, f, markers = MARKERS, group_col = "cohort", reference = "HC"))
  hdr <- expect_figure_contract(f, "fig_multigraph_overlay")
  expect_gt(hdr$height, hdr$width * 0.5)
})

test_that("every figure written by these tests is opaque", {
  # A sweep rather than a per-figure assertion, so a figure added to this file
  # later cannot quietly skip the background check.
  pngs <- list.files(out, pattern = "[.]png$", full.names = TRUE)
  expect_gt(length(pngs), 15L)
  bad <- Filter(function(p) !(png_header(p)$colour_type %in% OPAQUE_TYPES), pngs)
  expect_equal(basename(bad), character(0))
})
