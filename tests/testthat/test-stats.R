# =============================================================================
# The statistics: do they find the planted effect, and refuse the wrong ones?
# =============================================================================

test_that("differential state finds the planted marker shift", {
  res <- suppressMessages(
    stats_marker_state(synth_mfi(), synth_groups(), reference = "HC"))
  expect_false(is.null(res))

  m <- res[res$measure == "median_asinh", ]
  hit <- m[m$population == "CD4 T cells" & m$marker == "CD4" &
             m$comparison_group == "S1", ]
  expect_identical(nrow(hit), 1L)

  # The planted shift is +1.5, and delta_median is a DIFFERENCE, not a ratio,
  # median_asinh runs negative, where a ratio is uninterpretable.
  expect_gt(hit$delta_median, 1.0)
  expect_lt(hit$delta_median, 2.0)
  # Complete separation of 6 controls from 5 patients.
  expect_identical(hit$cliffs_delta, 1)

  # NOT which.min(p_value): a rank test returns the SAME minimum p for any
  # perfectly separated pair, whatever the effect size, so at n = 6 vs 5 the
  # smallest p is a large tie and picking one member of it is arbitrary. The
  # planted effect must be IN that tie, and must be the largest by effect size,
  # which is the quantity that actually distinguishes it.
  expect_equal(hit$p_value, min(m$p_value), tolerance = 1e-12)
  expect_identical(which.max(abs(m$delta_median)),
                   which(m$population == "CD4 T cells" & m$marker == "CD4" &
                           m$comparison_group == "S1"))
})

test_that("differential state counts DONORS, not cells", {
  # The whole point. Each sample contributes 500 cells; if n were cells the
  # reference n would be 3000.
  res <- suppressMessages(
    stats_marker_state(synth_mfi(), synth_groups(), reference = "HC"))
  expect_lte(max(res$n_reference), 6L)
  expect_lte(max(res$n_comparison), 5L)
  expect_gt(max(res$cells_reference), 1000)   # cells are reported, just not as n
})

test_that("differential state tests both measures and adjusts within each", {
  res <- suppressMessages(
    stats_marker_state(synth_mfi(), synth_groups(), reference = "HC"))
  expect_setequal(unique(res$measure), c("median_asinh", "pct_positive"))
  expect_true(all(res$p_adj_BH >= res$p_value - 1e-12))
  # Adjusting within each measure, not across both: the two are correlated views
  # of the same cells, so pooling them would over-penalise every hit.
  for (m in unique(res$measure)) {
    i <- res$measure == m
    expect_equal(res$p_adj_BH[i], stats::p.adjust(res$p_value[i], "BH"),
                 tolerance = 1e-12)
  }
})

test_that("controls and QC failures never reach a test", {
  mfi <- synth_mfi()
  mfi$is_control[mfi$sample_id == "S01"] <- TRUE
  mfi$qc_status[mfi$sample_id == "S02"] <- "failed"
  res <- suppressMessages(stats_marker_state(mfi, synth_groups(), reference = "HC"))
  expect_lte(max(res$n_reference), 4L)
})

test_that("a population too small to summarise is dropped, not down-weighted", {
  mfi <- synth_mfi(); mfi$n_cells <- 5L
  expect_null(suppressMessages(
    stats_marker_state(mfi, synth_groups(), reference = "HC")))
})

test_that("CLR is a proper compositional transform", {
  fclr <- clr_frequencies(synth_freq())
  expect_true("clr" %in% names(fclr))
  # A centred log-ratio sums to zero within each composition, by construction.
  sums <- tapply(fclr$clr, fclr$sample_id, sum)
  expect_true(all(abs(sums) < 1e-9))
})

test_that("CLR keeps genuinely-absent populations instead of dropping them", {
  # Zeros are common for rare subsets and are exactly the populations most likely
  # to differ, so dropping them would discard the signal.
  f0 <- synth_freq()
  f0$pct_of_cd45_pos[f0$population == "NK cells"] <- 0
  z <- clr_frequencies(f0)
  expect_identical(sum(z$population == "NK cells"), 16L)
  expect_true(all(is.finite(z$clr)))
})

test_that("concordance separates real effects from composition artefacts", {
  freq <- synth_freq(); grp <- synth_groups()
  raw <- stats_group_comparison(freq, grp, reference = "HC")
  cl  <- stats_group_comparison(clr_frequencies(freq), grp,
                                reference = "HC", value_col = "clr")
  conc <- compositional_concordance(raw, cl)

  expect_identical(nrow(conc), nrow(raw))
  expect_true(all(conc$verdict %in% c(
    "robust_to_composition", "raw_only__possible_composition_artefact",
    "clr_only__was_masked_by_composition", "not_significant_either")))
  # The planted CD4 drop is real, so it must survive the transform.
  expect_true(all(conc$verdict[conc$population == "CD4 T cells"] ==
                    "robust_to_composition"))

  # The closure's fingerprint, asserted on the ESTIMATES rather than on whether
  # they cross an alpha. Only CD4 was changed; every other population was held
  # fixed and merely reapportioned, so on raw percentages they must all appear to
  # move in the OPPOSITE direction to CD4, that is the artefact, while their
  # CLR effect sizes stay far smaller. Asserting "at least one raw_only verdict"
  # instead would be a coin flip at n = 6 vs 5, where significance is decided by
  # a handful of ranks.
  cd4 <- conc[conc$population == "CD4 T cells", ]
  rest <- conc[conc$population != "CD4 T cells", ]
  expect_true(all(cd4$cliffs_delta_raw < 0))
  expect_true(all(rest$cliffs_delta_raw > 0))
  expect_lt(mean(abs(rest$cliffs_delta_clr)), mean(abs(rest$cliffs_delta_raw)))
})

test_that("confounding needs BOTH facts before it cries wolf", {
  cf <- suppressMessages(stats_confounding(
    synth_freq(), synth_patients(), synth_smap(), synth_groups(),
    covariates = c("age", "sex")))
  expect_false(is.null(cf))

  bal <- cf[cf$comparison == "covariate vs cohort" & cf$covariate == "age", ]
  expect_identical(nrow(bal), 1L)
  expect_lt(bal$p_value, 0.05)                       # age IS unbalanced

  hi <- cf[grepl("^HIGH", cf$confounder_risk), ]
  # Anything flagged HIGH must be a covariate that is unbalanced AND associated.
  expect_true(all(hi$comparison == "covariate vs outcome"))
  expect_true(all(hi$p_value < 0.05))
})

test_that("rank ANCOVA labels itself and refuses when it cannot fit", {
  ra <- suppressMessages(stats_rank_ancova(
    synth_freq(), synth_patients(), synth_smap(), synth_groups()))
  expect_true(all(grepl("EXPLORATORY|NOT FITTED", ra$status)))
  # Where it declines, it must say why rather than emit a bare NA.
  nf <- ra[grepl("^NOT FITTED", ra$status), ]
  if (nrow(nf)) expect_true(all(is.na(nf$p_group_adjusted)))
})

test_that("paired analysis counts incomplete pairs rather than hiding them", {
  smap <- synth_smap()
  pair <- stats::setNames(smap$donor, smap$sample_id)
  cond <- stats::setNames(smap$timepoint, smap$sample_id)

  ps <- stats_paired_comparison(synth_freq(), pair, cond)
  expect_false(is.null(ps))
  expect_lte(max(ps$n_pairs), 8L)
  expect_true(all(abs(ps$rank_biserial) <= 1))

  cond["S02"] <- NA                                  # break one pair
  p2 <- stats_paired_comparison(synth_freq(), pair, cond)
  expect_gte(max(p2$n_incomplete_pairs_dropped), 1L)
})

test_that("threshold drift flags the marker whose gate moves, and only that one", {
  td <- stats_threshold_drift(synth_thresholds(), synth_groups(),
                              spec = default_population_spec())
  expect_false(is.null(td))
  cd4 <- td[td$marker == "CD4", ]
  expect_identical(nrow(cd4), 1L)
  expect_match(cd4$drift_flag, "^FLAGGED")
  expect_false(any(grepl("^FLAGGED", td$drift_flag[td$marker != "CD4"])))
})

test_that("a comparison with too few samples returns NULL, not a fake number", {
  freq <- synth_freq()
  g <- synth_groups()[1:3]                            # one group, 3 samples
  expect_null(stats_group_comparison(freq[freq$sample_id %in% names(g), ], g))
})
