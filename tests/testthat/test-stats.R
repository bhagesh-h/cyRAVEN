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

test_that("each paired comparison uses every donor seen at both of its conditions", {
  # 10 donors at d0 and d3, 7 of them also at d7: d3 vs d0 has 10 pairs, d7 vs d0
  # has 7, and Friedman runs on the 7 donors seen at all three.
  ids <- c(sprintf("P%02d_d0", 1:10), sprintf("P%02d_d3", 1:10), sprintf("P%02d_d7", 1:7))
  withr::local_seed(3)
  freq <- data.frame(sample_id = ids, population = "T cells",
                     pct_of_cd45_pos = stats::runif(length(ids), 10, 40),
                     qc_status = "pass", stringsAsFactors = FALSE)
  pair <- stats::setNames(sub("_d.$", "", ids), ids)
  cond <- stats::setNames(sub("^.*_", "", ids), ids)
  ps <- stats_paired_comparison(freq, pair, cond, reference = "d0")
  expect_identical(ps$n_pairs[ps$comparison_condition == "d3"], 10L)
  expect_identical(ps$n_pairs[ps$comparison_condition == "d7"], 7L)
  expect_identical(ps$n_incomplete_pairs_dropped[ps$comparison_condition == "d7"], 3L)
  expect_identical(ps$n_incomplete_pairs_dropped[ps$comparison_condition == "d3"], 0L)
  expect_true(all(is.finite(ps$p_omnibus_friedman)))
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

test_that("the paired effect size is the matched-pairs rank-biserial", {
  # Kerby 2014: (R+ - R-) / (R+ + R-) on the ranks of |d|. For d = 1, 2, 3, -4
  # the ranks are 1..4, so (1 + 2 + 3 - 4) / 10 = 0.2; the share of pairs up
  # minus down, which the column held before, would read 0.5.
  expect_equal(matched_rank_biserial(c(1, 2, 3, -4)), 0.2)
  expect_equal(matched_rank_biserial(c(0, 1, 2)), 1)
  expect_true(is.na(matched_rank_biserial(c(0, 0))))
})

test_that("a group effect does not make an unbalanced covariate a confounder", {
  # Two cohorts far apart in age, a real cohort difference in the outcome, and
  # within each cohort no relation to age. Pooled, age tracks the outcome
  # through the cohort alone; within cohorts it does not.
  sid <- sprintf("S%02d", 1:20); grp <- rep(c("A", "B"), each = 10)
  age <- c(21:30, 51:60)
  within <- c(5, 1, 9, 3, 7, 2, 10, 4, 8, 6)          # unrelated to age order
  y <- c(within, within + 20)
  expect_lt(stats::cor.test(age, y, method = "spearman", exact = FALSE)$p.value, 0.05)
  cf <- suppressMessages(stats_confounding(
    data.frame(sample_id = sid, population = "P", pct_of_cd45_pos = y,
               qc_status = "pass"),
    data.frame(patient_id = sid, age = age),
    data.frame(sample_id = sid, patient_id = sid), setNames(grp, sid),
    covariates = "age"))
  expect_false(any(grepl("^HIGH", cf$confounder_risk)))
  out <- cf[cf$comparison == "covariate vs outcome", ]
  expect_match(out$test, "within cohort")
  expect_gt(out$p_value, 0.05)

  # A real age effect inside each cohort is still a confounder.
  y2 <- c(age[1:10] + within / 10, age[11:20] + within / 10)
  cf2 <- suppressMessages(stats_confounding(
    data.frame(sample_id = sid, population = "P", pct_of_cd45_pos = y2,
               qc_status = "pass"),
    data.frame(patient_id = sid, age = age),
    data.frame(sample_id = sid, patient_id = sid), setNames(grp, sid),
    covariates = "age"))
  expect_true(any(grepl("^HIGH", cf2$confounder_risk)))
})

test_that("a covariate that explains nothing gets p = 1, not a missing row", {
  # anova() leaves Pr(>F) empty when the sums of squares are equal to rounding.
  # Ranks 1..4 split F = {1, 4}, M = {2, 3}: equal means, no sum of squares.
  expect_equal(rank_f_within(c(1, 2, 3, 4), c("F", "M", "M", "F"), rep("A", 4)), 1)
})

test_that("clinical marker association is tested within each population", {
  ids <- names(synth_groups())
  # Monocyte HLA-DR rises with severity, T-cell HLA-DR does not: one test per
  # population and marker keeps them apart, where a median over populations
  # would have mixed them.
  n <- length(ids)
  mfi <- do.call(rbind, lapply(seq_along(ids), function(i) data.frame(
    sample_id = ids[i], population = c("CD4 T cells", "Monocytes"),
    marker = "HLA-DR", median_asinh = c((i * 7) %% n, i),
    qc_status = "pass", stringsAsFactors = FALSE)))
  clin <- list(severity = stats::setNames(seq_len(n), ids))
  mk <- stats_clinical_association(NULL, mfi, clin)$markers
  expect_setequal(mk$population_marker, c("CD4 T cells, HLA-DR", "Monocytes, HLA-DR"))
  mono <- mk[mk$population == "Monocytes", ]
  expect_equal(mono$marker, "HLA-DR")
  expect_equal(mono$estimate, 1)
  expect_lt(abs(mk$estimate[mk$population == "CD4 T cells"]), 1)
})

test_that("the edgeR count model finds a planted change and not the others", {
  skip_if_not_installed("edgeR")
  # Six samples per group, five exclusive labels; label C doubles in group B.
  withr::local_seed(12)
  sids <- sprintf("S%02d", 1:12); grp <- rep(c("A", "B"), each = 6)
  base <- c(A = 4000, B = 2000, C = 300, D = 800, E = 2900)
  rows <- do.call(rbind, lapply(seq_along(sids), function(i) {
    mu <- base * stats::runif(1, 0.8, 1.2)
    if (grp[i] == "B") mu["C"] <- mu["C"] * 2
    data.frame(sample_id = sids[i], population = names(mu),
               count = stats::rnbinom(5, mu = mu, size = 30), qc_status = "pass",
               is_control = FALSE, stringsAsFactors = FALSE)
  }))
  er <- stats_edger_abundance(rows, setNames(grp, sids), reference = "A")
  expect_setequal(er$population, names(base))
  expect_true(er$significant_BH[er$population == "C"])
  expect_gt(er$logFC[er$population == "C"], 0.7)
  expect_lt(sum(er$significant_BH[er$population != "C"]), 2)
})
