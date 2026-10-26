# =============================================================================
# MFI summary and the composition of the embedding
#
# The claims: the cohort MFI is the median over donors, back-transformed into
# instrument units through the run's own transform; controls and QC failures
# are left out of it; the group rows carry a difference planted in one group;
# the composition shares sum to 100 with the catch-all listed last; and every
# figure satisfies the rendered-file contract.
# =============================================================================

test_that("defining markers read as a phenotype string", {
  spec <- list(
    "CD4 T" = list(CD3 = "above", CD4 = "above", CD127 = "below"),
    "NK" = list(CD3 = "below", any_of = list(CD16 = "above", CD56 = "above")),
    "Int mono" = list(CD14 = "intermediate", CD16 = "bright"))
  ph <- population_phenotype(spec)
  expect_identical(unname(ph), c("CD3+ CD4+ CD127-", "CD3- (CD16+ or CD56+)",
                                 "CD14int CD16hi"))
  expect_identical(names(ph), names(spec))
})

test_that("back-transformation inverts the run's transform", {
  tr <- make_transform("arcsinh", cofactor = 150)
  x <- c(-40, 0, 300, 25000)
  expect_equal(mfi_back_transform(tr$fn(x), rep("CD3", 4), tr), x)
  # Scatter that explore mode put on log10 comes back through 10^.
  expect_equal(mfi_back_transform(c(4, 5), c("SSC-A", "SSC-A"), tr,
                                  scatter = "SSC-A"), c(1e4, 1e5))
  expect_identical(format_mfi(c(840, 12345, 250000, 1.5e6, NA)),
                   c("840", "12.3k", "250k", "1.5M", ""))
})

test_that("the cohort MFI is the median over donors, in instrument units", {
  mfi <- synth_mfi()
  tr <- make_transform("arcsinh", cofactor = 150)
  ms <- summarise_mfi(mfi, transforms = list(P1 = tr), group_of = synth_groups())
  expect_setequal(unique(ms$group), c("All samples", "HC", "S1", "S2"))
  one <- ms[ms$group == "All samples" & ms$population == "B cells" &
              ms$marker == "CD19", ]
  v <- mfi$median_asinh[mfi$population == "B cells" & mfi$marker == "CD19"]
  expect_equal(one$median_transformed, round(stats::median(v), 4))
  expect_equal(one$mfi, signif(sinh(stats::median(v)) * 150, 5))
  expect_identical(one$n_samples, 16L)
  # PLANTED: CD4 in CD4 T cells is brighter in S1 only.
  cd4 <- ms[ms$population == "CD4 T cells" & ms$marker == "CD4", ]
  expect_gt(cd4$mfi[cd4$group == "S1"], 2 * cd4$mfi[cd4$group == "HC"])
  expect_lt(abs(log(cd4$mfi[cd4$group == "S2"] / cd4$mfi[cd4$group == "HC"])), 0.5)
})

test_that("controls and staining failures do not enter the cohort MFI", {
  mfi <- synth_mfi()
  mfi$median_asinh[mfi$sample_id == "S01"] <- 9
  mfi$is_control[mfi$sample_id == "S01"] <- TRUE
  tr <- make_transform("arcsinh", cofactor = 150)
  ms <- summarise_mfi(mfi, transforms = tr)
  expect_true(all(ms$n_samples == 15L))
  expect_lt(max(ms$median_transformed), 5)
})

test_that("composition shares sum to 100 with the catch-all last", {
  cells <- synth_cells()
  comp <- composition_table(cells, "population_label", group_col = "cohort")
  expect_equal(sum(comp$pct_of_embedded), 100, tolerance = 1e-3)
  expect_identical(comp$population_label[nrow(comp)], "Other CD45+")
  for (g in c("HC", "S1", "S2"))
    expect_equal(sum(comp[[paste0("pct_", g)]]), 100, tolerance = 1e-3)
  expect_identical(sum(comp$cells), nrow(cells))
  # Weighted shares follow sampling_weight: doubling one label's weight
  # doubles its odds against the rest.
  cells$sampling_weight <- ifelse(cells$population_label == "B cells", 2, 1)
  w <- composition_table(cells, "population_label")
  b0 <- comp$pct_of_embedded[comp$population_label == "B cells"]
  b1 <- w$pct_of_embedded[w$population_label == "B cells"]
  expect_gt(b1, b0 * 1.6)
})

test_that("the MFI and composition figures satisfy the figure contract", {
  d <- withr::local_tempdir()
  tr <- make_transform("arcsinh", cofactor = 150)
  ms <- summarise_mfi(synth_mfi(), transforms = tr, group_of = synth_groups())
  ms$panel <- NULL
  share <- c("CD4 T cells" = 20, "B cells" = 10)
  suppressMessages({
    fig_mfi_heatmap(ms, file.path(d, "h.png"), freq = share)
    fig_mfi_heatmap(ms, file.path(d, "hg.png"), by_group = TRUE)
    fig_mfi_dotplot(ms, file.path(d, "dp.png"))
    cells <- synth_cells()
    comp <- composition_table(cells, "population_label", group_col = "cohort",
                              colours = population_colours(
                                sort(unique(cells$population_label))))
    fig_umap_composition(cells, comp, file.path(d, "uc.png"))
  })
  for (f in c("h.png", "hg.png", "dp.png", "uc.png"))
    expect_figure_contract(file.path(d, f), f)
})

test_that("explore mode writes its MFI and composition outputs", {
  d <- withr::local_tempdir()
  cells <- synth_cells()
  cells$cluster <- paste0("k", match(cells$population_label,
                                     unique(cells$population_label)))
  cells$group <- cells$cohort
  feats <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14")
  ks <- sort(unique(cells$cluster))
  prof <- data.frame(cluster = ks, phenotype = paste0("pheno ", ks),
                     stringsAsFactors = FALSE)
  for (m in feats) prof[[paste0("frac_pos.", m)]] <- 0.5
  tr <- make_transform("arcsinh", cofactor = 150)
  w <- suppressMessages(explore_mfi_composition(cells, feats, prof, tr = tr,
                                                ex_dir = d))
  expect_true(all(c("explore_cluster_composition.csv", "explore_cluster_mfi.csv",
                    "explore_umap_composition.png",
                    "explore_cluster_mfi_heatmap.png",
                    "explore_cluster_mfi_dotplot.png",
                    "explore_cluster_mfi_heatmap_by_group.png") %in% w))
  ms <- utils::read.csv(file.path(d, "explore_cluster_mfi.csv"))
  # No thresholds were given, so percent positive falls back to the pooled call.
  expect_true(all(ms$pct_positive_median[ms$group == "All samples"] == 50))
  comp <- utils::read.csv(file.path(d, "explore_cluster_composition.csv"))
  expect_equal(sum(comp$pct_of_embedded), 100, tolerance = 1e-3)
})

test_that("percent positive per cell type is read from each cell's own call", {
  lab <- c("A", "A", "A", "A", "B", "B")
  pos <- cbind(CD3 = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
               CD19 = c(FALSE, FALSE, FALSE, FALSE, TRUE, NA))
  mp <- composition_marker_pct(lab, pos)
  expect_equal(mp["A", "CD3"], 75)
  expect_equal(mp["B", "CD3"], 0)
  # A cell with no call does not count against the cell type.
  expect_equal(mp["B", "CD19"], 100)
  # Weights follow sampling_weight.
  expect_equal(composition_marker_pct(lab, pos, w = c(3, 1, 1, 1, 1, 1))["A", "CD3"],
               100 * 5 / 6)
})

test_that("calls come from each sample's own threshold, at each cell's event", {
  pops <- list(
    s1 = list(tmat = cbind(CD3 = c(0.5, 2.5, 3.0)), thresholds = list(CD3 = 1)),
    s2 = list(tmat = cbind(CD3 = c(0.5, 2.5, 3.0)), thresholds = list(CD3 = 2.8)))
  cells <- data.frame(sample_id = c("s1", "s1", "s2", "s2"),
                      event_index = c(1L, 2L, 2L, 3L))
  pos <- cell_marker_positivity(cells, pops, c("CD3", "CD19"))
  expect_identical(pos[, "CD3"], c(FALSE, TRUE, FALSE, TRUE))
  # A marker with no threshold has no calls, rather than all negative ones.
  expect_true(all(is.na(pos[, "CD19"])))
})

test_that("a per-timepoint table keeps the pooled order and lists absent types at 0", {
  cells <- synth_cells()
  pooled <- composition_table(cells, "population_label")
  sub <- cells[cells$population_label != "B cells", ]
  mp <- composition_marker_pct(sub$population_label,
                               cbind(CD3 = sub$population_label == "CD4 T cells"))
  tp <- composition_table(sub, "population_label",
                          levels = pooled$population_label, marker_pct = mp)
  expect_identical(tp$population_label, pooled$population_label)
  expect_equal(tp$pct_of_embedded[tp$population_label == "B cells"], 0)
  expect_equal(sum(tp$pct_of_embedded), 100, tolerance = 1e-3)
  expect_equal(tp$pct_pos_CD3[tp$population_label == "CD4 T cells"], 100)
  expect_true(is.na(tp$pct_pos_CD3[tp$population_label == "B cells"]))
  # The marker columns are not mistaken for study-group columns in the figure.
  d <- withr::local_tempdir()
  tp$colour <- population_colours(tp$population_label)
  suppressMessages(fig_umap_composition(sub, tp, file.path(d, "tp.png"),
                                        xlim = range(cells$umap_1),
                                        ylim = range(cells$umap_2)))
  expect_figure_contract(file.path(d, "tp.png"), "per-timepoint composition")
})
