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
