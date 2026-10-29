# =============================================================================
# Cuts derived within a parent population
#
# A cut across every CD45+ cell sits wherever the lineages that dominate the
# mixture put a valley. A population declaring `parent:` has its own markers cut
# on the parent's events instead, as a gating template derives each gate on its
# parent (openCyto).
# =============================================================================

test_that("parent conditions are merged in and each marker's scope recorded", {
  spec <- resolve_population_parents(list(
    T = list(CD3 = "above"),
    CD4T = list(parent = "T", CD4 = "above"),
    Treg = list(parent = "CD4T", CD25 = "above", CD127 = "below")))
  expect_identical(spec$Treg, list(CD3 = "above", CD4 = "above",
                                   CD25 = "above", CD127 = "below"))
  sc <- attr(spec, "cut_scope")
  expect_identical(unname(sc$Treg[c("CD3", "CD4", "CD25", "CD127")]),
                   c(NA, "T", "CD4T", "CD4T"))
  expect_identical(attr(spec, "parents")[["Treg"]], "CD4T")
  # no parent anywhere: untouched
  plain <- list(A = list(CD3 = "above"))
  expect_identical(resolve_population_parents(plain), plain)
})

test_that("an undeclared parent, a loop and a contradiction are errors", {
  expect_error(resolve_population_parents(list(A = list(parent = "Z", CD3 = "above"))),
               "does not declare")
  expect_error(resolve_population_parents(list(A = list(parent = "B", CD3 = "above"),
                                               B = list(parent = "A", CD4 = "above"))),
               "loops back")
  expect_error(resolve_population_parents(list(A = list(CD3 = "above"),
                                               B = list(parent = "A", CD3 = "below"))),
               "requires CD3")
})

test_that("a cut within the parent finds the split the CD45+ cut misses", {
  # CD25 across all CD45+ cells: a large negative mode (B cells, monocytes and
  # the CD25- T cells) and a CD25-dim shoulder on activated conventional T
  # cells. Within CD4 T cells the CD25-high regulatory cells form their own mode.
  withr::local_seed(5)
  n_oth <- 30000; n_cd4 <- 10000
  cd3 <- c(rnorm(n_oth, 0.3, 0.3), rnorm(n_cd4, 3, 0.3))
  cd4 <- c(rnorm(n_oth, 0.3, 0.3), rnorm(n_cd4, 3, 0.3))
  cd25_oth <- rnorm(n_oth, 0.4, 0.35)
  cd25_cd4 <- c(rnorm(8800, 0.6, 0.3), rnorm(1200, 2.6, 0.25))
  tm <- cbind(CD3 = cd3, CD4 = cd4, CD25 = c(cd25_oth, cd25_cd4))
  par <- rep(TRUE, nrow(tm))
  spec <- resolve_population_parents(list(
    `CD4 T` = list(CD3 = "above", CD4 = "above"),
    Treg = list(parent = "CD4 T", CD25 = "above")))
  thr <- c(CD3 = 1.6, CD4 = 1.6, CD25 = 1.0)
  sc <- derive_scoped_cuts(tm, thr, par, spec)
  expect_named(sc, "CD4 T")
  expect_identical(sc[["CD4 T"]]$n_parent, sum(cd3 > 1.6 & cd4 > 1.6))
  cut <- sc[["CD4 T"]]$thr[["CD25"]]
  expect_gt(cut, 1.2); expect_lt(cut, 2.2)
  sp <- score_populations(tm, thr, par, spec, scoped = sc)
  # about the 1,200 planted regulatory cells, and not the CD25-dim tail the
  # CD45+ cut of 1.0 lets in
  expect_equal(sum(sp$masks$Treg), 1200, tolerance = 0.1)
  sp_glob <- score_populations(tm, thr, par, spec)
  expect_gt(sum(sp_glob$masks$Treg), 1.3 * sum(sp$masks$Treg))
  # the deeper label wins
  expect_identical(sp$labels[which(sp$masks$Treg)[1]], "Treg")
})

test_that("a parent too small to cut on leaves its children unavailable", {
  tm <- cbind(CD3 = c(rep(0, 990), rep(3, 10)), CD25 = rnorm(1000))
  spec <- resolve_population_parents(list(T = list(CD3 = "above"),
                                          Treg = list(parent = "T", CD25 = "above")))
  sc <- derive_scoped_cuts(tm, c(CD3 = 1, CD25 = 0), rep(TRUE, 1000), spec)
  expect_true(is.na(sc$T$thr[["CD25"]]))
  sp <- suppressMessages(score_populations(tm, c(CD3 = 1, CD25 = 0), rep(TRUE, 1000),
                                           spec, scoped = sc))
  expect_match(sp$unavailable$Treg, "within T")
})

test_that("the built-in regulatory T cell gate is cut within CD4 T cells", {
  sp <- default_population_spec()
  tr <- sp[["CD3 pos CD4 pos CD127 lo CD25 hi"]]
  expect_identical(tr[c("CD3", "CD4")], list(CD3 = "above", CD4 = "above"))
  expect_identical(attr(sp, "cut_scope")[["CD3 pos CD4 pos CD127 lo CD25 hi"]][["CD25"]],
                   "CD3 pos CD4 pos")
})
