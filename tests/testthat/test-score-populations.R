# Where a specification becomes numbers.
#
# WHY THIS FILE EXISTS. score_populations() turns the declared gating strategy
# into the cell labels every frequency, figure and test downstream is built
# from, and it had no test. Its `any_of` semantics are also mirrored by
# --check in R/input-check.R, and those two drifted apart once already: the
# checker counted the reserved `any_of` key as a marker name and reported every
# specification carrying a disjunction as naming a marker no file contains.
# Pinning the rules here is what stops the two diverging again.

# 800 cells: half CD3-bright, and within those half CD4-bright. CD16 and CD56
# mark overlapping but different halves, which is what a disjunction is for.
fixture <- function(seed = 5) {
  withr::local_seed(seed)
  n <- 800L
  cd3 <- rep(c(TRUE, FALSE), each = n / 2)
  cd4 <- rep(c(TRUE, FALSE, TRUE, FALSE), each = n / 4)
  tmat <- cbind(
    CD3  = ifelse(cd3, stats::rnorm(n, 5, 0.4), stats::rnorm(n, 0, 0.4)),
    CD4  = ifelse(cd4, stats::rnorm(n, 5, 0.4), stats::rnorm(n, 0, 0.4)),
    CD16 = ifelse(seq_len(n) %% 2 == 0, stats::rnorm(n, 5, 0.4),
                  stats::rnorm(n, 0, 0.4)),
    CD56 = ifelse(seq_len(n) %% 3 == 0, stats::rnorm(n, 5, 0.4),
                  stats::rnorm(n, 0, 0.4)))
  list(tmat = tmat, thr = c(CD3 = 2.5, CD4 = 2.5, CD16 = 2.5, CD56 = 2.5),
       parent = rep(TRUE, n), n = n)
}

test_that("a declared population is scored, and unmatched cells fall to Other", {
  f <- fixture()
  s <- suppressMessages(score_populations(
    f$tmat, f$thr, f$parent, list(`T cells` = list(CD3 = "above"))))
  expect_length(s$unavailable, 0L)
  expect_equal(sum(s$masks$`T cells`), f$n / 2, tolerance = 0.02 * f$n)
  expect_setequal(unique(s$labels), c("T cells", "Other CD45+"))
})

test_that("any_of is a disjunction: either member satisfies it", {
  f <- fixture()
  spec <- list(NK = list(CD3 = "below",
                         any_of = list(CD16 = "above", CD56 = "above")))
  s <- suppressMessages(score_populations(f$tmat, f$thr, f$parent, spec))
  expect_length(s$unavailable, 0L)

  # Every scored cell is CD3-negative and positive for at least one of the two.
  sel <- s$masks$NK
  expect_true(all(f$tmat[sel, "CD3"] < f$thr[["CD3"]]))
  expect_true(all(f$tmat[sel, "CD16"] > f$thr[["CD16"]] |
                  f$tmat[sel, "CD56"] > f$thr[["CD56"]]))

  # And it is a union, not an intersection: strictly more cells than requiring
  # both. This is the whole reason the key exists.
  both <- list(NK = list(CD3 = "below", CD16 = "above", CD56 = "above"))
  s2 <- suppressMessages(score_populations(f$tmat, f$thr, f$parent, both))
  expect_gt(sum(s$masks$NK), sum(s2$masks$NK))
})

test_that("a disjunction survives losing one member", {
  # CD16+ OR CD56+ is still a usable gate when only CD56 is in the panel. The
  # population must narrow, not disappear.
  f <- fixture()
  tmat <- f$tmat[, c("CD3", "CD4", "CD56")]
  spec <- list(NK = list(CD3 = "below",
                         any_of = list(CD16 = "above", CD56 = "above")))
  s <- suppressMessages(score_populations(
    tmat, f$thr[c("CD3", "CD4", "CD56")], f$parent, spec))
  expect_length(s$unavailable, 0L)
  expect_gt(sum(s$masks$NK), 0L)
  expect_true(all(tmat[s$masks$NK, "CD56"] > f$thr[["CD56"]]))
})

test_that("a disjunction with no member left makes the population unavailable", {
  f <- fixture()
  tmat <- f$tmat[, c("CD3", "CD4")]
  spec <- list(NK = list(CD3 = "below",
                         any_of = list(CD16 = "above", CD56 = "above")))
  s <- suppressMessages(score_populations(tmat, f$thr[c("CD3", "CD4")],
                                          f$parent, spec))
  expect_named(s$unavailable, "NK")
  expect_match(s$unavailable$NK, "no any_of marker available")
  expect_null(s$masks$NK)
})

test_that("a missing DIRECT marker makes the population unavailable", {
  f <- fixture()
  spec <- list(`T cells` = list(CD3 = "above", CD8 = "above"))
  s <- suppressMessages(score_populations(f$tmat, f$thr, f$parent, spec))
  expect_named(s$unavailable, "T cells")
  expect_match(s$unavailable$`T cells`, "marker\\(s\\) not in panel")
  expect_match(s$unavailable$`T cells`, "CD8")
})

test_that("a non-finite threshold counts as a missing marker", {
  f <- fixture()
  thr <- f$thr; thr[["CD4"]] <- NA_real_
  spec <- list(`CD4 T` = list(CD3 = "above", CD4 = "above"))
  s <- suppressMessages(score_populations(f$tmat, thr, f$parent, spec))
  expect_named(s$unavailable, "CD4 T")
  expect_match(s$unavailable$`CD4 T`, "CD4")
})

test_that("the most specific definition wins a cell", {
  # A cell satisfying both must be labelled by the deeper one, or the shallow
  # population silently absorbs the specific one.
  f <- fixture()
  spec <- list(`T cells`    = list(CD3 = "above"),
               `CD4 T cells` = list(CD3 = "above", CD4 = "above"))
  s <- suppressMessages(score_populations(f$tmat, f$thr, f$parent, spec))
  deep <- s$masks$`CD4 T cells`
  expect_true(all(s$labels[deep] == "CD4 T cells"))
  expect_true(any(s$labels == "T cells"))
})

test_that("an any_of group counts as ONE requirement for specificity", {
  # Otherwise a loose disjunction outranks a tight conjunction and takes its
  # cells. Two AND terms must beat one AND term plus a disjunction.
  f <- fixture()
  spec <- list(
    loose = list(CD3 = "above", any_of = list(CD16 = "above", CD56 = "above")),
    tight = list(CD3 = "above", CD4 = "above", CD16 = "above"))
  s <- suppressMessages(score_populations(f$tmat, f$thr, f$parent, spec))
  overlap <- s$masks$loose & s$masks$tight
  skip_if(!any(overlap), "fixture produced no overlapping cells")
  expect_true(all(s$labels[overlap] == "tight"))
})

test_that("an unknown direction is refused rather than guessed", {
  f <- fixture()
  expect_error(
    suppressMessages(score_populations(
      f$tmat, f$thr, f$parent, list(X = list(CD3 = "high")))),
    "unknown direction")
})

test_that("intermediate without an upper bound is unavailable, not silently wrong", {
  f <- fixture()
  spec <- list(`CD4 int` = list(CD4 = "intermediate"))
  s <- suppressMessages(score_populations(f$tmat, f$thr, f$parent, spec))
  expect_named(s$unavailable, "CD4 int")
  expect_match(s$unavailable$`CD4 int`, "no upper bound")

  s2 <- suppressMessages(score_populations(f$tmat, f$thr, f$parent, spec,
                                           hi_thr = list(CD4 = 4)))
  expect_length(s2$unavailable, 0L)
  sel <- s2$masks$`CD4 int`
  expect_true(all(f$tmat[sel, "CD4"] > 2.5 & f$tmat[sel, "CD4"] < 4))
})

test_that("the parent mask bounds every population", {
  f <- fixture()
  parent <- rep(c(TRUE, FALSE), length.out = f$n)
  s <- suppressMessages(score_populations(
    f$tmat, f$thr, parent, list(`T cells` = list(CD3 = "above"))))
  expect_true(all(parent[s$masks$`T cells`]))
  expect_true(all(is.na(s$labels[!parent])))
})
