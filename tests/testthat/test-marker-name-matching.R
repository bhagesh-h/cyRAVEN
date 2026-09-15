# The spelling a specification must use, and what --check says about it.
#
# This exists because the check previously asserted the opposite of what the run
# does, and told users to rewrite TCR-Vd1 as TCR.Vd1 -- the one spelling that
# matches nothing. The property under test is not "the message reads well", it
# is that the CHECK and the RUN agree about which spelling scores.

test_that("score_populations matches the marker symbol verbatim, hyphens included", {
  # colnames(tmat) are the resolved $PnS symbols; pipeline.R sets them with
  # colnames(tmat) <- avail, and avail comes from names(rd$marker_cols).
  tmat <- cbind(`CD45` = c(5, 5, 1), `CD3` = c(5, 5, 1), `TCR-Vd1` = c(5, 1, 1))
  thr  <- c(CD45 = 3, CD3 = 3, `TCR-Vd1` = 3)
  spec <- list(`Vd1 T cells` = list(CD45 = "above", CD3 = "above",
                                    `TCR-Vd1` = "above"))

  got <- score_populations(tmat, thr, parent = rep(TRUE, 3), spec = spec)
  expect_true(is.null(got$unavailable[["Vd1 T cells"]]))
  expect_equal(sum(got$masks[["Vd1 T cells"]]), 1L)
})

test_that("the syntactic spelling is the one that fails to score", {
  tmat <- cbind(`CD45` = c(5, 5), `CD3` = c(5, 5), `TCR-Vd1` = c(5, 1))
  thr  <- c(CD45 = 3, CD3 = 3, `TCR-Vd1` = 3)
  # A config written as make.names() would have it.
  spec <- list(`Vd1 T cells` = list(CD45 = "above", CD3 = "above",
                                    `TCR.Vd1` = "above"))

  got <- score_populations(tmat, thr, parent = rep(TRUE, 2), spec = spec)
  expect_false(is.null(got$unavailable[["Vd1 T cells"]]))
  expect_match(got$unavailable[["Vd1 T cells"]], "not in panel")
})

test_that("any_of drops the member this panel lacks and keeps the population", {
  # The two gamma-delta spellings across one cohort: no file carries both, so
  # any_of resolves to whichever is present and cannot double-count.
  spec <- list(`Vd1 T cells` = list(
    CD45 = "above", CD3 = "above",
    any_of = list(`TCR-Vd1` = "above", `Vd1` = "above")))

  # A file spelling it TCR-Vd1
  t1 <- cbind(CD45 = c(5, 5, 1), CD3 = c(5, 5, 1), `TCR-Vd1` = c(5, 1, 1))
  r1 <- score_populations(t1, c(CD45 = 3, CD3 = 3, `TCR-Vd1` = 3),
                          rep(TRUE, 3), spec)
  expect_true(is.null(r1$unavailable[["Vd1 T cells"]]))
  expect_equal(sum(r1$masks[["Vd1 T cells"]]), 1L)

  # A file spelling it Vd1 -- same population, same answer
  t2 <- cbind(CD45 = c(5, 5, 1), CD3 = c(5, 5, 1), `Vd1` = c(5, 1, 1))
  r2 <- score_populations(t2, c(CD45 = 3, CD3 = 3, `Vd1` = 3),
                          rep(TRUE, 3), spec)
  expect_true(is.null(r2$unavailable[["Vd1 T cells"]]))
  expect_equal(sum(r2$masks[["Vd1 T cells"]]), 1L)
})

test_that("any_of on a single present member equals the direct requirement", {
  # The guarantee that adding any_of cannot move the original cohort's numbers:
  # where only one member exists, the two specs must score identically.
  tmat <- cbind(CD45 = c(5, 5, 1, 5), CD3 = c(5, 5, 1, 1), `TCR-Vd1` = c(5, 1, 1, 5))
  thr  <- c(CD45 = 3, CD3 = 3, `TCR-Vd1` = 3)
  direct <- list(P = list(CD45 = "above", CD3 = "above", `TCR-Vd1` = "above"))
  anyof  <- list(P = list(CD45 = "above", CD3 = "above",
                          any_of = list(`TCR-Vd1` = "above", `Vd1` = "above")))

  a <- score_populations(tmat, thr, rep(TRUE, 4), direct)
  b <- score_populations(tmat, thr, rep(TRUE, 4), anyof)
  expect_identical(unname(a$masks[["P"]]), unname(b$masks[["P"]]))
})

test_that("any_of is unavailable only when every member is gone", {
  spec <- list(P = list(CD45 = "above",
                        any_of = list(`TCR-Vd1` = "above", `Vd1` = "above")))
  tmat <- cbind(CD45 = c(5, 1))
  got <- score_populations(tmat, c(CD45 = 3), rep(TRUE, 2), spec)
  expect_false(is.null(got$unavailable[["P"]]))
})
