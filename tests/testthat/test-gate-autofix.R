# --auto-fix-gates: which hierarchy gates are evidence and which are artefacts.
#
# Two observations trigger a skip, and the second exists because the first was
# not enough. A gate can carry a perfectly respectable "valley" source and still
# be placed somewhere absurd: on the cohort that motivated this, one sample of
# 20 was the only one where a viability minimum was found, and it kept 140 of
# 79,000 cells. Because the embedding equalises cells per sample, that single
# gate throttled all 20 samples to 140 cells and cost 98% of the UMAP. Judging a
# gate on its provenance alone let that through.

test_that("a quantile fallback is skipped: retention is the constant, not the data", {
  r <- gate_needs_autofix("quantile_fallback", enabled = TRUE,
                          kept = 90, parent = 100)
  expect_true(r$skip)
  expect_match(r$reason, "quantile fallback")
})

test_that("a valley that keeps almost nothing is skipped whatever its source says", {
  r <- gate_needs_autofix("valley", enabled = TRUE, kept = 140, parent = 79000)
  expect_true(r$skip)
  expect_match(r$reason, "0.18%")
  expect_match(r$reason, "misplaced")
})

test_that("a valley with plausible retention is left alone", {
  # 95% live is ordinary for PBMC; 30% is a poor preparation but real.
  expect_false(gate_needs_autofix("valley", TRUE, kept = 95, parent = 100)$skip)
  expect_false(gate_needs_autofix("valley", TRUE, kept = 30, parent = 100)$skip)
})

test_that("the floor is far below any plausible value, not near it", {
  # 6% survives, 4% does not: the floor catches broken gates, it does not
  # enforce an expectation about viability.
  expect_false(gate_needs_autofix("valley", TRUE, kept = 6, parent = 100)$skip)
  expect_true(gate_needs_autofix("valley", TRUE, kept = 4, parent = 100)$skip)
})

test_that("both observations at once are reported as both", {
  r <- gate_needs_autofix("quantile_fallback", TRUE, kept = 1, parent = 100)
  expect_true(r$skip)
  expect_match(r$reason, "AND")
})

test_that("without the flag nothing is skipped, but the reason still explains", {
  r <- gate_needs_autofix("valley", enabled = FALSE, kept = 140, parent = 79000)
  expect_false(r$skip)
  expect_match(r$reason, "--auto-fix-gates would skip it")
  # A gate that is fine says nothing at all, so the log stays quiet.
  expect_true(is.na(gate_needs_autofix("valley", FALSE, kept = 95, parent = 100)$reason))
})

test_that("missing counts fall back to judging the source alone", {
  # apply_gate_hierarchy always supplies these, but the function must not
  # invent a verdict from NA if a future caller forgets.
  expect_true(gate_needs_autofix("quantile_fallback", TRUE)$skip)
  expect_false(gate_needs_autofix("valley", TRUE)$skip)
  expect_false(gate_needs_autofix("valley", TRUE, kept = 5, parent = 0)$skip)
})
