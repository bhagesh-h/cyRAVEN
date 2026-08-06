# =============================================================================
# Conformance of a run against a baseline written from an earlier one
# =============================================================================

test_that("a specification fingerprint ignores order and catches redefinition", {
  a <- list(`CD4 T cells` = list(CD3 = "above", CD4 = "above", CD8 = "below"))
  b <- list(`CD4 T cells` = list(CD8 = "below", CD3 = "above", CD4 = "above"))
  c2 <- list(`CD4 T cells` = list(CD3 = "above", CD4 = "above"))
  expect_identical(spec_fingerprint(a), spec_fingerprint(b))
  expect_false(identical(spec_fingerprint(a), spec_fingerprint(c2)))
})

test_that("a baseline round-trips and a run compared against itself passes", {
  thr <- synth_thresholds()
  fq  <- synth_freq()
  spec <- list(`CD4 T cells` = list(CD3 = "above", CD4 = "above"),
               `B cells` = list(CD19 = "above"))
  p <- withr::local_tempfile(fileext = ".rds")

  write_spec_baseline(p, thr, fq, spec)
  bl <- read_spec_baseline(p)
  expect_identical(bl$spec, spec_fingerprint(spec))

  cf <- specification_conformance(thr, fq, spec, bl)
  expect_identical(cf$summary$n_fail, 0L)
  expect_null(cf$spec_changes)
  # A run against its own baseline has moved nowhere.
  expect_true(all(cf$markers$robust_z[is.finite(cf$markers$robust_z)] < 1e-6))
})

test_that("a cohort that moved as a whole fails, where the peer check cannot", {
  # THE CASE THIS FEATURE EXISTS FOR. Every sample's CD3 cut shifts by the same
  # amount. Within the run each sample still agrees with its peers, because the
  # peers moved too; only a fixed external reference sees it.
  thr <- synth_thresholds()
  fq  <- synth_freq()
  spec <- list(`B cells` = list(CD19 = "above"))
  p <- withr::local_tempfile(fileext = ".rds")
  write_spec_baseline(p, thr, fq, spec)

  drifted <- thr
  drifted$threshold[drifted$marker == "CD3"] <-
    drifted$threshold[drifted$marker == "CD3"] + 3
  cf <- specification_conformance(drifted, fq, spec, read_spec_baseline(p))

  cd3 <- cf$markers[cf$markers$marker == "CD3", ]
  expect_identical(cd3$verdict, "fail")
  expect_gt(cf$summary$n_fail, 0L)
  # Markers that did not move are not swept up with it.
  expect_true(all(cf$markers$verdict[cf$markers$marker != "CD3"] == "pass"))
})

test_that("a changed transform is reported once, not as fifty drifted markers", {
  thr <- synth_thresholds(); fq <- synth_freq()
  spec <- list(`B cells` = list(CD19 = "above"))
  p <- withr::local_tempfile(fileext = ".rds")
  write_spec_baseline(p, thr, fq, spec, opt = list(transform = "arcsinh"))

  cf <- specification_conformance(thr, fq, spec, read_spec_baseline(p),
                                  transform = "logicle")
  expect_true(cf$summary$transform_changed)
  expect_true(all(cf$markers$verdict == "not comparable"))
  expect_identical(cf$summary$n_fail, 0L)
})

test_that("a redefined population is not comparable rather than drifting", {
  thr <- synth_thresholds(); fq <- synth_freq()
  old <- list(`CD4 T cells` = list(CD3 = "above", CD4 = "above"))
  new <- list(`CD4 T cells` = list(CD3 = "above", CD4 = "above", CD8 = "below"),
              `NK cells` = list(CD56 = "above"))
  p <- withr::local_tempfile(fileext = ".rds")
  write_spec_baseline(p, thr, fq, old)

  cf <- specification_conformance(thr, fq, new, read_spec_baseline(p))
  expect_true(cf$summary$spec_changed)
  expect_true("CD4 T cells" %in% cf$spec_changes$population)
  expect_identical(cf$spec_changes$change[cf$spec_changes$population == "CD4 T cells"],
                   "redefined")
  expect_true("NK cells" %in% cf$spec_changes$population[cf$spec_changes$change == "added"])
  # Its frequency comparison is withdrawn rather than reported as drift.
  r <- cf$populations[cf$populations$population == "CD4 T cells", ]
  expect_identical(r$verdict, "not comparable")
  expect_true(is.na(r$robust_z))
})

test_that("reading something that is not a baseline fails clearly", {
  p <- withr::local_tempfile(fileext = ".rds")
  saveRDS(list(a = 1), p)
  expect_error(read_spec_baseline(p), "not a cyRAVEN baseline")
  expect_error(read_spec_baseline(tempfile()), "baseline not found")
})
