# =============================================================================
# Fluorescence-minus-one reference controls
#
# The claim: a per-marker control can be declared, it is used in preference to an
# unstained tube for the markers it covers, and the run reports how far the
# derived cut sits from it in units the reader can act on.
# =============================================================================

SM <- data.frame(
  file = c("s1.fcs", "s2.fcs", "fmo_ccr7.fcs", "fmo_cd45ra.fcs", "unst.fcs"),
  sample_id = c("S1", "S2", "FMO_CCR7", "FMO_CD45RA", "UNST"),
  fmo_for = c(NA, NA, "CCR7", "CD45RA, CD27", NA),
  is_control = c(FALSE, FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE)

test_that("the fmo_for column is parsed into one row per controlled marker", {
  m <- parse_fmo_map(SM)
  expect_s3_class(m, "data.frame")
  expect_identical(nrow(m), 3L)
  expect_setequal(m$marker, c("CCR7", "CD45RA", "CD27"))
  # One file may control several channels, and whitespace after the comma is
  # what a person actually types.
  expect_identical(sort(m$sample_id[m$marker %in% c("CD45RA", "CD27")]),
                   c("FMO_CD45RA", "FMO_CD45RA"))
})

test_that("a sample map with no FMO column or no entries yields nothing", {
  expect_null(parse_fmo_map(NULL))
  expect_null(parse_fmo_map(SM[, c("file", "sample_id")]))
  blank <- SM; blank$fmo_for <- NA
  expect_null(parse_fmo_map(blank))
  empty <- SM; empty$fmo_for <- ""
  expect_null(parse_fmo_map(empty))
})

test_that("a control is selected for the markers it covers and no others", {
  m <- parse_fmo_map(SM)
  expect_identical(fmo_for_sample(m, "S1", "CCR7"), "FMO_CCR7")
  expect_identical(fmo_for_sample(m, "S1", "CD27"), "FMO_CD45RA")
  # A marker nobody controls falls through, which is what keeps every other
  # marker on the derived path.
  expect_true(is.na(fmo_for_sample(m, "S1", "CD3")))
  expect_true(is.na(fmo_for_sample(NULL, "S1", "CCR7")))
})

test_that("a control never controls itself", {
  # Its own channel is the one left out, so its distribution there is the
  # negative population by construction and cannot also be the sample.
  m <- parse_fmo_map(SM)
  expect_true(is.na(fmo_for_sample(m, "FMO_CCR7", "CCR7")))
})

test_that("control_group confines a control to the batch it was acquired in", {
  # Reagent lots and instrument settings change between batches, so an FMO from
  # one batch is not the negative population of another.
  sm <- SM
  sm$control_group <- c("b1", "b2", "b1", NA, NA)
  m <- parse_fmo_map(sm)
  grp <- c(S1 = "b1", S2 = "b2", FMO_CCR7 = "b1")

  expect_identical(fmo_for_sample(m, "S1", "CCR7", grp), "FMO_CCR7")
  # S2 is in b2 and the only CCR7 control is b1, so it gets none rather than a
  # control from the wrong batch.
  expect_true(is.na(fmo_for_sample(m, "S2", "CCR7", grp)))
  # A control with no group still applies everywhere.
  expect_identical(fmo_for_sample(m, "S2", "CD45RA", grp), "FMO_CD45RA")
})

test_that("the source records which control the cut came from", {
  withr::local_seed(2)
  x <- c(stats::rnorm(8000, 0, 0.6), stats::rnorm(8000, 5, 0.6))
  flat <- stats::runif(20000)          # no valley, so the control decides
  ctrl <- stats::rnorm(5000, 0.5, 0.3)

  unst <- resolve_threshold("M", flat, control_x = ctrl)
  fmo  <- resolve_threshold("M", flat, control_x = ctrl, control_kind = "fmo_q995")
  expect_identical(unst$source, "control_q995")
  expect_identical(fmo$source, "fmo_q995")
  # Same arithmetic; the two are named apart because they are different
  # experiments, not because they compute differently.
  expect_equal(unst$threshold, fmo$threshold)

  # And the rejection branch carries the kind through too.
  low <- resolve_threshold("M", x, control_x = stats::rnorm(5000, 8, 0.2),
                           control_kind = "fmo_q995")
  expect_identical(low$source, "fmo_q995_valley_rejected")
})

test_that("an FMO-derived threshold reports NA uncertainty with an FMO reason", {
  # The control's events are not passed to the uncertainty code, so neither
  # component is computable. A small number here would read as confidence.
  u <- threshold_uncertainty(stats::rnorm(4000), source = "fmo_q995")
  expect_true(is.na(u$u_combined))
  expect_match(u$basis, "FMO")
})

test_that("agreement is reported in units of the threshold's own uncertainty", {
  # THE POINT OF THE FEATURE. A gap of 0.3 units means nothing until you know
  # whether the cut moves by 0.05 or by 0.5 under resampling.
  thr <- data.frame(sample_id = c("S1", "S2", "S3"), panel = "P1",
                    marker = "CCR7", threshold = c(2.0, 3.5, 2.05),
                    source = "valley", stringsAsFactors = FALSE)
  fmo <- data.frame(sample_id = c("S1", "S2", "S3"), marker = "CCR7",
                    fmo_threshold = c(1.95, 2.0, 2.0),
                    fmo_sample = "FMO_CCR7", stringsAsFactors = FALSE)
  unc <- data.frame(sample_id = c("S1", "S2", "S3"), marker = "CCR7",
                    u_combined = c(0.10, 0.10, 0.50), stringsAsFactors = FALSE)

  out <- fmo_agreement(thr, fmo, unc)
  rownames(out) <- out$sample_id

  # S1: half an uncertainty away. Corroborated.
  expect_identical(out["S1", "verdict"], "corroborated by the FMO")
  # S2: fifteen uncertainties above. The derived cut is discarding signal.
  expect_gt(out["S2", "distance_in_u"], 3)
  expect_match(out["S2", "verdict"], "above the FMO")
  # S3: the same 0.05 gap as S1 but a cut five times less certain, so it is
  # comfortably inside what that threshold can explain.
  expect_equal(out["S3", "distance"], 0.05)
  expect_identical(out["S3", "verdict"], "corroborated by the FMO")

  # Worst disagreement sorts first, so the reader meets it before the rest.
  expect_identical(out$sample_id[1], "S2")
})

test_that("agreement without an uncertainty table says so rather than guessing", {
  thr <- data.frame(sample_id = "S1", panel = "P1", marker = "CCR7",
                    threshold = 2.0, source = "valley", stringsAsFactors = FALSE)
  fmo <- data.frame(sample_id = "S1", marker = "CCR7", fmo_threshold = 1.5,
                    fmo_sample = "F", stringsAsFactors = FALSE)
  out <- fmo_agreement(thr, fmo, unc = NULL)
  expect_identical(out$verdict, "no uncertainty available")
  expect_equal(out$distance, 0.5)
  expect_true(is.na(out$distance_in_u))

  expect_null(fmo_agreement(NULL, fmo))
  expect_null(fmo_agreement(thr, NULL))
})
