# Options with a fixed set of values stop the run on anything else. They were
# read with a fallback, so `--subsample Rare` ran the uniform draw while the
# manifest recorded "Rare".

test_that("an unknown value for a fixed-choice option stops the run", {
  d <- withr::local_tempdir()
  suppressMessages(write_demo_cohort(file.path(d, "demo")))
  run <- function(...) suppressMessages(run_cyraven(list(
    dir = file.path(d, "demo", "fcs"), outdir = file.path(d, "out"),
    no_session = TRUE, ...)))
  expect_error(run(subsample = "Rare"), "--subsample must be one of uniform, rare")
  expect_error(run(p_adjust_display = "bh"), "--p-adjust-display must be one of raw, BH")
  expect_error(run(transform = "asinh"), "--transform must be one of arcsinh, logicle, none")
})
