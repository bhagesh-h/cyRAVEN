# The failure report
# =============================================================================
# A run that stops still has to leave a record that explains itself, because the
# person opening the results directory afterwards may never have seen stderr.

test_that("known failure modes are interpreted rather than echoed", {
  d <- diagnose_failure("these input files are not in the sample sheet: a.fcs")
  expect_false(is.null(d))
  expect_match(d$meaning, "no row")
  expect_match(d$action, "--write-samples")

  d2 <- diagnose_failure("cannot allocate vector of size 4.1 Gb")
  expect_match(d2$action, "--max-events-per-file")

  d3 <- diagnose_failure("the sample sheet gives conflicting values for the same subject")
  expect_match(d3$action, "--check")
})

test_that("an unrecognised message yields no invented interpretation", {
  expect_null(diagnose_failure("something nobody has seen before"))
})

test_that("every catalogue pattern is a valid regex and has both fields", {
  for (e in failure_catalogue()) {
    expect_silent(grepl(e$pattern, "probe", perl = TRUE))
    expect_true(nzchar(e$meaning))
    expect_true(nzchar(e$action))
  }
})

test_that("a failed run writes a report carrying the diagnosis", {
  d <- file.path(tempdir(), paste0("failrep", as.integer(runif(1, 1, 1e6))))
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  # An output written before the failure, which must survive into the report.
  write.csv(data.frame(sample_id = "S1", qc_status = "pass"),
            file.path(d, "staining_qc.csv"), row.names = FALSE)

  p <- suppressMessages(write_run_report(
    d, failure = simpleError("these input files are not in the sample sheet: x.fcs")))
  expect_true(file.exists(p))
  h <- paste(readLines(p, warn = FALSE), collapse = "\n")

  expect_match(h, "did not complete")
  expect_match(h, "Why this run failed")
  expect_match(h, "not in the sample sheet")   # the message itself
  expect_match(h, "no row")                    # the interpretation
  expect_match(h, "--write-samples")           # the action
  # The partial output is embedded, not discarded.
  expect_match(h, "staining_qc.csv")
})

test_that("a failed run with no outputs at all still gets a report", {
  d <- file.path(tempdir(), paste0("failempty", as.integer(runif(1, 1, 1e6))))
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- suppressMessages(write_run_report(d, failure = simpleError("boom")))
  expect_true(file.exists(p))
  expect_match(paste(readLines(p, warn = FALSE), collapse = "\n"), "boom")
})

test_that("the report records the stage the run reached", {
  d <- file.path(tempdir(), paste0("failstage", as.integer(runif(1, 1, 1e6))))
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  cyRAVEN:::log_reset()
  suppressMessages(cyRAVEN:::log_step("STEP 4 - scoring populations"))
  p <- suppressMessages(write_run_report(d, failure = simpleError("boom")))
  expect_match(paste(readLines(p, warn = FALSE), collapse = "\n"),
               "STEP 4 - scoring populations")
})

test_that("run_cyraven re-raises the original error after reporting", {
  d <- file.path(tempdir(), paste0("failraise", as.integer(runif(1, 1, 1e6))))
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  # A directory with no FCS files fails early, which is enough to exercise the
  # path. The caller must still see the real error, not one about reporting.
  expect_error(suppressMessages(
    run_cyraven(list(dir = file.path(d, "no-such-dir"), outdir = d))))
})
