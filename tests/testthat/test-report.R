# =============================================================================
# Run report
#
# The claim: the outputs are presented in the order the documentation says they
# must be read, a run that excluded samples says so before anything else, and a
# section with nothing behind it is omitted rather than shown empty.
# =============================================================================

#' A results directory holding the named files. Deliberately NOT
#' withr::local_tempdir(): that is scoped to the frame it is called in, so the
#' directory would be removed the moment this helper returned.
fake_results <- function(files = character(0)) {
  d <- file.path(tempdir(), paste0("rep-", as.integer(stats::runif(1, 1, 1e9))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (f in files) {
    if (grepl("\\.csv$", f))
      utils::write.csv(data.frame(a = 1:3, b = c("x", "y", "z")),
                       file.path(d, f), row.names = FALSE)
    else writeLines("not really a png", file.path(d, f))
  }
  d
}

test_that("sections appear in the documented reading order", {
  # The order is the product. A directory listing sorts alphabetically and so
  # presents a failed QC and a headline p-value as equals.
  d <- fake_results(c("gating_qc.png", "staining_qc.csv",
                      "population_frequencies.csv", "group_comparison_stats.csv"))
  write_run_report(d)
  txt <- paste(readLines(file.path(d, "report.html")), collapse = "\n")

  pos <- function(s) regexpr(s, txt, fixed = TRUE)
  expect_gt(pos("Did the gates land"), 0)
  expect_lt(pos("Did the gates land"), pos("Is the staining usable"))
  expect_lt(pos("Is the staining usable"), pos("What are the populations"))
  expect_lt(pos("What are the populations"), pos("Do the groups differ"))
})

test_that("a section with no files behind it is omitted", {
  d <- fake_results("staining_qc.csv")
  write_run_report(d)
  txt <- paste(readLines(file.path(d, "report.html")), collapse = "\n")
  expect_match(txt, "Is the staining usable", fixed = TRUE)
  # Nothing was written for these, so claiming them would be a false statement
  # about what the run produced.
  expect_false(grepl("Do the groups differ", txt, fixed = TRUE))
  expect_false(grepl("Was the acquisition stable", txt, fixed = TRUE))
})

test_that("excluded samples are announced before any result", {
  d <- fake_results(c("staining_qc.csv", "group_comparison_stats.csv"))
  verdicts <- list(
    S1 = list(qc_status = "pass",   is_control = FALSE),
    S2 = list(qc_status = "failed", is_control = FALSE),
    S3 = list(qc_status = "pass",   is_control = TRUE))
  write_run_report(d, verdicts = verdicts)
  txt <- paste(readLines(file.path(d, "report.html")), collapse = "\n")

  expect_match(txt, "1 of 3 sample\\(s\\) failed staining QC")
  # And it comes before the first result section, not after it.
  expect_lt(regexpr("failed staining QC", txt), regexpr("Do the groups differ", txt))
})

test_that("a clean run says so rather than staying silent", {
  d <- fake_results("staining_qc.csv")
  write_run_report(d, verdicts = list(S1 = list(qc_status = "pass", is_control = FALSE),
                                      S2 = list(qc_status = "pass", is_control = FALSE)))
  txt <- paste(readLines(file.path(d, "report.html")), collapse = "\n")
  expect_match(txt, "All 2 declared sample\\(s\\) passed staining QC")
})

test_that("figures are referenced and tables are linked", {
  d <- fake_results(c("gating_qc.png", "staining_qc.csv"))
  write_run_report(d)
  txt <- paste(readLines(file.path(d, "report.html")), collapse = "\n")
  expect_match(txt, "<img src='gating_qc.png'", fixed = TRUE)
  expect_match(txt, "href='staining_qc.csv'", fixed = TRUE)
  # Values from the CSV are shown, so the reader does not have to open it to see
  # whether it is worth opening.
  expect_match(txt, "<td>x</td>", fixed = TRUE)
})

test_that("a long table is truncated with the count stated", {
  d <- withr::local_tempdir()
  utils::write.csv(data.frame(a = 1:100), file.path(d, "population_frequencies.csv"),
                   row.names = FALSE)
  write_run_report(d)
  txt <- paste(readLines(file.path(d, "report.html")), collapse = "\n")
  expect_match(txt, "of 100 rows")
})

test_that("HTML special characters in the data are escaped", {
  # A marker named with an angle bracket must not close a tag.
  d <- withr::local_tempdir()
  utils::write.csv(data.frame(marker = "<script>", verdict = "a & b"),
                   file.path(d, "staining_qc.csv"), row.names = FALSE)
  write_run_report(d)
  txt <- paste(readLines(file.path(d, "report.html")), collapse = "\n")
  expect_match(txt, "&lt;script&gt;", fixed = TRUE)
  expect_false(grepl("<td><script></td>", txt, fixed = TRUE))
  expect_match(txt, "a &amp; b", fixed = TRUE)
})

test_that("an empty or missing directory produces nothing rather than failing", {
  expect_null(write_run_report(file.path(tempdir(), "does-not-exist-at-all")))
  d <- withr::local_tempdir()
  # No outputs at all: the header still writes, and no section claims anything.
  p <- write_run_report(d)
  expect_true(file.exists(p))
  txt <- paste(readLines(p), collapse = "\n")
  expect_false(grepl("Did the gates land", txt, fixed = TRUE))
})
