# =============================================================================
# --list-channels
#
# The flag answers one question: what do I write in the config. On a spectral
# panel the answer is not what the instrument recorded, because $PnS carries the
# detector and the pulse statistic alongside the marker. Naming the raw string
# produces a full set of tables containing zeros, so the resolution shown here
# has to be the same one the run applies.
# =============================================================================

#' Write a small FCS file with a chosen parameter block.
#' @noRd
fcs_with <- function(path, nm, desc, n = 40L) {
  set.seed(3)
  mat <- matrix(stats::runif(n * length(nm), 1, 1000), nrow = n,
                dimnames = list(NULL, nm))
  ff <- flowCore::flowFrame(mat)
  pd <- flowCore::pData(flowCore::parameters(ff))
  pd$desc <- desc
  flowCore::pData(flowCore::parameters(ff)) <- pd
  suppressWarnings(flowCore::write.FCS(ff, path))
  path
}

test_that("a spectral $PnS resolves to the short marker symbol", {
  dir <- withr::local_tempdir()
  f <- fcs_with(
    file.path(dir, "a.fcs"),
    nm   = c("FSC-A", "SSC-A", "SparkUV-387-A", "BUV496-A", "Time"),
    desc = c("FSC - Area", "SSC - Area",
             "CD45 : SparkUV-387 - Area", "CD19 : BUV496 - Area", NA))

  res <- capture.output(out <- cyRAVEN:::list_channels(f))

  # The detector and pulse statistic are stripped; the symbol is what a
  # population specification has to name.
  expect_true("CD45" %in% out$symbol)
  expect_true("CD19" %in% out$symbol)
  expect_false(any(grepl(" : ", out$symbol)))

  # A parameter with no description falls back to the channel name rather than
  # resolving to an empty string.
  expect_true(any(grepl("^Time", out$symbol)))
})

test_that("roles separate markers from scatter and time", {
  dir <- withr::local_tempdir()
  f <- fcs_with(
    file.path(dir, "a.fcs"),
    nm   = c("FSC-A", "FSC-H", "SSC-A", "SparkUV-387-A", "Time"),
    desc = c("FSC - Area", "FSC - Height", "SSC - Area",
             "CD45 : SparkUV-387 - Area", NA))

  invisible(capture.output(out <- cyRAVEN:::list_channels(f)))

  # Keyed on $PnN, because a scatter channel resolves its SYMBOL from $PnS and
  # so is called "FSC - Area" rather than "FSC-A". That is the same name the
  # run uses and the same one --check reports.
  role <- setNames(out$role, out$name)

  expect_match(role[["SparkUV-387-A"]], "^marker")
  expect_match(role[["FSC-A"]], "scatter")
  expect_match(role[["SSC-A"]], "scatter")
  expect_match(role[["Time"]], "acquisition time")

  # Height is not an area channel, so it carries no marker.
  expect_match(role[["FSC-H"]], "scatter")
  expect_equal(out$symbol[out$name == "SparkUV-387-A"], "CD45")
})

test_that("an identical panel across files is reported as identical", {
  dir <- withr::local_tempdir()
  nm   <- c("FSC-A", "SSC-A", "SparkUV-387-A", "BUV496-A")
  desc <- c("FSC - Area", "SSC - Area",
            "CD45 : SparkUV-387 - Area", "CD19 : BUV496 - Area")
  f1 <- fcs_with(file.path(dir, "a.fcs"), nm, desc)
  f2 <- fcs_with(file.path(dir, "b.fcs"), nm, desc)

  msg <- capture.output(invisible(cyRAVEN:::list_channels(c(f1, f2))), type = "message")
  expect_true(any(grepl("every file carries the same", msg)))
  expect_false(any(grepl("PANEL DIFFERS", msg)))
})

test_that("a panel that differs between files is named, with the missing marker", {
  dir <- withr::local_tempdir()
  f1 <- fcs_with(file.path(dir, "a.fcs"),
                 nm   = c("FSC-A", "SparkUV-387-A", "BUV496-A"),
                 desc = c("FSC - Area", "CD45 : SparkUV-387 - Area",
                          "CD19 : BUV496 - Area"))
  # Second file is missing CD19 and carries CD3 instead.
  f2 <- fcs_with(file.path(dir, "b.fcs"),
                 nm   = c("FSC-A", "SparkUV-387-A", "BV711-A"),
                 desc = c("FSC - Area", "CD45 : SparkUV-387 - Area",
                          "CD3 : BV711 - Area"))

  msg <- capture.output(invisible(cyRAVEN:::list_channels(c(f1, f2))), type = "message")
  expect_true(any(grepl("PANEL DIFFERS", msg)))
  expect_true(any(grepl("b\\.fcs", msg)))
  expect_true(any(grepl("missing: CD19", msg)))
  expect_true(any(grepl("extra: CD3", msg)))
})

test_that("an unreadable file does not stop the listing", {
  dir <- withr::local_tempdir()
  f1 <- fcs_with(file.path(dir, "a.fcs"),
                 nm = c("FSC-A", "SparkUV-387-A"),
                 desc = c("FSC - Area", "CD45 : SparkUV-387 - Area"))
  bad <- file.path(dir, "broken.fcs")
  writeLines("not an FCS file", bad)

  msg <- capture.output(invisible(cyRAVEN:::list_channels(c(f1, bad))), type = "message")
  expect_true(any(grepl("CD45", msg)))
  expect_true(any(grepl("unreadable header", msg)))
})
