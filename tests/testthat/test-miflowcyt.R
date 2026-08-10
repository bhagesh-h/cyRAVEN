# =============================================================================
# MIFlowCyt report
#
# The rule this file enforces throughout: the report may state what the files
# and the run establish, and must mark everything else as outstanding. A blank
# that reads as "not applicable" is worse than no report, because a reviewer
# cannot tell it from a completed section.
# =============================================================================

#' A keyword block resembling what an instrument writes.
fake_keywords <- function(par = 3L, spillover = FALSE, drop = character(0)) {
  kw <- list(`$PAR` = as.character(par), `$TOT` = "50000",
             `$CYT` = "FACSCanto II", `$CYTSN` = "V96500123",
             `$DATE` = "12-Mar-2024", `$BTIM` = "10:14:02", `$ETIM` = "10:19:41",
             `$OP` = "operator1", `$SYS` = "BD FACSDiva 8.0.1",
             `$P1N` = "FSC-A", `$P1S` = "", `$P1V` = "310", `$P1R` = "262144",
             `$P1E` = "0,0", `$P1B` = "18",
             `$P2N` = "FITC-A", `$P2S` = "CD3", `$P2V` = "425", `$P2R` = "262144",
             `$P2E` = "0,0", `$P2B` = "18",
             `$P3N` = "APC-A", `$P3S` = "CD4", `$P3V` = "500", `$P3R` = "262144",
             `$P3E` = "4,0", `$P3B` = "18")
  if (spillover) kw$`$SPILLOVER` <- "2,FITC-A,APC-A,1,0.02,0.01,1"
  kw[drop] <- NULL
  kw
}

fake_reads <- function(n = 2L, ...) {
  stats::setNames(lapply(seq_len(n), function(i) list(
    file = paste0("/data/donor", i, ".fcs"),
    sample_id = paste0("S", i),
    n_events = 50000L,
    keywords = fake_keywords(...))), paste0("S", seq_len(n)))
}

test_that("a keyword is found under any of its spellings, and absence is NA", {
  kw <- list(`$CYT` = "Aurora", CYT = "ignored", `$OP` = "  ")
  expect_identical(kw_get(kw, "$CYT", "CYT"), "Aurora")
  expect_identical(kw_get(kw, "CYT", "$CYT"), "ignored")
  # Whitespace-only is empty, not a value.
  expect_true(is.na(kw_get(kw, "$OP")))
  expect_true(is.na(kw_get(kw, "$NOPE")))
})

test_that("a missing keyword is named rather than silently dropped", {
  # The absence of a keyword is a fact about the acquisition and is often the
  # answer to why something could not be computed.
  expect_identical(or_absent(NA), "not recorded in the FCS file")
  expect_identical(or_absent(""), "not recorded in the FCS file")
  expect_identical(or_absent(NULL, "none"), "none")
  expect_identical(or_absent("FACSCanto II"), "FACSCanto II")
})

test_that("the detector table reports detector and marker separately", {
  p <- miflowcyt_parameters(fake_keywords())
  expect_identical(nrow(p), 3L)
  expect_identical(p$detector, c("FSC-A", "FITC-A", "APC-A"))
  # $PnS is optional. Where it is empty the row must say so rather than repeat
  # the detector name, because "which channel has no marker assigned" is exactly
  # what a reader of this table is checking.
  expect_identical(p$marker[1], "none ($PnS absent)")
  expect_identical(p$marker[2:3], c("CD3", "CD4"))
  expect_identical(p$voltage, c("310", "425", "500"))
  # $PnE is "f1,f2" with f1 = 0 meaning linear.
  expect_match(p$amplification[1], "^linear")
  expect_match(p$amplification[3], "^logarithmic")
})

test_that("an unreadable parameter block returns NULL rather than a bad table", {
  expect_null(miflowcyt_parameters(fake_keywords(drop = "$PAR")))
  expect_null(miflowcyt_parameters(list()))
})

test_that("markdown tables escape the pipe character", {
  d <- data.frame(a = "x|y", b = 1, stringsAsFactors = FALSE)
  out <- md_table(d)
  expect_length(out, 3L)
  expect_match(out[3], "x\\\\|y")
  expect_identical(md_table(d[0, ]), character(0))
})

test_that("the report completes what it can and marks the rest outstanding", {
  f <- tempfile(fileext = ".md")
  spec <- list(`CD4 T cells` = list(CD3 = "above", CD4 = "above", CD8 = "below"))
  suppressMessages(
    write_miflowcyt(f, fake_reads(2L), opt = list(transform = "arcsinh"),
                    spec = spec))
  txt <- paste(readLines(f), collapse = "\n")

  # All four MIFlowCyt parts are present as sections.
  for (h in c("## 1. Experiment overview", "## 2. Flow sample and specimen",
              "## 3. Instrumentation", "## 4. Data analysis"))
    expect_match(txt, h, fixed = TRUE)

  # Intent and biology cannot come from the files and must be flagged.
  expect_match(txt, "TO BE COMPLETED")
  # Instrument facts that CAN come from the files must actually be there.
  expect_match(txt, "FACSCanto II", fixed = TRUE)
  expect_match(txt, "V96500123", fixed = TRUE)
  expect_match(txt, "operator1", fixed = TRUE)
  # And the gating definition, which is the section this package owns outright.
  expect_match(txt, "CD4 T cells", fixed = TRUE)
  expect_match(txt, "CD8: below", fixed = TRUE)
})

test_that("compensation is reported as applied, absent or mixed, never assumed", {
  spec <- list(P = list(CD3 = "above"))

  f1 <- tempfile(fileext = ".md")
  suppressMessages(write_miflowcyt(f1, fake_reads(2L, spillover = TRUE), spec = spec))
  expect_match(paste(readLines(f1), collapse = " "), "was applied to every file")

  f2 <- tempfile(fileext = ".md")
  suppressMessages(write_miflowcyt(f2, fake_reads(2L, spillover = FALSE), spec = spec))
  expect_match(paste(readLines(f2), collapse = " "), "No `\\$SPILLOVER` keyword")

  # A cohort where only some files carry a matrix is a real and dangerous state,
  # and saying "applied" or "absent" would both be false.
  mixed <- c(fake_reads(1L, spillover = TRUE),
             stats::setNames(fake_reads(1L, spillover = FALSE), "S2"))
  f3 <- tempfile(fileext = ".md")
  suppressMessages(write_miflowcyt(f3, mixed, spec = spec))
  expect_match(paste(readLines(f3), collapse = " "), "Mixed: 1 of 2")
})

test_that("a threshold source tally is reported when the table is supplied", {
  f <- tempfile(fileext = ".md")
  thr <- data.frame(sample_id = rep(c("S1", "S2"), each = 2),
                    marker = c("CD3", "CD4"),
                    source = c("valley", "valley", "valley", "quantile_fallback"),
                    stringsAsFactors = FALSE)
  suppressMessages(write_miflowcyt(f, fake_reads(2L), spec = list(P = list(CD3 = "above")),
                                   thresholds = thr))
  txt <- paste(readLines(f), collapse = "\n")
  expect_match(txt, "`valley`: 3", fixed = TRUE)
  expect_match(txt, "`quantile_fallback`: 1", fixed = TRUE)
})

test_that("no reads means no report rather than an empty one", {
  expect_null(write_miflowcyt(tempfile(), NULL))
  expect_null(write_miflowcyt(tempfile(), list()))
})
