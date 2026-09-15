# --total-counts: parsing the blocked layout, and the joins that must not be
# allowed to go wrong quietly.

wide_sheet <- function(path) {
  # The layout as it comes off Excel: a block-label row, a header row, spacer
  # columns between blocks, "-" for not measured, and a blank cell for the same.
  writeLines(c(
    "d0,,,,,d3,,,,,d7,",
    "Patient,PBMC count x 10^6,,,,Patient,PBMC count x 10^6,,,,Patient,PBMC count x 10^6",
    "P01,15,,,,P01,-,,,,P01,-",
    "P02,9,,,,P02,,,,,P02,26",
    "P03,25,,,,P03,15,,,,P03,11"), path)
  path
}

test_that("the blocked layout parses into one row per measured draw", {
  f <- wide_sheet(tempfile(fileext = ".csv"))
  tc <- parse_total_counts_csv(f)

  # 5 measured values: P01 d0, P02 d0, P02 d7, P03 d0, P03 d3, P03 d7 = 6
  expect_equal(nrow(tc), 6L)
  expect_setequal(names(tc), c("patient_raw", "timepoint", "total_cells"))
  expect_setequal(unique(tc$timepoint), c("d0", "d3", "d7"))

  # "x 10^6" in the header is applied, so 15 becomes fifteen million.
  expect_equal(tc$total_cells[tc$patient_raw == "P01" & tc$timepoint == "d0"], 15e6)
  expect_equal(attr(tc, "scale_label"), "10^6")
})

test_that("'-' and blank are dropped rather than read as zero", {
  f <- wide_sheet(tempfile(fileext = ".csv"))
  tc <- parse_total_counts_csv(f)
  # P01 has a "-" at d3 and d7; P02 has a blank at d3.
  expect_false(any(tc$patient_raw == "P01" & tc$timepoint %in% c("d3", "d7")))
  expect_false(any(tc$patient_raw == "P02" & tc$timepoint == "d3"))
  expect_true(all(tc$total_cells > 0))
})

test_that("a missing multiplier is taken at face value, not guessed", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Patient,PBMC count", "P01,1500000", "P02,900000"), f)
  tc <- parse_total_counts_csv(f)
  expect_equal(tc$total_cells, c(1500000, 900000))
  expect_true(is.na(attr(tc, "scale_label")))
  # No timepoint row above the header, so the join key is the patient alone.
  expect_true(all(is.na(tc$timepoint)))
})

test_that("the join keeps a patient's timepoints apart", {
  f <- wide_sheet(tempfile(fileext = ".csv"))
  tc <- parse_total_counts_csv(f)
  smap <- data.frame(
    sample_id  = c("P03_d0", "P03_d3", "P03_d7"),
    patient_id = c("P03", "P03", "P03"),
    timepoint  = c("d0", "d3", "d7"),
    stringsAsFactors = FALSE)
  m <- match_total_counts_samples(tc, smap)

  expect_equal(nrow(m), 3L)
  expect_equal(m$total_cells[m$sample_id == "P03_d0"], 25e6)
  expect_equal(m$total_cells[m$sample_id == "P03_d3"], 15e6)
  expect_equal(m$total_cells[m$sample_id == "P03_d7"], 11e6)
})

test_that("an unlabelled yield facing several draws is unmatched, not broadcast", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Patient,PBMC count x 10^6", "P03,25"), f)
  tc <- parse_total_counts_csv(f)
  smap <- data.frame(sample_id = c("P03_d0", "P03_d3"),
                     patient_id = c("P03", "P03"),
                     timepoint = c("d0", "d3"), stringsAsFactors = FALSE)
  # Two acquisitions, one undated yield: guessing which draw it belongs to
  # would put a fabricated number into every absolute figure.
  expect_equal(nrow(match_total_counts_samples(tc, smap)), 0L)
})

test_that("two yields for one acquisition is an error, not a silent pick", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Patient,PBMC count x 10^6", "P01,15", "P01,22"), f)
  tc <- parse_total_counts_csv(f)
  smap <- data.frame(sample_id = "P01_d0", patient_id = "P01",
                     timepoint = "d0", stringsAsFactors = FALSE)
  expect_error(match_total_counts_samples(tc, smap), "more than one total")
})

test_that("a sheet with no header row names the fix", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("P01,15", "P02,9"), f)
  expect_error(parse_total_counts_csv(f), "no header row")
})

test_that("attach_absolute_cells adds columns without touching the share", {
  ab <- data.frame(
    sample_id = c("S1", "S1", "S2", "S2", "S3"),
    population = c("1", "2", "1", "2", "1"),
    pct_of_gated = c(25, 75, 10, 90, 50),
    stringsAsFactors = FALSE)
  tc <- data.frame(sample_id = c("S1", "S2"), total_cells = c(1e6, 2e6),
                   stringsAsFactors = FALSE)
  out <- attach_absolute_cells(ab, tc)

  expect_identical(out$pct_of_gated, ab$pct_of_gated)   # never overwritten
  expect_equal(out$cells_absolute[1], 0.25 * 1e6)
  expect_equal(out$cells_absolute[3], 0.10 * 2e6)
  # S3 has no external total: NA, and still present, so the two tables line up.
  expect_true(is.na(out$cells_absolute[5]))
  expect_equal(nrow(out), nrow(ab))
  expect_true(all(out$count_basis[1:4] == "dual-platform: share x external total"))
})

test_that("the multiplier is read from several spellings", {
  expect_equal(total_counts_scale("PBMC count x 10^6")$scale, 1e6)
  expect_equal(total_counts_scale("cells x10^9/L")$scale, 1e9)
  expect_equal(total_counts_scale("count (millions)")$scale, 1e6)
  expect_equal(total_counts_scale("PBMC count")$scale, 1)
  expect_true(is.na(total_counts_scale("PBMC count")$label))
})
