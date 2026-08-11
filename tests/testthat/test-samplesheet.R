# The unified sample sheet
# =============================================================================
# The sheet is a different way to supply the same facts as the three separate
# tables, so the tests that matter are the ones asserting equivalence and the
# ones asserting the hazards the one-file shape introduces are caught.

fake_fcs <- function(n = 3L) paste0("s", seq_len(n), ".fcs")

write_tmp <- function(txt, ext = ".csv") {
  p <- tempfile(fileext = ext)
  writeLines(txt, p)
  p
}

test_that("a minimal sheet resolves the acquisition columns", {
  p <- write_tmp(c("file,sample_id,patient_id,is_control",
                   "s1.fcs,A,P1,FALSE",
                   "s2.fcs,B,P2,TRUE",
                   "s3.fcs,C,P2,FALSE"))
  r <- read_samplesheet(p, fake_fcs())
  expect_equal(nrow(r$smap), 3L)
  expect_equal(r$smap$sample_id, c("A", "B", "C"))
  # is_control is coerced to logical from any of the accepted spellings.
  expect_identical(r$smap$is_control, c(FALSE, TRUE, FALSE))
  expect_null(r$patients)
  expect_null(r$counts)
})

test_that("a file with no row is a fatal error, not a guess", {
  p <- write_tmp(c("file,sample_id", "s1.fcs,A", "s2.fcs,B"))
  expect_error(read_samplesheet(p, fake_fcs(3L)),
               "not in the sample sheet")
})

test_that("two rows naming the same file are rejected", {
  p <- write_tmp(c("file,sample_id", "s1.fcs,A", "s1.fcs,B"))
  expect_error(read_samplesheet(p, "s1.fcs"), "duplicate file entries")
})

test_that("a sheet with no file column is rejected", {
  p <- write_tmp(c("sample_id,cohort", "A,x"))
  expect_error(read_samplesheet(p, fake_fcs(1L)), "must contain a 'file' column")
})

test_that("subject columns are lifted to one row per patient", {
  p <- write_tmp(c("file,sample_id,patient_id,cohort,sex,age_years",
                   "s1.fcs,A,P1,cases,m,40",
                   "s2.fcs,B,P1,cases,m,40",   # same patient, agreeing
                   "s3.fcs,C,P2,controls,w,35"))
  r <- read_samplesheet(p, fake_fcs())
  expect_equal(nrow(r$patients), 2L)
  expect_equal(sort(r$patients$patient_id), c("P1", "P2"))
  # The shared normaliser must have run: sex translated, age numeric.
  expect_equal(sort(r$patients$sex), c("female", "male"))
  expect_type(r$patients$age_years, "double")
})

test_that("rows of one subject that disagree are a fatal error", {
  # THE hazard of a per-file sheet: a subject attribute repeats, so the copies
  # can contradict each other. Taking the first would silently pick one.
  p <- write_tmp(c("file,sample_id,patient_id,sex",
                   "s1.fcs,A,P1,m",
                   "s2.fcs,B,P1,w",
                   "s3.fcs,C,P2,m"))
  expect_error(read_samplesheet(p, fake_fcs()),
               "conflicting values for the same subject")
  # The message must name the subject and the column, or it cannot be acted on.
  err <- tryCatch(read_samplesheet(p, fake_fcs()), error = conditionMessage)
  expect_match(err, "P1")
  expect_match(err, "sex")
})

test_that("a blank in one row does not count as a disagreement", {
  p <- write_tmp(c("file,sample_id,patient_id,sex",
                   "s1.fcs,A,P1,m",
                   "s2.fcs,B,P1,",
                   "s3.fcs,C,P2,w"))
  r <- read_samplesheet(p, fake_fcs())
  expect_equal(r$patients$sex[r$patients$patient_id == "P1"], "male")
})

test_that("count columns become a long table in cells per uL", {
  p <- write_tmp(c("file,sample_id,count.Granulocytes,count.Monocytes",
                   "s1.fcs,A,3810,420",
                   "s2.fcs,B,2990,395",
                   "s3.fcs,C,,510"))
  r <- read_samplesheet(p, fake_fcs())
  expect_equal(sort(unique(r$counts$population)), c("Granulocytes", "Monocytes"))
  # The blank is dropped rather than becoming a zero count.
  expect_equal(nrow(r$counts), 5L)
  expect_equal(r$counts$cells_per_ul[r$counts$sample_id == "A" &
                                     r$counts$population == "Granulocytes"], 3810)
})

test_that("cells per mL is converted rather than assumed", {
  p <- write_tmp(c("file,sample_id,count.Lymphocytes",
                   "s1.fcs,A,1650000", "s2.fcs,B,1880000", "s3.fcs,C,1200000"))
  r <- read_samplesheet(p, fake_fcs(), count_unit = "cells/mL")
  expect_equal(r$counts$cells_per_ul[1], 1650)
})

test_that("unreserved columns are kept and reported as study variables", {
  p <- write_tmp(c("file,sample_id,batch,treatment_arm",
                   "s1.fcs,A,b1,placebo", "s2.fcs,B,b1,drug", "s3.fcs,C,b2,drug"))
  r <- read_samplesheet(p, fake_fcs())
  expect_true(all(c("batch", "treatment_arm") %in% r$study_columns))
  expect_equal(r$smap$treatment_arm, c("placebo", "drug", "drug"))
})

test_that("German headers resolve through the alias map", {
  p <- write_tmp(c("file,sample_id,Patient ID,Geschlecht,Kohorte",
                   "s1.fcs,A,P1,m,Kontrollen",
                   "s2.fcs,B,P2,w,Kontrollen",
                   "s3.fcs,C,P3,m,Patienten"))
  r <- read_samplesheet(p, fake_fcs())
  expect_equal(nrow(r$patients), 3L)
  expect_true(all(r$patients$sex %in% c("male", "female")))
})

test_that("a semicolon-separated sheet is detected", {
  p <- write_tmp(c("file;sample_id;cohort",
                   "s1.fcs;A;cases", "s2.fcs;B;cases", "s3.fcs;C;controls"))
  r <- read_samplesheet(p, fake_fcs())
  expect_equal(r$smap$sample_id, c("A", "B", "C"))
})

test_that("the template covers every input file", {
  d <- file.path(tempdir(), paste0("tmpl", as.integer(runif(1, 1, 1e6))))
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- file.path(d, "samples.csv")
  cyRAVEN:::write_samplesheet_template(fake_fcs(4L), p)
  t <- read.csv(p, stringsAsFactors = FALSE)
  expect_equal(nrow(t), 4L)
  expect_true("file" %in% names(t))
  # A template must be readable by the reader it is a template for.
  expect_silent(suppressMessages(read_samplesheet(p, fake_fcs(4L))))
})

test_that("the sheet and the three separate tables carry the same facts", {
  # Equivalence at the level this can be tested without running the pipeline:
  # the sheet's smap must match what load_sample_map() produces from the same
  # rows, and its patients must match load_patient_table()'s output.
  sheet <- write_tmp(c("file,sample_id,patient_id,is_control,cohort,sex,age_years",
                       "s1.fcs,A,P1,FALSE,cases,m,40",
                       "s2.fcs,B,P2,FALSE,controls,w,35",
                       "s3.fcs,C,P3,TRUE,controls,m,50"))
  smap_f <- write_tmp(c("file,sample_id,patient_id,is_control",
                        "s1.fcs,A,P1,FALSE",
                        "s2.fcs,B,P2,FALSE",
                        "s3.fcs,C,P3,TRUE"))
  pat_f <- write_tmp(c("patient_id,cohort,sex,age_years",
                       "P1,cases,m,40", "P2,controls,w,35", "P3,controls,m,50"))

  r <- read_samplesheet(sheet, fake_fcs())
  sm <- load_sample_map(smap_f, fake_fcs())
  pt <- suppressMessages(load_patient_table(pat_f))

  for (cn in names(sm)) expect_equal(r$smap[[cn]], sm[[cn]], info = cn)
  for (cn in names(pt))
    expect_equal(r$patients[[cn]][order(r$patients$patient_id)],
                 pt[[cn]][order(pt$patient_id)], info = cn)
})
