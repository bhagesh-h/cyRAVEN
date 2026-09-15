# =============================================================================
# Run report
#
# The claims: the outputs are presented in the order the documentation says they
# must be read; a run that excluded samples says so before anything else; a
# section with nothing behind it is omitted rather than shown empty; the file is
# self-contained, so it references nothing; and no output the run wrote is left
# out of it.
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

read_report <- function(d) paste(readLines(file.path(d, "report.html"),
                                           warn = FALSE), collapse = "\n")

#' A section heading as it appears in the markup.
#'
#' Headings used to carry an ordinal ("1. Gate placement") and the assertions
#' below matched that string. Without the number a bare title is an ordinary
#' phrase that also occurs in prose -- "Provenance" and "Population abundance"
#' both do -- so match the element instead. This is what the reader sees as a
#' heading, and nothing else in the document produces it.
sec_title <- function(t) paste0("sec-title'>", t, "<")

test_that("sections appear in the documented reading order", {
  # The order is the product. A directory listing sorts alphabetically and so
  # presents a failed QC and a headline p-value as equals.
  d <- fake_results(c("gating_qc.png", "staining_qc.csv",
                      "population_frequencies.csv", "group_comparison_stats.csv"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)

  pos <- function(s) regexpr(sec_title(s), txt, fixed = TRUE)
  expect_gt(pos("Gate placement"), 0)
  expect_lt(pos("Gate placement"), pos("Staining quality control"))
  expect_lt(pos("Staining quality control"),
            pos("Population abundance and the shared embedding"))
  expect_lt(pos("Population abundance and the shared embedding"),
            pos("Between-group differences"))
})

test_that("headings state what the section reports rather than asking", {
  d <- fake_results(c("gating_qc.png", "staining_qc.csv"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, sec_title("Gate placement"), fixed = TRUE)
  # The old question form must not survive anywhere in the headings.
  expect_false(grepl("Did the gates land", txt, fixed = TRUE))
  expect_false(grepl("Is the staining usable", txt, fixed = TRUE))
})

test_that("a section with no files behind it is omitted", {
  d <- fake_results("staining_qc.csv")
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, sec_title("Staining quality control"), fixed = TRUE)
  # Nothing was written for these, so claiming them would be a false statement
  # about what the run produced.
  expect_false(grepl(sec_title("Between-group differences"), txt, fixed = TRUE))
  expect_false(grepl(sec_title("Acquisition stability"), txt, fixed = TRUE))
})

test_that("excluded samples are announced before any result", {
  d <- fake_results(c("staining_qc.csv", "group_comparison_stats.csv"))
  verdicts <- list(
    S1 = list(qc_status = "pass",   is_control = FALSE),
    S2 = list(qc_status = "failed", is_control = FALSE),
    S3 = list(qc_status = "pass",   is_control = TRUE))
  suppressMessages(write_run_report(d, verdicts = verdicts))
  txt <- read_report(d)

  expect_match(txt, "1 of 3 sample\\(s\\) failed staining QC")
  # And it comes before the first result section, not after it. Anchored on the
  # section element rather than its title, because the title also appears in the
  # sidebar, which precedes everything.
  expect_lt(regexpr("failed staining QC", txt, fixed = TRUE),
            regexpr("<section class='sec' id='s8'", txt, fixed = TRUE))
})

test_that("a clean run says so rather than staying silent", {
  d <- fake_results("staining_qc.csv")
  suppressMessages(write_run_report(
    d, verdicts = list(S1 = list(qc_status = "pass", is_control = FALSE),
                       S2 = list(qc_status = "pass", is_control = FALSE))))
  expect_match(read_report(d), "All 2 declared sample\\(s\\) passed staining QC")
})

test_that("the report is self-contained: figures embed, nothing is linked", {
  d <- fake_results(c("gating_qc.png", "staining_qc.csv"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  # The figure travels inside the file as a data URI.
  expect_match(txt, "src='data:image/png;base64,", fixed = TRUE)
  # A report that points at a neighbouring file breaks the moment it is moved,
  # so no src or href may name one.
  expect_false(grepl("src='gating_qc.png'", txt, fixed = TRUE))
  expect_false(grepl("href='staining_qc.csv'", txt, fixed = TRUE))
  # No external resource of any kind: no CDN, no webfont, no remote script.
  expect_false(grepl("src=\"http", txt, fixed = TRUE))
  expect_false(grepl("src='http", txt, fixed = TRUE))
  expect_false(grepl("@import", txt, fixed = TRUE))
})

test_that("table data travels with the report and is not truncated", {
  d <- withr::local_tempdir()
  utils::write.csv(data.frame(a = 1:100), file.path(d, "population_frequencies.csv"),
                   row.names = FALSE)
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, "100 rows", fixed = TRUE)
  # Every row is carried, because the browser does the paging: the last row must
  # be present even though only 50 are shown by default.
  expect_match(txt, "\"100\"", fixed = TRUE)
  # And the controls that page it exist.
  expect_match(txt, "All rows", fixed = TRUE)
  expect_match(txt, "Export CSV", fixed = TRUE)
  expect_match(txt, "Search this table", fixed = TRUE)
})

test_that("HTML special characters in the data cannot break out", {
  # A marker named with an angle bracket must not close a tag, in the JSON
  # payload any more than in markup.
  d <- withr::local_tempdir()
  utils::write.csv(data.frame(marker = "</script><b>", verdict = "a & b"),
                   file.path(d, "staining_qc.csv"), row.names = FALSE)
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  # The closing tag is neutralised inside the JSON block: it is carried as
  # <\/script>, which JSON reads back as the original text but which no HTML
  # parser sees as a closing tag.
  expect_false(grepl("</script><b>", txt, fixed = TRUE))
  expect_match(txt, "<\\/script>", fixed = TRUE)
})

test_that("a quote in the data does not corrupt the JSON payload", {
  d <- withr::local_tempdir()
  utils::write.csv(data.frame(note = 'he said "no"', n = 1),
                   file.path(d, "staining_qc.csv"), row.names = FALSE)
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, '\\\\"no\\\\"')
})

test_that("no output the run wrote is left out of the report", {
  # The sections name their files explicitly, so a file no section names would
  # be silently absent. The sweep is what stops that.
  d <- fake_results(c("staining_qc.csv", "some_future_output.csv",
                      "some_future_figure.png"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, "some_future_output.csv", fixed = TRUE)
  expect_match(txt, "some_future_figure.png", fixed = TRUE)
  expect_match(txt, sec_title("Further outputs"), fixed = TRUE)
})

test_that("the sweep section is absent when every output is placed", {
  d <- fake_results(c("gating_qc.png", "staining_qc.csv"))
  suppressMessages(write_run_report(d))
  expect_false(grepl(sec_title("Further outputs"), read_report(d), fixed = TRUE))
})

test_that("the sidebar indexes every section, figure and table", {
  d <- fake_results(c("gating_qc.png", "staining_qc.csv",
                      "group_comparison_stats.csv"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, "nav class='side'", fixed = TRUE)
  expect_match(txt, "href='#fig-gating-qc'", fixed = TRUE)
  expect_match(txt, "href='#tab-staining-qc'", fixed = TRUE)
  expect_match(txt, "href='#s1'", fixed = TRUE)
})

test_that("figures are collapsible, zoomable and downloadable", {
  d <- fake_results("gating_qc.png")
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, "<details", fixed = TRUE)     # collapsible sections
  expect_match(txt, "cyZoom(", fixed = TRUE)      # zoom
  expect_match(txt, "Full resolution PNG", fixed = TRUE)
  expect_match(txt, "id='lb'", fixed = TRUE)      # the lightbox itself
})

test_that("an empty or missing directory produces nothing rather than failing", {
  expect_null(write_run_report(file.path(tempdir(), "does-not-exist-at-all")))
  d <- withr::local_tempdir()
  # No outputs at all and no failure: there is nothing to report.
  expect_null(suppressMessages(write_run_report(d)))
})

test_that("a figure's split-out variants follow the pooled one, ordered by number", {
  d <- withr::local_tempdir()
  file.create(file.path(d, c("umap_markers_d10.png", "umap_markers_d3.png",
                             "umap_markers.png", "umap_markers_d0.png",
                             "umap_markers_d7.png")))
  # Pooled FIRST. The collation order on the usual locales ignores punctuation
  # when comparing, so a plain sort() puts umap_markers_d0.png ahead of
  # umap_markers.png and the tab strip opens on a single day rather than on the
  # pooled view the section tells the reader to take first.
  expect_identical(pooled_then_split(d, "umap_markers"),
                   c("umap_markers.png", "umap_markers_d0.png",
                     "umap_markers_d3.png", "umap_markers_d7.png",
                     "umap_markers_d10.png"))
})

test_that("a panel suffix does not stand in for the timepoint it precedes", {
  d <- withr::local_tempdir()
  file.create(file.path(d, c("umap_markers_panel_1_d10.png",
                             "umap_markers_panel_1_d3.png",
                             "umap_markers_panel_1.png",
                             "umap_markers_panel_1_d0.png")))
  # Ordering reads the LAST digit run, so the 1 of _panel_1 is not mistaken for
  # the visit number and every panel-1 figure would not collapse to one rank.
  expect_identical(pooled_then_split(d, "umap_markers"),
                   c("umap_markers_panel_1.png", "umap_markers_panel_1_d0.png",
                     "umap_markers_panel_1_d3.png",
                     "umap_markers_panel_1_d10.png"))
})

test_that("the by-timepoint twin is named beside its pooled figure, not left to the catch-all", {
  d <- fake_results(c("functional_markers.png", "functional_markers_by_timepoint.png",
                      "population_ratios.png", "population_ratios_by_timepoint.png",
                      "absolute_vs_share.png", "absolute_vs_share_by_timepoint.png"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  # Section 8 holds them; the catch-all exists to be empty. Its whole purpose is
  # to surface files no section names, so a day-separated figure reaching it is
  # the bug this guards: present on disk, invisible where anyone would look.
  # Split on the section opener rather than slicing with a greedy regex: the
  # report is one long string with newlines and several ids sharing the "s8"
  # prefix, which makes a regex slice easy to get quietly wrong.
  chunks <- strsplit(txt, "<section class='sec' id='", fixed = TRUE)[[1]]
  s8 <- chunks[startsWith(chunks, "s8'")]
  expect_length(s8, 1L)
  for (f in c("functional_markers_by_timepoint", "population_ratios_by_timepoint",
              "absolute_vs_share_by_timepoint"))
    expect_match(s8, f, fixed = TRUE)
})

test_that("figures and tables are separate tab groups", {
  d <- fake_results(c("gating_qc.png", "staining_qc.csv"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  expect_match(txt, ">Figures <span", fixed = TRUE)
  expect_match(txt, ">Tables <span", fixed = TRUE)
  # Each group is wrapped, which is what keeps cyTab() from blanking the other
  # group's open pane: it searches the button's grandparent for panes.
  expect_gte(length(gregexpr("class='tabgroup'", txt, fixed = TRUE)[[1]]), 2L)
})

test_that("provenance is the first section", {
  d <- fake_results(c("gating_qc.png"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  chunks <- strsplit(txt, "<section class='sec' id='", fixed = TRUE)[[1]]
  expect_true(startsWith(chunks[2], "s0'"))
  expect_match(chunks[2], "Provenance", fixed = TRUE)
})

test_that("no two sections share an id", {
  # A duplicate id makes every link to the second of them land on the first, so
  # that section becomes unreachable from the contents while still rendering.
  d <- fake_results(c("gating_qc.png", "umap_overview.png", "marker_state.png",
                      "clinical_landscape.png", "timepoint_subset_balance.png",
                      "batch_diagnostic.png", "staining_qc.csv"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  ids <- regmatches(txt, gregexpr("<section class='sec' id='[^']+'", txt))[[1]]
  expect_identical(anyDuplicated(ids), 0L)
})

test_that("sections are ordered by the step that produced them", {
  # Acquisition stability is computed at STEP 1b, before the gates at STEP 2,
  # so it reads first. Constructed in the opposite order in the source.
  d <- fake_results(c("acquisition_qc.png", "gating_qc.png", "staining_qc.csv",
                      "population_frequencies.csv", "group_comparison_stats.csv"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  pos <- function(s) regexpr(sec_title(s), txt, fixed = TRUE)
  expect_lt(pos("Acquisition stability"), pos("Gate placement"))
  expect_lt(pos("Gate placement"), pos("Staining quality control"))
  expect_lt(pos("Staining quality control"),
            pos("Population abundance and the shared embedding"))
  expect_lt(pos("Population abundance and the shared embedding"),
            pos("Between-group differences"))
})

test_that("unsupervised discovery is pulled ahead of the declared analysis", {
  # It reaches the reader through the catch-all at the foot of the report
  # otherwise, which is where the one view computed WITHOUT the specification
  # was landing in a report built entirely around it.
  d <- fake_results(c("gating_qc.png", "population_frequencies.csv"))
  dir.create(file.path(d, "explore"))
  writeLines("not really a png",
             file.path(d, "explore", "explore_umap_clusters.png"))
  suppressMessages(write_run_report(d))
  txt <- read_report(d)
  pos <- function(s) regexpr(sec_title(s), txt, fixed = TRUE)
  expect_gt(pos("Unsupervised discovery"), 0)
  expect_lt(pos("Provenance"), pos("Unsupervised discovery"))
  expect_lt(pos("Unsupervised discovery"), pos("Gate placement"))
})
