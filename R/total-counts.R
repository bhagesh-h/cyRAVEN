# SECTION 6c -- EXTERNAL TOTAL CELL COUNTS PER SAMPLE (--total-counts)
# =============================================================================
#
# WHY THIS EXISTS ALONGSIDE --absolute-counts (Section 6b) AND wbc_per_ul:
# the three are different measurements and are deliberately not merged.
#
#   --absolute-counts  one directly measured number PER POPULATION per sample,
#                      from a clinical lab's own gating. Independent evidence;
#                      never multiplied by anything this pipeline computed.
#   wbc_per_ul         one whole-blood concentration per PATIENT, from a
#                      haemogram. Keyed by patient_id, so it cannot vary
#                      between that patient's timepoints.
#   --total-counts     one TOTAL cell yield per ACQUISITION -- the number of
#                      cells the tube was loaded with (PBMC yield, viable count
#                      off a haemocytometer or an automated counter). Keyed by
#                      patient AND timepoint, because a yield is a property of
#                      the draw, not of the person.
#
# WHY IT IS WORTH HAVING AT ALL. Everything this pipeline reports about
# abundance is compositional: pct_of_gated is a share of the events acquired,
# and shares are constrained to sum to 100. That constraint means a frequency
# table cannot, even in principle, distinguish "this population expanded" from
# "everything else contracted" -- the composition is identical either way. The
# centred log-ratio in stats-compositional.R fixes the geometry of the test but
# not this; as its own header says, CLR does not recover absolute cell numbers,
# and if every population doubles the composition is unchanged and CLR sees
# nothing.
#
# A per-acquisition total is the missing scale factor. Multiplying it by a
# population share re-expresses that population as a cell number, and cell
# numbers are unconstrained -- they can all move the same way. In sepsis, where
# the characteristic finding is a global lymphopenia rather than a
# redistribution between subsets, that is the difference between seeing the
# effect and seeing nothing at all.
#
# WHAT IT COSTS. This is a DUAL-PLATFORM derivation: one number from this
# pipeline (the share) multiplied by one number from a separate instrument (the
# total). Both errors propagate into the product, and the interlaboratory
# literature puts dual-platform CVs at roughly 20-33 percent against 10-16 for
# single-platform bead counting. The derived numbers are therefore reported in
# their own columns and their own test, NEVER by overwriting pct_of_gated, and
# every table carrying them records which route produced them.
#
# EXPECTED LAYOUT (wide, blocked by timepoint -- how these come off Excel):
#
#     d0      |                   |  | d3      |          ...  <- block labels
#     Patient | PBMC count x 10^6 |  | Patient | PBMC ...      <- block headers
#     HS5_016 | 15                |  | HS5_016 | -             <- one row/patient
#
# A single unblocked block (just an identifier column and a count column, no
# timepoint row above) is also accepted, and then the join is on patient alone.
#
# "-", blanks and any other non-numeric cell all mean not measured, and are
# dropped rather than read as zero. A zero yield and an unperformed count are
# not the same thing and must not become the same number.

#' Normalize any --total-counts input down to one CSV on disk
#'
#' Mirrors `absolute_counts_to_csv()` deliberately: same reasons (one
#' irregular-layout scan, written once, into outdir so the user can open
#' exactly what was parsed), and the same CSV convention every other tabular
#' input to this pipeline already uses.
#' @param path File path.
#' @param outdir Directory to write outputs to.
#' @keywords internal
total_counts_to_csv <- function(path, outdir) {
  ext <- tolower(tools::file_ext(path))
  grid <- switch(ext,
    xlsx = ,
    xls = {
      if (!requireNamespace("readxl", quietly = TRUE))
        stop("--total-counts points at an Excel file, which needs the 'readxl' ",
             "package.\n  run:  install.packages(\"readxl\")", call. = FALSE)
      # col_types = "text": the block-label row, the header row and the numeric
      # body all share the same columns, so letting readxl infer one type per
      # column would turn two of the three into NA.
      as.data.frame(readxl::read_excel(path, col_names = FALSE, col_types = "text",
                                       .name_repair = "minimal"),
                    stringsAsFactors = FALSE)
    },
    tsv = ,
    txt = utils::read.delim(path, header = FALSE, colClasses = "character",
                            check.names = FALSE, blank.lines.skip = FALSE,
                            fill = TRUE, quote = "", stringsAsFactors = FALSE),
    NULL)
  if (is.null(grid) && ext != "csv")
    stop("--total-counts: unrecognized extension '.", ext,
         "' (expected .xlsx, .xls, .csv or .tsv)", call. = FALSE)
  csv_path <- file.path(outdir, "total_counts_raw.csv")
  if (is.null(grid)) {
    file.copy(path, csv_path, overwrite = TRUE)
  } else {
    grid[] <- lapply(grid, function(col) ifelse(is.na(col), "", trimws(as.character(col))))
    utils::write.table(grid, csv_path, sep = ",", row.names = FALSE, col.names = FALSE,
                       na = "", qmethod = "double", fileEncoding = "UTF-8")
  }
  log_msg("  --total-counts: normalized ", basename(path), " -> ", csv_path)
  csv_path
}

#' Read the multiplier out of a count column's own header text
#'
#' WHY THIS IS NEITHER OPTIONAL NOR GUESSED SILENTLY. "PBMC count x 10^6" and
#' "PBMC count" differ by a factor of a million. Getting it wrong produces a
#' table that is internally consistent, passes every downstream check, and is
#' wrong by six orders of magnitude. So the multiplier is read from the header
#' wherever one is stated, and where none is stated the values are taken at
#' face value with a NOTE naming the header it looked at, which the user can
#' confirm against a file they can open.
#' @param txt Header text of the count column.
#' @return list(scale, label)
#' @keywords internal
total_counts_scale <- function(txt) {
  # The multiplication sign is written as an escape rather than literally:
  # R CMD check requires package R code to be ASCII, and a spreadsheet header
  # reading "PBMC count x 10^6" may use either character for the x.
  t1 <- gsub("\u00d7", "x", tolower(txt %||% ""))
  exp10 <- if (grepl("10\\s*\\^?\\s*12|e12", t1)) 12L
    else if (grepl("10\\s*\\^?\\s*9|e9|billion", t1)) 9L
    else if (grepl("10\\s*\\^?\\s*6|e6|million|mio", t1)) 6L
    else if (grepl("10\\s*\\^?\\s*3|e3|thousand", t1)) 3L
    else NA_integer_
  if (is.na(exp10)) list(scale = 1, label = NA_character_)
  else list(scale = 10^exp10, label = paste0("10^", exp10))
}

#' Parse the normalized CSV into long format: patient_raw, timepoint, total_cells
#'
#' Blocks are located by their header cell rather than by a fixed column
#' offset, because the spacer columns between blocks are a formatting choice
#' that varies between two exports of the same sheet.
#' @param csv_path The csv path.
#' @keywords internal
parse_total_counts_csv <- function(csv_path) {
  grid <- as.matrix(utils::read.csv(csv_path, header = FALSE, colClasses = "character",
                                    check.names = FALSE, na.strings = character(0),
                                    blank.lines.skip = FALSE, fill = TRUE))
  if (nrow(grid) < 2L)
    stop("--total-counts: fewer than two rows; expected a header row and at ",
         "least one data row", call. = FALSE)
  # trimws() on a matrix goes through sub(), which drops dim. Restore it.
  dm <- dim(grid)
  grid[is.na(grid)] <- ""
  grid <- trimws(grid)
  dim(grid) <- dm

  # The header row is the one carrying the subject-key label. Searched over the
  # first few rows rather than assumed to be row 1, because a blocked sheet
  # puts the timepoint labels above it.
  key_pat <- "^(patient|patient.?id|subject|sample|sample.?id|id|donor)$"
  hdr_row <- NA_integer_
  for (i in seq_len(min(6L, nrow(grid)))) {
    if (any(grepl(key_pat, tolower(grid[i, ])))) { hdr_row <- i; break }
  }
  if (is.na(hdr_row))
    stop("--total-counts: no header row found. Expected a row containing a cell ",
         "reading 'Patient' (or Subject/Sample/Donor/ID) above the identifiers.",
         call. = FALSE)
  if (hdr_row == nrow(grid))
    stop("--total-counts: the header row is the last row; there are no data rows",
         call. = FALSE)

  key_cols <- which(grepl(key_pat, tolower(grid[hdr_row, ])))
  blk_row  <- if (hdr_row > 1L) hdr_row - 1L else NA_integer_
  body     <- grid[seq.int(hdr_row + 1L, nrow(grid)), , drop = FALSE]
  is_num   <- function(x) suppressWarnings(!is.na(as.numeric(x)) & nzchar(trimws(x)))

  blocks <- lapply(seq_along(key_cols), function(b) {
    kc  <- key_cols[b]
    end <- if (b < length(key_cols)) key_cols[b + 1L] - 1L else ncol(grid)
    if (kc + 1L > end) return(NULL)
    # The value column is the first column right of the key, inside this block,
    # carrying a header AND at least one number beneath it. "First with data"
    # rather than "kc + 1" so a sheet with an annotation column between the
    # identifier and the count still parses.
    vc <- NA_integer_
    for (j in seq.int(kc + 1L, end)) {
      if (nzchar(grid[hdr_row, j]) && any(is_num(body[, j]))) { vc <- j; break }
    }
    if (is.na(vc)) return(NULL)
    tp <- if (!is.na(blk_row)) {
      lab <- grid[blk_row, seq_len(kc)]
      lab <- lab[nzchar(lab)]
      if (length(lab)) lab[length(lab)] else NA_character_
    } else NA_character_
    list(key_col = kc, val_col = vc, timepoint = tp, header = grid[hdr_row, vc])
  })
  blocks <- Filter(Negate(is.null), blocks)
  if (!length(blocks))
    stop("--total-counts: found the header row but no column beneath it holds a ",
         "number. Check that a count column sits beside each identifier column.",
         call. = FALSE)

  sc <- total_counts_scale(blocks[[1]]$header)
  if (is.na(sc$label)) {
    log_msg("  NOTE --total-counts: no multiplier found in the count header ('",
            blocks[[1]]$header, "'); values are taken as they stand. If the ",
            "sheet means millions of cells, VERIFY that before quoting any ",
            "cells_absolute number.")
  } else {
    log_msg("  --total-counts: header '", blocks[[1]]$header,
            "' -> multiplying by ", sc$label)
  }

  rows <- lapply(blocks, function(bk) {
    key <- body[, bk$key_col]
    raw <- body[, bk$val_col]
    keep <- nzchar(key) & is_num(raw)
    if (!any(keep)) return(NULL)
    data.frame(patient_raw = unname(key[keep]),
               timepoint   = if (is.na(bk$timepoint)) NA_character_ else bk$timepoint,
               total_cells = unname(as.numeric(raw[keep])) * sc$scale,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out))
    stop("--total-counts: every data cell was blank, '-' or non-numeric",
         call. = FALSE)

  # A yield of zero is almost always a placeholder that survived as a number.
  # Dropped rather than propagated: it would divide the sample out of every
  # absolute figure without any error being raised anywhere.
  z <- out$total_cells <= 0
  if (any(z)) {
    log_msg("  NOTE --total-counts: ", sum(z), " row(s) carry a total of zero or ",
            "less and are dropped (a zero yield and an unperformed count are ",
            "not the same measurement)")
    out <- out[!z, , drop = FALSE]
  }
  if (!nrow(out))
    stop("--total-counts: every row was zero or negative", call. = FALSE)
  attr(out, "scale_label") <- sc$label
  out
}

#' Match parsed total counts to this run's sample_id, on patient AND timepoint
#'
#' The timepoint half of the key is what separates this from
#' `match_absolute_counts_samples()`. A yield belongs to one draw, so HS5_018
#' at d0 and HS5_018 at d7 are different measurements and must not be allowed
#' to collapse onto each other. Where the sheet carries no timepoint row the
#' join falls back to patient alone, and that is accepted only while it stays
#' unambiguous: a patient with several acquisitions and one unlabelled yield is
#' reported as unmatched rather than broadcast across their timepoints.
#' @param tc Parsed total counts.
#' @param smap Sample map joining sample_id to patient_id.
#' @keywords internal
match_total_counts_samples <- function(tc, smap) {
  if (!"patient_id" %in% names(smap))
    stop("--total-counts needs a patient_id column in the sample sheet to join on",
         call. = FALSE)
  norm <- function(x) norm_id(trimws(as.character(x)))
  s_pat <- norm(smap$patient_id)
  s_tp  <- if ("timepoint" %in% names(smap)) norm(smap$timepoint) else rep("", nrow(smap))

  tc$sample_id <- NA_character_
  for (i in seq_len(nrow(tc))) {
    hit <- which(s_pat == norm(tc$patient_raw[i]))
    if (length(hit) > 1L && !is.na(tc$timepoint[i]))
      hit <- hit[s_tp[hit] == norm(tc$timepoint[i])]
    if (length(hit) == 1L) tc$sample_id[i] <- smap$sample_id[hit]
  }

  bad <- tc[is.na(tc$sample_id), , drop = FALSE]
  if (nrow(bad)) {
    lab <- unique(paste0(bad$patient_raw,
                         ifelse(is.na(bad$timepoint), "",
                                paste0(" @", bad$timepoint))))
    log_msg("  NOTE --total-counts: ", length(lab), " row(s) matched no single ",
            "acquisition in the sample sheet and are dropped (either that ",
            "patient is not in this run, or the patient/timepoint pair is ",
            "ambiguous): ", paste(utils::head(lab, 12), collapse = "; "),
            if (length(lab) > 12) " ..." else "")
  }
  tc <- tc[!is.na(tc$sample_id), , drop = FALSE]

  dup <- unique(tc$sample_id[duplicated(tc$sample_id)])
  if (length(dup))
    stop("--total-counts: ", length(dup), " acquisition(s) received more than one ",
         "total (", paste(dup, collapse = ", "),
         "). A sample cannot have two yields; fix the source sheet.", call. = FALSE)
  tc
}

#' Load, normalize, parse and sample-match a --total-counts input
#'
#' @param path File path.
#' @param smap Sample map joining sample_id to patient_id.
#' @param outdir Directory to write outputs to.
#' @return data.frame(sample_id, patient_raw, timepoint, total_cells), or NULL
#' @export
load_total_counts <- function(path, smap, outdir) {
  if (is.null(smap)) {
    log_msg("NOTE --total-counts requires a sample sheet (--samples or ",
            "--sample-map) for the sample_id join; skipped")
    return(NULL)
  }
  csv_path <- total_counts_to_csv(path, outdir)
  tc <- parse_total_counts_csv(csv_path)
  lab <- attr(tc, "scale_label")
  tc <- match_total_counts_samples(tc, smap)
  if (!nrow(tc)) {
    log_msg("NOTE --total-counts: no rows matched any acquisition; skipped")
    return(NULL)
  }
  attr(tc, "scale_label") <- lab
  log_msg("  --total-counts: ", nrow(tc), " acquisition(s) carry a total cell count",
          " (median ", format(stats::median(tc$total_cells), big.mark = ",",
                              scientific = FALSE), " cells)")
  tc
}

#' Attach absolute cell numbers to a per-sample abundance table
#'
#' Multiplies each population share by that acquisition's total yield. Rows for
#' samples with no total are left NA rather than dropped, so the relative table
#' and the absolute table always have the same shape and the same row order and
#' can be read side by side.
#'
#' WHY `cells_absolute` AND NOT A REPLACEMENT FOR pct_of_gated: the two answer
#' different questions and have different error. The share is measured once, by
#' this pipeline, on events it counted itself. The product carries the yield's
#' error as well, from an instrument this pipeline never saw. Overwriting the
#' share would hide that, so both travel together and the provenance column
#' says so.
#' @param ab Abundance table with sample_id and a share column.
#' @param tc Matched total counts from `load_total_counts()`.
#' @param share_col Column holding the share, in percent.
#' @return `ab` with total_cells, cells_absolute and count_basis added.
#' @export
attach_absolute_cells <- function(ab, tc, share_col = "pct_of_gated") {
  if (is.null(tc) || !nrow(tc) || is.null(ab) || !nrow(ab)) return(ab)
  if (!share_col %in% names(ab)) return(ab)
  tot <- setNames(tc$total_cells, tc$sample_id)
  ab$total_cells <- unname(tot[ab$sample_id])
  ab$cells_absolute <- ab[[share_col]] / 100 * ab$total_cells
  ab$count_basis <- ifelse(is.na(ab$cells_absolute), NA_character_,
                           "dual-platform: share x external total")
  ab
}
