# SECTION 6b -- EXTERNAL ABSOLUTE CELL COUNTS (--absolute-counts)
# =============================================================================
#
# WHY THIS EXISTS ALONGSIDE wbc_per_ul (Section 6, abundance_measure() below):
# wbc_per_ul derives an absolute count by multiplying THIS pipeline's own
# pct_of_cd45_pos by a single whole-blood WBC number -- one absolute anchor
# spread across every population under this run's own gating. A directly
# measured absolute count per population (dual-platform counting, counting
# beads, or a clinical flow lab's own analysis) is independent evidence and
# does not depend on this pipeline's gates agreeing with theirs, so it is
# ingested and reported as its own dataset rather than forced into `freq`.
#
# EXPECTED LAYOUT (wide, matches how these are exported from Excel/FlowJo):
#   first row    = population names, one per column, starting anywhere after
#                  column 1 (a units-label spacer column, e.g. "in 1 mL
#                  blood", is auto-detected and skipped: it has a header but
#                  no numeric data in any row)
#   first column = a sample identifier, one row per sample
#   blank rows   = tolerated (cohort block separators); ignored
#
# The identifier is matched against --sample-map by EITHER patient_id OR the
# acquisition filename (case/whitespace-insensitive) because real exports mix
# both conventions row to row -- not every collaborator names every row the
# same way the FCS files were named.

#' Normalize any --absolute-counts input down to one CSV on disk
#'
#' WHY funnel every format through a written CSV rather than branch the parser
#' on format: the source layout is irregular (blank separator rows, a
#' units-label spacer column) and that scan should exist exactly once,
#' regardless of whether the input arrived as .xlsx, .csv, or .tsv. CSV, not
#' TSV, because every other tabular input this pipeline reads (--sample-map,
#' --patient-table) is already CSV -- one delimiter convention throughout
#' rather than a second one just for this file. Writing the intermediate into
#' outdir (not a tempfile) also gives the user something to open when a match
#' or a number looks wrong.
#' @param path File path.
#' @param outdir Directory to write outputs to.
#' @keywords internal
absolute_counts_to_csv <- function(path, outdir) {
  ext <- tolower(tools::file_ext(path))
  grid <- switch(ext,
    xlsx = {
      if (!requireNamespace("readxl", quietly = TRUE))
        stop("--absolute-counts points at an .xlsx file, which needs the 'readxl' ",
             "package.\n  run:  install.packages(\"readxl\")", call. = FALSE)
      # col_types = "text" everywhere: with col_names = FALSE the header row
      # (text) sits in the same column as the data rows (numeric), and letting
      # readxl guess a column type per column risks it guessing "numeric" and
      # turning the header cell into NA.
      as.data.frame(readxl::read_excel(path, col_names = FALSE, col_types = "text",
                                       .name_repair = "minimal"),
                    stringsAsFactors = FALSE)
    },
    tsv = ,
    txt = utils::read.delim(path, header = FALSE, colClasses = "character", check.names = FALSE,
                     blank.lines.skip = FALSE, fill = TRUE, quote = "",
                     stringsAsFactors = FALSE),
    NULL)
  if (is.null(grid) && ext != "csv")
    stop("--absolute-counts: unrecognized extension '.", ext,
         "' (expected .xlsx, .csv, or .tsv)", call. = FALSE)
  csv_path <- file.path(outdir, "absolute_counts_raw.csv")
  if (is.null(grid)) {
    file.copy(path, csv_path, overwrite = TRUE)
  } else {
    grid[] <- lapply(grid, function(col) ifelse(is.na(col), "", trimws(as.character(col))))
    utils::write.table(grid, csv_path, sep = ",", row.names = FALSE, col.names = FALSE,
               na = "", qmethod = "double", fileEncoding = "UTF-8")
  }
  log_msg("  --absolute-counts: normalized ", basename(path), " -> ", csv_path)
  csv_path
}

#' Parse the normalized CSV into long format: sample_raw, population, cells_per_ul
#'
#' Auto-detects which columns are real population data (header text AND at
#' least one numeric value somewhere below it) versus a units-label spacer
#' column (header text, never any numeric value below it -- e.g. "in 1 mL
#' blood"), and uses that spacer text, if found, to convert per-mL counts to
#' the cells_per_ul convention this pipeline's own wbc_per_ul route uses.
#' Falls back to assuming cells/uL already, LOUDLY, when no unit is found --
#' silently guessing a 1000x factor wrong is worse than asking the user to
#' check.
#' @param csv_path The csv path.
#' @keywords internal
parse_absolute_counts_csv <- function(csv_path) {
  grid <- as.matrix(read.csv(csv_path, header = FALSE, colClasses = "character",
                             check.names = FALSE, na.strings = character(0),
                             blank.lines.skip = FALSE, fill = TRUE))
  if (nrow(grid) < 2L)
    stop("--absolute-counts: no data rows found below the header row", call. = FALSE)

  header <- trimws(grid[1, ])
  body   <- grid[-1, , drop = FALSE]
  is_num <- function(x) suppressWarnings(!is.na(as.numeric(x)) & nzchar(trimws(x)))

  has_data <- vapply(seq_len(ncol(body)), function(j) any(is_num(body[, j])), logical(1))
  pop_cols <- which(has_data & nzchar(header))
  if (!length(pop_cols))
    stop("--absolute-counts: no column has both a header label and a numeric value",
         call. = FALSE)

  spacer <- header[nzchar(header)]
  spacer <- spacer[!seq_along(header)[nzchar(header)] %in% pop_cols]
  spacer_txt <- paste(spacer, collapse = " ")
  scale_to_ul <- if (grepl("(?i)micro.?lit|\\bu[l]\\b|\u00b5\\s*l", spacer_txt, perl = TRUE)) 1
    else if (grepl("(?i)\\bm[l]\\b|milli.?lit", spacer_txt, perl = TRUE)) 1 / 1000
    else NA_real_
  if (is.na(scale_to_ul)) {
    log_msg("  NOTE --absolute-counts: no mL/\u00b5L unit label found (looked at: '",
            spacer_txt, "'); assuming values are already cells/\u00b5L. VERIFY this",
            " before trusting absolute_counts.png/csv.")
    scale_to_ul <- 1
  } else if (scale_to_ul != 1) {
    log_msg("  --absolute-counts: header label '", trimws(spacer_txt),
            "' -> treating values as cells/mL, converting to cells/\u00b5L (\u00f7 1000)")
  }

  rows <- lapply(seq_len(nrow(body)), function(i) {
    # unname(): a single-cell matrix extraction still carries a name from the
    # column dimnames (R quirk), which data.frame() below would otherwise
    # warn about and discard anyway ("row names were found from a short
    # variable")
    sample_raw <- unname(trimws(body[i, 1]))
    if (!nzchar(sample_raw)) return(NULL)          # blank cohort-separator row
    vals <- suppressWarnings(as.numeric(body[i, pop_cols]))
    keep <- !is.na(vals)
    if (!any(keep)) return(NULL)
    data.frame(sample_raw = sample_raw, population = unname(header[pop_cols][keep]),
              cells_per_ul = unname(vals[keep]) * scale_to_ul, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out))
    stop("--absolute-counts: every data row was blank or non-numeric", call. = FALSE)
  out
}

#' Resolve sample_id -> an arbitrary patient-table column, via the sample map
#'
#' Shared by the between-group figures (grouping column, e.g. "cohort") and by
#' --absolute-counts sample matching (see below): both need patient_id as the
#' join key between the patient table and the sample map's per-file rows.
#' @param patients Patient metadata table, keyed by patient_id.
#' @param smap Sample map joining sample_id to patient_id.
#' @param gcol The gcol.
#' @keywords internal
resolve_group_of <- function(patients, smap, gcol) {
  if (is.null(patients) || !gcol %in% names(patients) ||
      is.null(smap) || !"patient_id" %in% names(smap)) return(NULL)
  k <- smap[, intersect(c("sample_id", "well", "patient_id"), names(smap)), drop = FALSE]
  k$sample_id <- k$sample_id %||% k$well
  pg <- patients[, c("patient_id", gcol)]
  pg <- pg[!is.na(pg[[gcol]]) & nzchar(as.character(pg[[gcol]])), , drop = FALSE]
  if (!nrow(pg)) return(NULL)
  k$.j <- norm_id(k$patient_id); pg$.j <- norm_id(pg$patient_id)
  m <- merge(k, pg[, c(".j", gcol)], by = ".j", all.x = TRUE)
  m <- m[!is.na(m[[gcol]]) & !is.na(m$sample_id), , drop = FALSE]
  if (!nrow(m)) return(NULL)
  setNames(as.character(m[[gcol]]), m$sample_id)
}

#' Match each --absolute-counts row to this run's sample_id
#'
#' Tries patient_id first, then the acquisition filename (both
#' case/whitespace-insensitive, extension and a trailing " copy" stripped),
#' because real exports are not internally consistent about which one labels
#' a given row (seen in practice: most rows use the filename, a handful use
#' the patient_id instead). Whichever key matches wins; rows matching neither
#' are dropped with a NOTE naming them, rather than silently guessed at.
#' @param ac The ac.
#' @param smap Sample map joining sample_id to patient_id.
#' @keywords internal
match_absolute_counts_samples <- function(ac, smap) {
  strip <- function(x) norm_id(sub("\\s+copy$", "", sub("[.](fcs|FCS)$", "", trimws(x))))
  key <- strip(ac$sample_raw)

  by_patient <- if ("patient_id" %in% names(smap))
    setNames(smap$sample_id, strip(smap$patient_id)) else NULL

  # File match is a SUFFIX match, not exact: real acquisition filenames carry
  # a batch/folder prefix the counts sheet's row label never repeats (e.g.
  # "Blood samples_HC-13;2.fcs" vs a row labelled "HC-13;2.fcs"). A delimiter
  # (_, space or -) must sit immediately before the matched suffix, so "HC-1"
  # cannot accidentally match a file ending in "...HC-11"; an ambiguous key
  # matching more than one file is treated as no match rather than guessed.
  file_stripped <- strip(basename(smap$file))
  find_file_match <- function(k) {
    hit <- which(file_stripped == k)
    if (!length(hit)) {
      ends <- endsWith(file_stripped, k)
      if (any(ends)) {
        cand <- which(ends)
        ok <- vapply(cand, function(i) {
          h <- file_stripped[i]
          if (nchar(h) == nchar(k)) return(TRUE)
          substr(h, nchar(h) - nchar(k), nchar(h) - nchar(k)) %in% c("_", " ", "-")
        }, logical(1))
        hit <- cand[ok]
      }
    }
    if (length(hit) == 1L) smap$sample_id[hit] else NA_character_
  }

  ac$sample_id <- if (!is.null(by_patient)) unname(by_patient[key]) else NA_character_
  miss <- is.na(ac$sample_id)
  if (any(miss)) ac$sample_id[miss] <- vapply(key[miss], find_file_match, character(1))

  unmatched <- unique(ac$sample_raw[is.na(ac$sample_id)])
  if (length(unmatched))
    log_msg("  NOTE --absolute-counts: ", length(unmatched), " label(s) matched no sample",
            " in --sample-map (check for typos) and are dropped: ",
            paste(unmatched, collapse = "; "))
  ac[!is.na(ac$sample_id), , drop = FALSE]
}

#' Load, normalize, parse and sample-match an --absolute-counts input
#' @param path File path.
#' @param smap Sample map joining sample_id to patient_id.
#' @param outdir Directory to write outputs to.
#' @return data.frame(sample_id, sample_raw, population, cells_per_ul) or NULL
#' @export
load_absolute_counts <- function(path, smap, outdir) {
  if (is.null(smap)) {
    log_msg("NOTE --absolute-counts requires --sample-map (for the sample_id join); skipped")
    return(NULL)
  }
  csv_path <- absolute_counts_to_csv(path, outdir)
  ac <- parse_absolute_counts_csv(csv_path)
  ac <- match_absolute_counts_samples(ac, smap)
  if (!nrow(ac)) {
    log_msg("NOTE --absolute-counts: no rows matched any sample; skipped")
    return(NULL)
  }
  log_msg("  --absolute-counts: ", length(unique(ac$sample_id)), " sample(s), ",
          length(unique(ac$population)), " population(s) matched")
  ac
}

# =============================================================================
