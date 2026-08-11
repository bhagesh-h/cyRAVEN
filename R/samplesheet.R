# SECTION 2b -- THE UNIFIED SAMPLE SHEET
# =============================================================================
#
# WHY THIS FILE EXISTS. A run used to need up to three separate tables: a sample
# map keyed by filename, a patient table keyed by patient_id, and an absolute-count
# export keyed by whichever of the two the counting instrument happened to write.
# Each is a different shape with a different key, and the joins between them are
# where cohorts get mislabelled: a patient present in one table and absent from
# another produces no error, only an empty covariate panel or a silently dropped
# count. The user assembling them has to hold three schemas in mind at once and
# has no single place to look to answer "what is sample X".
#
# The sample sheet is one row per acquired file, carrying every fact about that
# file: what it is, whose it is, what group and batch it belongs to, and any
# externally measured counts for it. One key, one file, no joins.
#
# WHAT IT DOES NOT ABSORB. The population specification, thresholds, colours and
# translations stay in the YAML. Those describe the ANALYSIS and are shared by
# every sample; folding them into a per-sample table would repeat one decision
# on every row and invite the rows to disagree about it. The pair is therefore
# one CSV describing the samples and one YAML describing the analysis, and
# between them they specify a run completely.
#
# HOW IT RELATES TO THE THREE-FILE ROUTE. This reader splits the sheet into
# exactly the three structures the pipeline already consumed, with the same
# columns and the same coercions, by calling the same helpers. It is a different
# way to supply the same facts, not a different analysis. The equivalence is
# asserted by a test that runs one study both ways and compares every output.
#
# THE ONE HAZARD THE THREE-FILE ROUTE DID NOT HAVE. Subject attributes are
# properties of a patient, but the sheet has a row per FILE, so a patient with
# several acquisitions repeats them. Two rows of the same patient can therefore
# disagree about that patient's sex. A reader that took the first value would
# silently pick one. read_samplesheet() reports every such disagreement and
# stops, because there is no defensible way to choose.

#' Reserved column names of the unified sample sheet
#'
#' `acquisition` names the file and how it was acquired; `subject` names
#' properties of the patient, which repeat across that patient's rows. Any
#' column not listed here and not carrying the count prefix is kept as a
#' study variable, usable as `--group-column`, `--batch-column` or a covariate.
#'
#' @return list(acquisition, subject, required, count_prefix)
#' @export
samplesheet_columns <- function() {
  list(
    required    = "file",
    acquisition = c("file", "well", "sample_id", "patient_id", "panel",
                    "timepoint", "is_control", "fmo_for", "control_group"),
    # Recognised by name here so that the sheet needs no translation step. A
    # German header still works: it is resolved through default_column_map()
    # exactly as in a standalone patient table.
    subject     = c("patient_id", "date_of_birth", "sex", "age_years",
                    "height_cm", "weight_kg", "infection_focus", "cohort",
                    "wbc_per_ul"),
    count_prefix = "count.")
}

#' Resolve one canonical column name against a sheet's headers
#'
#' Accepts the canonical spelling and every alias in `column_map`, so a sheet
#' exported in German resolves without being translated first.
#' @param nms character vector of the sheet's column names.
#' @param canon canonical name to look for.
#' @param column_map The column map. Default `default_column_map()`.
#' @return index into `nms`, or NA.
#' @keywords internal
match_sheet_column <- function(nms, canon, column_map = default_column_map()) {
  cands <- unique(c(canon, column_map[[canon]]))
  for (cd in cands) {
    j <- which(tolower(trimws(nms)) == tolower(cd))
    if (length(j)) return(j[1])
  }
  NA_integer_
}

#' Check that a subject attribute does not disagree between a patient's rows
#'
#' @param df the sheet.
#' @param key character vector of patient identifiers, one per row.
#' @param cols subject columns present in the sheet.
#' @return invisible TRUE, or stops naming every conflict.
#' @keywords internal
assert_subject_consistent <- function(df, key, cols) {
  bad <- character(0)
  for (cn in cols) {
    v <- trimws(as.character(df[[cn]]))
    v[is.na(v) | !nzchar(v)] <- NA_character_
    sp <- split(v, key)
    for (p in names(sp)) {
      u <- unique(stats::na.omit(sp[[p]]))
      if (length(u) > 1L)
        bad <- c(bad, sprintf("  %s / %s: %s", p, cn,
                              paste(sort(u), collapse = " vs ")))
    }
  }
  if (length(bad))
    stop("the sample sheet gives conflicting values for the same subject.\n",
         "  A subject attribute repeats on every row of that subject, so the ",
         "rows must agree.\n",
         paste(utils::head(bad, 20), collapse = "\n"),
         if (length(bad) > 20) sprintf("\n  ... and %d more", length(bad) - 20) else "",
         call. = FALSE)
  invisible(TRUE)
}

#' Read the unified sample sheet
#'
#' One row per acquired file. Splits into the three structures the pipeline
#' consumes, applying the same coercions as the separate loaders so that a study
#' analysed either way produces the same numbers.
#'
#' @param path Path to the sheet (CSV or TSV; delimiter and encoding detected).
#' @param fcs_files The fcs files being analysed, for the coverage check.
#' @param column_map The column map. Default `default_column_map()`.
#' @param value_map The value map. Default `default_value_map()`.
#' @param reference_date The reference date. Default `Sys.Date()`.
#' @param count_unit "cells/uL" or "cells/mL"; the latter is divided by 1000.
#' @return list(smap, patients, counts, study_columns)
#' @export
read_samplesheet <- function(path, fcs_files, column_map = default_column_map(),
                             value_map = default_value_map(),
                             reference_date = Sys.Date(),
                             count_unit = "cells/uL") {
  S <- samplesheet_columns()

  # Same delimiter and encoding detection as the patient table: a sheet
  # assembled in a German Excel is semicolon-separated latin1, and failing on it
  # would push the user back to hand-editing the very file this format exists to
  # simplify.
  raw <- readLines(path, warn = FALSE, n = 5L)
  delim <- if (mean(lengths(regmatches(raw, gregexpr(";", raw)))) >
                mean(lengths(regmatches(raw, gregexpr(",", raw))))) ";" else ","
  if (grepl("[.]tsv$|[.]tab$", path, ignore.case = TRUE)) delim <- "\t"
  test <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
  # U+FFFD is the replacement character. Its presence means the bytes were not
  # valid UTF-8 and something has already substituted for them, which is the case
  # validUTF8() below can no longer see. Built with intToUtf8() rather than
  # written literally: a portable package must have ASCII-only R source.
  repl_char <- intToUtf8(0xFFFDL)
  enc <- if (is.null(test) || any(grepl(repl_char, test, fixed = TRUE)) ||
             any(!validUTF8(test)))
    "latin1" else "UTF-8"
  df <- read.csv(path, sep = delim, fileEncoding = enc, check.names = FALSE,
                 stringsAsFactors = FALSE, na.strings = c("", "NA", "na"))
  names(df) <- trimws(names(df))
  if (!nrow(df)) stop("the sample sheet has no data rows", call. = FALSE)
  log_msg("  sample sheet: ", nrow(df), " row(s), ", ncol(df), " column(s) (",
          delim, "-separated, ", enc, ")")

  # ---- acquisition half, which becomes the sample map -----------------------
  jf <- match_sheet_column(names(df), "file", list(file = c("file", "filename",
                                                           "File", "Dateiname")))
  if (is.na(jf))
    stop("the sample sheet must contain a 'file' column naming each FCS file",
         call. = FALSE)
  smap <- data.frame(file = basename(trimws(as.character(df[[jf]]))),
                     stringsAsFactors = FALSE)
  for (cn in setdiff(S$acquisition, "file")) {
    j <- which(tolower(names(df)) == cn)
    if (length(j)) smap[[cn]] <- df[[j[1]]]
  }
  # patient_id is the join key for everything subject-level, and is resolved
  # through the alias map so a sheet headed "Patient ID" works unchanged.
  if (!"patient_id" %in% names(smap)) {
    jp <- match_sheet_column(names(df), "patient_id", column_map)
    if (!is.na(jp)) smap$patient_id <- df[[jp]]
  }
  if (any(duplicated(smap$file)))
    stop("duplicate file entries in the sample sheet: ",
         paste(unique(smap$file[duplicated(smap$file)]), collapse = ", "),
         call. = FALSE)
  unmapped <- setdiff(basename(fcs_files), smap$file)
  if (length(unmapped))
    stop("these input files are not in the sample sheet: ",
         paste(unmapped, collapse = ", "),
         "\n  Add a row for each, or run --write-samples to regenerate the ",
         "sheet from this directory.", call. = FALSE)
  if ("is_control" %in% names(smap))
    smap$is_control <- tolower(trimws(as.character(smap$is_control))) %in%
                       c("true", "t", "yes", "y", "1")

  sid <- as.character(smap$sample_id %||% smap$well %||% smap$file)

  # ---- subject half, which becomes the patient table ------------------------
  patients <- NULL
  subj <- setdiff(S$subject, "patient_id")
  present <- stats::setNames(vapply(subj, function(cn)
    match_sheet_column(names(df), cn, column_map), integer(1)), subj)
  present <- present[!is.na(present)]
  if (length(present)) {
    if (!"patient_id" %in% names(smap)) {
      # Without a patient key every row is its own subject. That is correct for
      # a one-acquisition-per-subject study and is stated rather than assumed,
      # because it also silently disables the donor-level aggregation the
      # statistics rely on when it is NOT true.
      log_msg("  NOTE the sheet carries subject columns but no patient_id; ",
              "treating each file as its own subject")
      smap$patient_id <- sid
    }
    pkey <- as.character(smap$patient_id)
    sub_df <- df[, unname(present), drop = FALSE]
    names(sub_df) <- names(present)
    assert_subject_consistent(sub_df, pkey, names(sub_df))
    # One row per patient: take each attribute's first non-missing value, which
    # is unambiguous only because the check above proved the rows agree.
    first_present <- function(x) {
      x <- trimws(as.character(x)); x[!nzchar(x)] <- NA_character_
      u <- stats::na.omit(x); if (length(u)) u[1] else NA_character_
    }
    ids <- unique(pkey)
    patients <- data.frame(patient_id = ids, stringsAsFactors = FALSE)
    for (cn in names(sub_df))
      patients[[cn]] <- vapply(ids, function(p)
        first_present(sub_df[[cn]][pkey == p]), character(1))
    patients <- normalise_patient_columns(patients, value_map = value_map,
                                          reference_date = reference_date)
    log_msg("  ", nrow(patients), " subject(s) with ", ncol(patients) - 1L,
            " attribute(s): ", paste(setdiff(names(patients), "patient_id"),
                                     collapse = ", "))
  }

  # ---- count columns, which become the absolute-count table -----------------
  counts <- NULL
  cc <- grep(paste0("^", S$count_prefix), names(df), ignore.case = TRUE)
  if (length(cc)) {
    scale_to_ul <- if (grepl("(?i)m[l]|milli", count_unit, perl = TRUE)) 1 / 1000 else 1
    rows <- lapply(cc, function(j) {
      pop <- sub(paste0("(?i)^", S$count_prefix), "", names(df)[j], perl = TRUE)
      v <- suppressWarnings(as.numeric(gsub(",", ".",
             trimws(as.character(df[[j]])))))
      ok <- !is.na(v)
      if (!any(ok)) return(NULL)
      data.frame(sample_id = sid[ok], sample_raw = sid[ok], population = pop,
                 cells_per_ul = v[ok] * scale_to_ul, stringsAsFactors = FALSE)
    })
    counts <- do.call(rbind, rows)
    if (!is.null(counts))
      log_msg("  absolute counts: ", length(unique(counts$population)),
              " population(s) across ", length(unique(counts$sample_id)),
              " sample(s), read as ", count_unit,
              if (scale_to_ul != 1) " (converted to cells/uL)" else "")
  }

  # ---- everything else is a study variable ----------------------------------
  used <- unique(c(names(df)[unname(present)], names(df)[cc], names(df)[jf],
                   names(df)[tolower(names(df)) %in% S$acquisition]))
  extra <- setdiff(names(df), used)
  if (length(extra)) {
    for (cn in extra) if (!cn %in% names(smap)) smap[[cn]] <- df[[cn]]
    log_msg("  study variable(s) available to --group-column / --batch-column: ",
            paste(extra, collapse = ", "))
  }

  list(smap = smap, patients = patients, counts = counts, study_columns = extra)
}

#' Write a sample-sheet template covering every input file
#'
#' Emits the reserved columns with the filename-derived identifiers filled in
#' and the rest blank, so the user edits a sheet that already accounts for every
#' file rather than assembling one and discovering an omission at run time.
#'
#' @param fcs_files The fcs files.
#' @param path File path to write.
#' @param sample_ids Optional identifiers; derived from the filenames otherwise.
#' @param populations Optional population names; one `count.` column each.
#' @keywords internal
write_samplesheet_template <- function(fcs_files, path, sample_ids = NULL,
                                       populations = character(0)) {
  S <- samplesheet_columns()
  ids <- sample_ids %||% vapply(basename(fcs_files), derive_sample_id, "", kw = NULL)
  tmpl <- data.frame(
    file = basename(fcs_files), sample_id = ids, patient_id = ids,
    is_control = "FALSE", timepoint = "", panel = "",
    cohort = "", sex = "", age_years = "",
    stringsAsFactors = FALSE)
  for (p in populations) tmpl[[paste0(S$count_prefix, p)]] <- ""
  write.csv(tmpl, path, row.names = FALSE, na = "")
  log_msg("wrote sample sheet template: ", path, " (", nrow(tmpl), " row(s))")
  log_msg("  fill in patient_id, cohort and any study variable, then re-run ",
          "with --samples ", path)
  log_msg("  every column is documented in the Inputs article; ",
          "--check validates the sheet without running the analysis")
  invisible(tmpl)
}
