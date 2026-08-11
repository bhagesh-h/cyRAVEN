# SECTION 2c -- INPUT VALIDATION WITHOUT ANALYSIS
# =============================================================================
#
# WHY THIS FILE EXISTS. Everything that makes an input wrong is knowable before
# a single event is read: whether the sheet covers every file, whether the
# specification names markers the panel contains, whether the declared group
# column exists and has more than one level, whether a subject's rows agree.
# None of that needs the data, yet all of it used to surface only partway
# through a run that costs many minutes on a real cohort. A misspelt marker name
# is the single most common cause of an empty frequency table, and the run
# reports it after the reading stage rather than before it.
#
# --check reads FCS HEADERS only, which is a keyword block per file rather than
# a matrix, and reports what the run would do. The point is that it is cheap
# enough to run every time the sheet is edited.

#' Report what a run would do, from headers and metadata alone
#'
#' @param fcs Paths to the FCS files.
#' @param smap Sample map or sheet-derived equivalent; may be NULL.
#' @param sheet Result of [read_samplesheet()], or NULL on the three-file route.
#' @param spec Population specification.
#' @param opt Parsed options.
#' @return invisible TRUE when nothing fatal was found, FALSE otherwise.
#' @keywords internal
report_input_check <- function(fcs, smap, sheet, spec, opt) {
  problems <- character(0)
  note <- function(...) problems <<- c(problems, paste0(...))

  log_msg(length(fcs), " FCS file(s) in ", opt$dir %||% "the file list")

  # ---- markers, from the keyword block of each file -------------------------
  per_file <- lapply(fcs, function(f) {
    kw <- tryCatch(flowCore::read.FCSheader(f)[[1]], error = function(e) NULL)
    if (is.null(kw)) { note("unreadable FCS header: ", basename(f)); return(NULL) }
    np <- suppressWarnings(as.integer(kw[["$PAR"]]))
    if (is.na(np)) return(character(0))
    s <- vapply(seq_len(np), function(i) {
      v <- kw[[paste0("$P", i, "S")]]
      if (is.null(v)) "" else trimws(as.character(v))
    }, character(1))
    s[nzchar(s)]
  })
  names(per_file) <- basename(fcs)
  markers <- sort(unique(unlist(per_file)))
  log_msg("markers resolved from $PnS across all files (", length(markers), "): ",
          paste(markers, collapse = ", "))
  shared <- Reduce(intersect, Filter(length, per_file))
  if (length(shared) && length(shared) < length(markers))
    log_msg("  present in EVERY file (", length(shared), "): ",
            paste(sort(shared), collapse = ", "))

  # ---- does the specification match the panel -------------------------------
  need <- unique(unlist(lapply(spec, names)))
  need <- setdiff(need, c("FSC-A", "FSC-H", "SSC-A", "SSC-H"))
  absent <- setdiff(need, markers)
  log_msg("population specification: ", length(spec), " population(s) needing ",
          length(need), " marker(s)")
  if (length(absent)) {
    note("the specification names ", length(absent),
         " marker(s) no file contains: ", paste(absent, collapse = ", "),
         ". Those populations would be reported UNAVAILABLE. Names must match ",
         "$PnS exactly.")
  } else {
    log_msg("  every named marker is present")
  }
  partial <- setdiff(need, shared)
  partial <- setdiff(partial, absent)
  if (length(partial))
    log_msg("  NOTE present in some files but not all: ",
            paste(partial, collapse = ", "),
            " (those populations are scored only where the marker exists)")

  # ---- coverage and study design -------------------------------------------
  if (is.null(smap)) {
    log_msg("no sample sheet: identifiers will be derived from filenames, and ",
            "no between-group comparison is possible")
  } else {
    log_msg("sheet covers all ", length(fcs), " input file(s)")
    extra_rows <- setdiff(smap$file, basename(fcs))
    if (length(extra_rows))
      log_msg("  NOTE ", length(extra_rows), " sheet row(s) name files not in ",
              "this run (harmless): ", paste(utils::head(extra_rows, 5),
                                             collapse = ", "))
    if ("patient_id" %in% names(smap)) {
      np <- length(unique(smap$patient_id))
      log_msg("  ", np, " subject(s) across ", nrow(smap), " acquisition(s)")
      if (np == nrow(smap))
        log_msg("    one acquisition per subject, so sample-level and ",
                "donor-level tests coincide")
    }
    if (isTRUE(any(smap$is_control)))
      log_msg("  ", sum(smap$is_control), " file(s) declared as controls")
  }

  gcol <- opt$group_column
  if (!is.null(gcol)) {
    src <- if (!is.null(sheet$patients) && gcol %in% names(sheet$patients))
      sheet$patients else smap
    if (is.null(src) || !gcol %in% names(src)) {
      note("--group-column '", gcol, "' is not a column of the sheet. ",
           "Available: ", paste(setdiff(names(src %||% smap), "file"),
                                collapse = ", "))
    } else {
      lv <- table(trimws(as.character(src[[gcol]])))
      lv <- lv[names(lv) != "" & !is.na(names(lv))]
      log_msg("--group-column ", gcol, ": ",
              paste(names(lv), lv, sep = " n=", collapse = "; "))
      if (length(lv) < 2L)
        note("--group-column '", gcol, "' has fewer than two levels, so no ",
             "between-group comparison can be made.")
      if (!is.null(opt$reference_group) && !opt$reference_group %in% names(lv))
        note("--reference-group '", opt$reference_group, "' is not one of the ",
             "levels of ", gcol, ": ", paste(names(lv), collapse = ", "))
      if (any(lv < 3L))
        log_msg("  NOTE group(s) with fewer than 3 samples: ",
                paste(names(lv)[lv < 3L], collapse = ", "),
                " - a rank test on them has little power whatever the effect")
    }
  }

  bcol <- opt[["batch_column", exact = TRUE]]
  if (!is.null(bcol)) {
    src <- if (!is.null(smap) && bcol %in% names(smap)) smap else sheet$patients
    if (is.null(src) || !bcol %in% names(src))
      note("--batch-column '", bcol, "' is not a column of the sheet.")
    else
      log_msg("--batch-column ", bcol, ": ",
              length(unique(src[[bcol]])), " batch(es)")
  }

  if (!is.null(sheet$counts))
    log_msg("absolute counts: ", length(unique(sheet$counts$population)),
            " population(s) for ", length(unique(sheet$counts$sample_id)),
            " sample(s)")
  if (length(sheet$study_columns))
    log_msg("study variable(s) usable as --group-column or --batch-column: ",
            paste(sheet$study_columns, collapse = ", "))

  # ---- verdict --------------------------------------------------------------
  if (length(problems)) {
    log_msg("")
    log_msg(length(problems), " PROBLEM(S) that would affect the run:")
    for (p in problems) log_msg("  - ", p)
    log_msg("")
    log_msg("--check made no changes and analysed nothing.")
    return(invisible(FALSE))
  }
  log_msg("")
  log_msg("no problems found. Remove --check to run the analysis.")
  invisible(TRUE)
}
