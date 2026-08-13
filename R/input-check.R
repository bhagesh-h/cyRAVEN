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
    # Resolve the name the SAME WAY the run will, via marker_symbol(): a
    # spectral $PnS reads "CD45 : SparkUV-387 - Area" and the run strips from
    # " : " onward. Reporting the raw string here made --check announce that
    # every population named a marker no file contained, on data where the run
    # resolved all of them -- a false alarm from the one command whose purpose
    # is to catch real ones.
    # kw_value(), not kw[[...]]: the guard that used to be here could not fire.
    # See kw_value() for why a missing $PnS aborted the command instead.
    s <- vapply(seq_len(np), function(i) kw_value(kw, i, "S"), character(1))
    n <- vapply(seq_len(np), function(i) kw_value(kw, i, "N"), character(1))
    s <- marker_symbol(s, n)
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

#' One `$Pn*` keyword, or "" when the file does not carry it
#'
#' read.FCSHeader returns a NAMED CHARACTER VECTOR rather than a list, and `[[`
#' on one raises "subscript out of bounds" for a name that is absent, where the
#' same call on a list returns NULL. An `if (is.null(v))` guard therefore never
#' runs: the error is raised before it is reached.
#'
#' A parameter with no `$PnS` is ordinary rather than exceptional. Time channels
#' routinely carry none, and some instruments omit it for scatter. Reading a
#' cohort containing one used to abort the whole command.
#'
#' @param kw Keyword vector from [flowCore::read.FCSheader()].
#' @param i Parameter index.
#' @param suffix "N" or "S".
#' @return Trimmed character scalar, "" when absent or NA.
#' @keywords internal
kw_value <- function(kw, i, suffix) {
  key <- paste0("$P", i, suffix)
  if (!key %in% names(kw)) return("")
  v <- kw[[key]]
  if (is.null(v) || is.na(v)) "" else trimws(as.character(v))
}

# --list-channels ------------------------------------------------------------
#
# WHY THIS IS A FLAG AND NOT A SNIPPET. Writing a specification requires knowing
# the name to put in it, and on a spectral panel that name is not what the
# instrument recorded: $PnS reads "CD45 : SparkUV-387 - Area" and the run uses
# "CD45". Getting it wrong produces a full set of tables containing zeros, which
# is the most common way a run fails.
#
# That left every user dumping the keyword block with an ad-hoc Rscript one-liner
# against the FIRST file. Two things are wrong with that. It shows the raw
# keyword rather than the name the run resolves, so it does not answer the
# question being asked; and reading one file cannot detect the panel differing
# between files, which is the failure worth catching before a run rather than
# during one.

#' List the acquisition parameters and the names the run resolves them to
#'
#' One row per parameter: index, the detector name (`$PnN`), the description
#' (`$PnS`), the symbol the run resolves it to, and how the run will use it. The
#' resolved symbol is the string a population specification has to name.
#'
#' Reads headers only, so it costs nothing on large files, and reads every file
#' rather than the first so a panel that differs between files is reported.
#'
#' @param fcs Paths to the FCS files.
#' @return invisible data.frame of the reference file's parameters.
#' @keywords internal
list_channels <- function(fcs) {
  read_params <- function(f) {
    kw <- tryCatch(flowCore::read.FCSheader(f)[[1]], error = function(e) NULL)
    if (is.null(kw)) return(NULL)
    np <- suppressWarnings(as.integer(kw[["$PAR"]]))
    if (is.na(np)) return(NULL)
    get_kw <- function(suffix) vapply(seq_len(np), function(i) kw_value(kw, i, suffix),
                                      character(1))
    n <- get_kw("N")
    s <- get_kw("S")
    data.frame(idx = seq_len(np), name = n, desc = s,
               symbol = marker_symbol(s, n),
               stringsAsFactors = FALSE, check.names = FALSE)
  }

  # The role column mirrors the rules in read_fcs_resolved(): area channels
  # carry the markers, height and width are redundant duplicates of area that
  # would double-weight a marker in the embedding, and scatter is matched by its
  # own patterns because it is used for gating rather than as a marker.
  role_of <- function(nm, ds) {
    scatter <- c("FSC.*A(rea)?$", "FSC.*H(eight)?$", "SSC.*A(rea)?$", "SSC.*H(eight)?$")
    is_sc <- Reduce(`|`, lapply(scatter, function(p) grepl(p, nm) | grepl(p, ds)))
    is_area <- grepl("(^|[^A-Za-z])A$|\\.A$|-A$|Area$", nm) | grepl("- Area$", ds)
    is_time <- grepl("^Time", nm, ignore.case = TRUE) |
               grepl("^Time", ds, ignore.case = TRUE)
    out <- rep("ignored, not an area channel", length(nm))
    out[is_area] <- "marker"
    out[is_sc]   <- "scatter, used for gating"
    out[is_time] <- "ignored, acquisition time"
    out
  }

  ref <- NULL
  for (f in fcs) { ref <- read_params(f); if (!is.null(ref)) { ref_file <- f; break } }
  if (is.null(ref)) {
    log_msg("no readable FCS header in ", length(fcs), " file(s)")
    return(invisible(NULL))
  }
  ref$role <- role_of(ref$name, ref$desc)

  log_msg(length(fcs), " file(s); parameters listed from ", basename(ref_file))
  log_msg("")
  cat(sprintf("%-5s %-22s %-38s %-16s %s\n",
              "idx", "$PnN", "$PnS", "resolves to", "role"))
  cat(strrep("-", 110), "\n", sep = "")
  for (i in seq_len(nrow(ref)))
    cat(sprintf("P%-4d %-22s %-38s %-16s %s\n", ref$idx[i],
                substr(ref$name[i], 1, 22), substr(ref$desc[i], 1, 38),
                substr(ref$symbol[i], 1, 16), ref$role[i]))
  cat("\n")

  usable <- sort(ref$symbol[ref$role == "marker"])
  log_msg(length(usable), " marker name(s) to use in a specification:")
  log_msg("  ", paste(usable, collapse = ", "))

  # ---- does every file carry the same panel ---------------------------------
  if (length(fcs) > 1L) {
    others <- lapply(fcs, function(f) {
      p <- read_params(f)
      if (is.null(p)) return(NULL)
      sort(p$symbol[role_of(p$name, p$desc) == "marker"])
    })
    names(others) <- basename(fcs)
    bad <- names(others)[vapply(others, is.null, logical(1))]
    others <- Filter(Negate(is.null), others)
    same <- vapply(others, function(x) identical(x, usable), logical(1))
    log_msg("")
    if (all(same)) {
      log_msg("every file carries the same ", length(usable), " marker(s)")
    } else {
      log_msg("PANEL DIFFERS between files. ", sum(!same), " of ", length(others),
              " do not match ", basename(ref_file), ":")
      for (nm in names(others)[!same]) {
        miss <- setdiff(usable, others[[nm]])
        extra <- setdiff(others[[nm]], usable)
        log_msg("  ", nm,
                if (length(miss)) paste0("  missing: ", paste(miss, collapse = ", ")) else "",
                if (length(extra)) paste0("  extra: ", paste(extra, collapse = ", ")) else "")
      }
      log_msg("A population naming a marker some files lack is UNAVAILABLE in ",
              "those files, not an error.")
      # NAME THE FIX, NOT JUST THE SYMPTOM.
      #
      # A differing panel is not itself the problem. The problem is what it
      # causes: the fingerprint is the set of resolved marker names, so files
      # that differ become separate panels, each with its own cofactor,
      # thresholds and embedding, and a cohort quietly stops being one cohort.
      # Twelve comparable files have become seven panels this way.
      #
      # Two causes are mechanical enough to detect without judgement, so they
      # are named here with the exact flag that merges the cohort. Anything else
      # is reported as a genuine panel difference and left alone, because that
      # is a decision about the experiment rather than about the file format.
      .all <- unique(unlist(others))
      .n_with <- vapply(.all, function(s)
        sum(vapply(others, function(x) s %in% x, logical(1))), integer(1))
      .varies <- names(.n_with)[.n_with < length(others)]
      # 1. Spectral unmixing writes extracted autofluorescence back as extra
      #    channels, and how many appear varies per acquisition. Not stains.
      .af <- grep("^\\[?AF([ _-]|$)|autofluor|unmix", .varies,
                  ignore.case = TRUE, value = TRUE)
      # 2. A real marker stained in only a minority of files. Legitimate, but it
      #    still splits the cohort, so the choice has to be made deliberately.
      .minority <- setdiff(.varies[.n_with[.varies] < length(others) / 2], .af)
      if (length(.af) || length(.minority)) {
        log_msg("")
        log_msg("SUGGESTED FIX. These channels are why the files split:")
        if (length(.af))
          log_msg("  autofluorescence/unmixing artefacts, not stains: ",
                  paste(.af, collapse = ", "))
        for (s in .minority)
          log_msg("  '", s, "' is present in only ", .n_with[[s]], " of ",
                  length(others), " file(s)")
        .pat <- c(if (length(.af)) "[AF color*", .minority)
        log_msg("  Dropping them merges the cohort into one panel:")
        log_msg("    --ignore-channels '", paste(.pat, collapse = ","), "'")
        log_msg("  Quote it, or the shell expands [ and * itself. Dropped ",
                "before the fingerprint is computed, which is the only place ",
                "it works.")
        if (length(.minority))
          log_msg("  Keep a minority marker only if you accept that the answer ",
                  "covers those files alone, which is not a cohort-wide result.")
      }
    }
    if (length(bad))
      log_msg("unreadable header(s): ", paste(bad, collapse = ", "))
  }

  log_msg("")
  log_msg("--list-channels read headers only and analysed nothing.")
  invisible(ref)
}
