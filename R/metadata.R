# SECTION 6 -- GERMAN -> ENGLISH PATIENT METADATA
# =============================================================================

#' Default German -> English column mapping (the user's requested seven)
#' @keywords internal
default_column_map <- function() {
  # Each canonical name maps to the header spellings seen in the wild, in BOTH
  # source languages AND the canonical name itself.
  #
  # WHY include the canonical names: the script's own output uses them, and a
  # table that has already been translated (by an earlier run, by a collaborator,
  # or written in English from the start) must round-trip. Without the canonical
  # spellings such a table loses every column silently -- the run continues, the
  # covariate panels come out empty, and nothing reports why.
  list(patient_id = c("patient_id", "ID", "Patient", "PatientID", "Patienten-ID",
                      "Patient ID", "SubjectID", "Subject"),
       date_of_birth = c("date_of_birth", "Geburtsdatum", "GebDatum",
                         "Geburtstag", "DOB", "Birthdate", "Date of Birth"),
       sex = c("sex", "Geschlecht", "Sex", "Gender"),
       infection_focus = c("infection_focus", "Fokus", "Infektionsfokus",
                           "Focus", "Infection focus"),
       age_years = c("age_years", "Alter", "Age", "Age (years)", "age"),
       height_cm = c("height_cm", "Groesse", "Gr\u00f6sse", "Gr\u00f6\u00dfe",
                     "Grosse", "Height", "Height (cm)"),
       weight_kg = c("weight_kg", "Gewicht", "Weight", "Weight (kg)"),
       # Study-design grouping. WHY it belongs in the DEFAULT map: a cohort,
       # genotype or treatment-arm column is the independent variable of most
       # case-control designs, and omitting it from the map means the one
       # comparison the experiment exists to make is the one panel not drawn.
       cohort = c("cohort", "Cohort", "Kohorte", "Group", "Gruppe", "Genotype",
                  "Diagnosis", "Diagnose", "Arm"),
       # Absolute leukocyte concentration from an external blood count, in
       # cells/uL. OPTIONAL, and the only route to absolute cell numbers.
       #
       # WHY IT CANNOT BE DERIVED FROM THE FCS FILES: a cytometer reports how many
       # events it recorded, not the volume it sampled. Event counts therefore
       # measure how long the operator ran the tube, and differ several-fold
       # between samples of identical composition. Frequencies are also
       # COMPOSITIONAL: one population rising necessarily pushes the others down,
       # so "% of CD45+" cannot distinguish a real expansion of one lineage from a
       # contraction of another. Supplying WBC/uL from a haemogram (or a counting-
       # bead calibration) converts each frequency into cells/uL, which is the
       # quantity a clinical claim about lymphopenia or NK-cell deficiency needs.
       wbc_per_ul = c("wbc_per_ul", "WBC", "WBC/uL", "WBC_per_uL", "wbc",
                      "Leukozyten", "Leukozyten/uL", "Leukocytes",
                      "leukocytes_per_ul", "absolute_wbc", "Leukos"))
}

#' Default value translations; extend in the YAML config without touching code
#' @keywords internal
default_value_map <- function() {
  # Identity entries for the already-English target values are deliberate: the
  # translation is idempotent, so a table that is already English (or was produced
  # by an earlier run of this script) passes through silently instead of emitting
  # an "untranslated value" warning for every row. A warning that fires on correct
  # input trains the reader to ignore it.
  list(sex = list(m = "male", w = "female", f = "female", weiblich = "female",
                  "m\u00e4nnlich" = "male", maennlich = "male",
                  male = "male", female = "female",
                  man = "male", woman = "female"),
       # Single-term German clinical vocabulary. Compound free-text values are
       # handled token-wise by load_patient_table(), so only single terms belong
       # here. Extend via `metadata: value_translations:` in the YAML config.
       infection_focus = list(
         Lunge = "lung", pulmonal = "pulmonary", Pneumonie = "pneumonia",
         Abdomen = "abdomen", abdominell = "abdominal", Peritonitis = "peritonitis",
         Cholangitis = "cholangitis", Cholezystitis = "cholecystitis",
         Harnwege = "urinary tract", Harnwegsinfekt = "urinary tract infection",
         Urosepsis = "urosepsis", Niere = "kidney",
         Haut = "skin", Weichteile = "soft tissue", Wunde = "wound",
         Abszess = "abscess", Erysipel = "erysipelas", Fasziitis = "fasciitis",
         Knochen = "bone", Gelenk = "joint", Knie = "knee", "H\u00fcfte" = "hip",
         Osteomyelitis = "osteomyelitis", Spondylodiszitis = "spondylodiscitis",
         "Spondylitis" = "spondylitis", "infizierte H\u00fcft-TEP" = "infected hip prosthesis",
         "H\u00fcft-TEP" = "hip prosthesis", "Knie-TEP" = "knee prosthesis",
         "Protheseninfekt" = "prosthetic joint infection",
         Endokarditis = "endocarditis", Myokarditis = "myocarditis",
         Meningitis = "meningitis", Encephalitis = "encephalitis",
         Katheter = "catheter", ZVK = "central venous catheter",
         Port = "port catheter", ZNS = "central nervous system",
         Weichgewebe = "soft tissue", Mediastinitis = "mediastinitis",
         Sinusitis = "sinusitis", Tonsillitis = "tonsillitis",
         "Bakteri\u00e4mie" = "bacteraemia", Sepsis = "sepsis",
         unbekannt = "unknown", unklar = "unclear", keiner = "none",
         kein = "none", Sonstige = "other", andere = "other"))
}

#' Parse a whole column of dates, resolving day/month order from the column
#'
#' WHY NOT try formats one at a time until one sticks: R's `%Y` happily accepts a
#' two-digit year, so `as.Date("12/11/58", "%d/%m/%Y")` returns 0058-11-12 rather
#' than failing. A trial-and-error loop therefore locks onto the first permissive
#' format and silently produces dates ~1900 years wrong. Two-digit years must be
#' parsed as such, and the field order decided from the column as a whole.
#'
#' WHY column-wise: "12/11/58" is ambiguous in isolation (12 Nov or 11 Dec?), but
#' a column of clinical dates almost always contains at least one value whose
#' first component exceeds 12, which settles the order for every row. Deciding
#' per-row would silently mix conventions within one column.
#'
#' @param v character vector of raw date strings.
#' @param reference_date Date used to place two-digit years and to flag values in
#'   the future.
#' @param label column name, used in messages.
#' @return Date vector, NA where unparseable.
#' @keywords internal
parse_dates_column <- function(v, reference_date = Sys.Date(), label = "date") {
  v <- trimws(as.character(v))
  has <- !is.na(v) & nzchar(v)
  d <- as.Date(rep(NA_real_, length(v)), origin = "1970-01-01")
  if (!any(has)) return(d)

  # ISO first: unambiguous, and must not be reinterpreted as D/M or M/D.
  iso <- has & grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", v)
  if (any(iso)) d[iso] <- as.Date(v[iso], format = "%Y-%m-%d")

  # Numeric triplets separated by . / or -
  tri <- has & is.na(d) & grepl("^[0-9]{1,4}[./-][0-9]{1,2}[./-][0-9]{2,4}$", v)
  if (any(tri)) {
    parts <- do.call(rbind, strsplit(v[tri], "[./-]"))
    p1 <- suppressWarnings(as.integer(parts[, 1]))
    p2 <- suppressWarnings(as.integer(parts[, 2]))
    p3 <- suppressWarnings(as.integer(parts[, 3]))

    # Decide field order from evidence across the whole column: a component
    # greater than 12 can only be a day.
    first_is_day  <- any(p1 > 12, na.rm = TRUE)
    second_is_day <- any(p2 > 12, na.rm = TRUE)
    if (first_is_day && second_is_day)
      warning(label, ": leading and second components both exceed 12 in ",
              "different rows - the column mixes D/M and M/D orders; ",
              "parsed dates are unreliable")
    # Fall back to day-first (European convention) only when nothing decides it.
    day_first <- !second_is_day
    day   <- if (day_first) p1 else p2
    month <- if (day_first) p2 else p1
    year  <- p3

    # Two-digit years: choose the century that does not place the value in the
    # future. A birth date cannot follow the reference date, so "58" read at
    # 2026 means 1958, not 2058.
    two <- nchar(parts[, 3]) <= 2 & !is.na(year)
    if (any(two)) {
      ref_y <- as.integer(format(reference_date, "%Y"))
      cand <- year[two] + (ref_y %/% 100) * 100
      cand[cand > ref_y] <- cand[cand > ref_y] - 100
      year[two] <- cand
    }
    ok <- !is.na(day) & !is.na(month) & !is.na(year) &
          day >= 1 & day <= 31 & month >= 1 & month <= 12
    parsed <- as.Date(rep(NA_real_, length(day)), origin = "1970-01-01")
    parsed[ok] <- as.Date(sprintf("%04d-%02d-%02d", year[ok], month[ok], day[ok]))
    d[tri] <- parsed
    message("[meta] ", label, ": parsed as ",
            if (day_first) "day/month/year" else "month/day/year",
            if (second_is_day || first_is_day) " (decided from the column)"
            else " (no row disambiguates; European default)")
  }

  # Remaining values: month-name formats, both orders.
  todo <- has & is.na(d)
  for (fmt in c("%d %b %Y", "%b %d %Y", "%d %B %Y", "%B %d %Y", "%Y/%m/%d")) {
    if (!any(todo)) break
    d[todo] <- as.Date(v[todo], format = fmt)
    todo <- has & is.na(d)
  }

  bad <- has & is.na(d)
  if (any(bad))
    warning(label, ": ", sum(bad), " unparseable value(s), e.g. ",
            paste(utils::head(unique(v[bad]), 3), collapse = ", "))
  future <- !is.na(d) & d > reference_date
  if (any(future))
    warning(label, ": ", sum(future), " value(s) after the reference date ",
            reference_date, " - check the source column")
  d
}

#' Read, filter, rename and translate the patient table
#'
#' Detects the delimiter and the encoding (German exports are often
#' latin1), trims whitespace, tolerates umlaut spelling variants, converts
#' decimal commas, and PASSES THROUGH any value not in the dictionary while
#' reporting it -- silently mangling an unmapped clinical value would be worse
#' than leaving it in German for the user to map.
#' @param path File path.
#' @param column_map The column map. Default `default_column_map()`.
#' @param value_map The value map. Default `default_value_map()`.
#' @param reference_date The reference date. Default `Sys.Date()`.
#' @export
load_patient_table <- function(path, column_map = default_column_map(),
                               value_map = default_value_map(),
                               reference_date = Sys.Date()) {
  raw <- readLines(path, warn = FALSE, n = 5L)
  delim <- if (mean(lengths(regmatches(raw, gregexpr(";", raw)))) >
                mean(lengths(regmatches(raw, gregexpr(",", raw))))) ";" else ","
  enc <- "UTF-8"
  test <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(test) || any(grepl("\ufffd", test)) ||
      any(!validUTF8(test))) enc <- "latin1"
  df <- read.csv(path, sep = delim, fileEncoding = enc, check.names = FALSE,
                 stringsAsFactors = FALSE, na.strings = c("", "NA", "na", "-", "?"))
  names(df) <- trimws(names(df))
  log_msg("  patient table: ", nrow(df), " rows, ", ncol(df), " columns (",
          delim, "-separated, ", enc, ")")

  out <- data.frame(row.names = seq_len(nrow(df)))
  for (eng in names(column_map)) {
    cands <- column_map[[eng]]
    hit <- NA_integer_
    for (cd in cands) {
      j <- which(tolower(names(df)) == tolower(cd))
      if (length(j)) { hit <- j[1]; break }
    }
    if (is.na(hit)) { warning("column for '", eng, "' not found, skipped"); next }
    out[[eng]] <- df[[hit]]
  }
  if (!ncol(out)) stop("none of the requested columns were found in the patient table")

  # numeric coercion with German decimal commas
  for (nc in intersect(c("age_years", "height_cm", "weight_kg"), names(out))) {
    v <- gsub(",", ".", trimws(as.character(out[[nc]])))
    num <- suppressWarnings(as.numeric(v))
    bad <- !is.na(v) & nzchar(v) & is.na(num)
    if (any(bad)) warning(sum(bad), " unparseable value(s) in ", nc, ": ",
                          paste(unique(v[bad])[1:min(5, sum(bad))], collapse = ", "))
    out[[nc]] <- num
  }
  # ISO dates
  # Physically impossible values are placeholders, not data.
  #
  # WHY: registries pre-create rows at enrolment and fill clinical fields later,
  # leaving numeric columns at 0 rather than blank. A zero height or weight cannot
  # occur in a living patient, so it is unambiguously "not recorded". Left as 0
  # these values are silently averaged, plotted, and used to stratify -- a cohort
  # of "0 kg, age 0" patients that a reader would never notice in a UMAP legend.
  # Converting to NA makes the gap visible instead.
  for (cn in intersect(c("height_cm", "weight_kg"), names(out))) {
    x <- suppressWarnings(as.numeric(out[[cn]]))
    zero <- !is.na(x) & x <= 0
    if (any(zero)) {
      log_msg("  ", cn, ": ", sum(zero),
              " non-positive value(s) treated as missing (placeholder, not data)")
      x[zero] <- NA_real_
    }
    out[[cn]] <- x
  }

  if ("date_of_birth" %in% names(out)) {
    v <- trimws(as.character(out$date_of_birth))
    d <- parse_dates_column(v, reference_date = reference_date,
                            label = "date_of_birth")
    out$date_of_birth <- format(d, "%Y-%m-%d")
    if ("age_years" %in% names(out)) {
      a <- suppressWarnings(as.numeric(out$age_years))
      # An age of 0 is legitimate for a neonate, so it is only rejected when the
      # date of birth is also absent -- i.e. nothing corroborates it. Where a DOB
      # exists it is authoritative and the age is recomputed from it.
      bogus <- !is.na(a) & a <= 0 & is.na(d)
      if (any(bogus)) {
        log_msg("  age_years: ", sum(bogus), " value(s) <= 0 with no date of ",
                "birth treated as missing (placeholder, not a neonate)")
        a[bogus] <- NA_real_
      }
      out$age_years <- a
      derive <- (is.na(out$age_years) | out$age_years <= 0) & !is.na(d)
      if (any(derive)) {
        out$age_years[derive] <- as.numeric(
          floor(as.numeric(difftime(reference_date, d[derive], units = "days")) / 365.25))
        log_msg("  derived age for ", sum(derive), " row(s) from date of birth")
      }
    }
  }
  # value translation
  #
  # Two-pass: exact match first, then TOKEN-WISE for compound entries. Clinical
  # free-text fields routinely hold several findings joined by "/", "+" or ","
  # (e.g. "Pneumonie/Knie/Endokarditis"), and a trailing "?" marks an uncertain
  # diagnosis. Translating token-wise means the dictionary stays a simple list of
  # single terms while still handling arbitrary combinations, and the uncertainty
  # marker is preserved rather than silently dropped.
  untranslated <- list()
  for (cn in intersect(names(value_map), names(out))) {
    dict <- value_map[[cn]]
    keys <- tolower(trimws(names(dict)))
    vals <- unlist(dict, use.names = FALSE)
    v <- trimws(as.character(out[[cn]]))
    newv <- v; miss <- character(0)
    for (i in seq_along(v)) {
      s <- v[i]
      if (is.na(s) || !nzchar(s)) next
      j <- match(tolower(s), keys)
      if (!is.na(j)) { newv[i] <- vals[j]; next }
      uncertain <- grepl("\\?\\s*$", s)
      core <- trimws(sub("\\?\\s*$", "", s))
      toks <- trimws(strsplit(core, "\\s*[/,+]\\s*|\\s+und\\s+")[[1]])
      toks <- toks[nzchar(toks)]
      tj <- match(tolower(toks), keys)
      if (length(toks) && all(!is.na(tj))) {
        newv[i] <- paste0(paste(vals[tj], collapse = " / "),
                          if (uncertain) " (uncertain)" else "")
      } else if (length(toks) > 1L && any(!is.na(tj))) {
        # partial: translate what we can, keep the rest verbatim
        out_tok <- ifelse(is.na(tj), toks, vals[tj])
        newv[i] <- paste0(paste(out_tok, collapse = " / "),
                          if (uncertain) " (uncertain)" else "")
        miss <- c(miss, toks[is.na(tj)])
      } else {
        miss <- c(miss, s)
      }
    }
    if (length(miss)) untranslated[[cn]] <- unique(miss)
    out[[cn]] <- newv
  }
  if (length(untranslated)) {
    for (cn in names(untranslated))
      warning("untranslated value(s) in '", cn, "' (passed through unchanged; ",
              "add them to value_translations in the config): ",
              paste(untranslated[[cn]], collapse = ", "))
  }
  # Report completeness PER COLUMN, not per row. A row-level count is misleading:
  # partially-enrolled registries commonly carry an ID with every clinical field
  # still blank, which a row-level test scores as "has data". Anyone joining this
  # table to cytometry results needs to know which covariates are actually
  # populated before using them for colouring or stratification.
  comp <- vapply(out, function(x) sum(!is.na(x) & trimws(as.character(x)) != ""),
                 integer(1))
  log_msg("  ", nrow(out), " patient rows; non-empty values per column:")
  for (cn in names(comp))
    log_msg(sprintf("    %-18s %4d/%d (%.0f%%)", cn, comp[[cn]], nrow(out),
                    100 * comp[[cn]] / max(1L, nrow(out))))
  sparse <- names(comp)[comp < 0.5 * nrow(out) & names(comp) != "patient_id"]
  if (length(sparse))
    log_msg("  NOTE <50% populated (unreliable for stratification): ",
            paste(sparse, collapse = ", "))
  attr(out, "untranslated") <- untranslated
  out
}

#' Read and validate the sample map
#'
#' Schema (see README): file, well, sample_id, patient_id, timepoint,
#' is_control, panel, fmo_for, control_group. Only `file` is required.
#'
#' `fmo_for` names the markers a file is the fluorescence-minus-one control for,
#' comma separated; `control_group` confines a control to the samples sharing its
#' value. Both are optional and inert when absent. See [parse_fmo_map()].
#' A supplied map that does not cover every input file is a FATAL error: guessing
#' well->patient assignment from plate order would silently mislabel patients.
#' @param path File path.
#' @param fcs_files The fcs files.
#' @export
load_sample_map <- function(path, fcs_files) {
  sm <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                 na.strings = c("", "NA"))
  names(sm) <- trimws(tolower(names(sm)))
  if (!"file" %in% names(sm)) stop("sample map must contain a 'file' column")
  sm$file <- basename(trimws(sm$file))
  if (any(duplicated(sm$file)))
    stop("duplicate file entries in the sample map: ",
         paste(unique(sm$file[duplicated(sm$file)]), collapse = ", "))
  have <- basename(fcs_files)
  unmapped <- setdiff(have, sm$file)
  if (length(unmapped))
    stop("these input files are not in the sample map: ",
         paste(unmapped, collapse = ", "),
         "\n  Add them (or omit --sample-map to use filename-derived labels).")
  if ("is_control" %in% names(sm))
    sm$is_control <- tolower(trimws(as.character(sm$is_control))) %in%
                     c("true", "t", "yes", "y", "1")
  sm
}

#' Write a sample-map template the user can fill in
#' @param fcs_files The fcs files.
#' @param path File path.
#' @param sample_ids The sample ids.
#' @keywords internal
write_sample_map_template <- function(fcs_files, path, sample_ids = NULL) {
  tmpl <- data.frame(
    file = basename(fcs_files),
    well = sample_ids %||% vapply(basename(fcs_files), derive_sample_id, "", kw = NULL),
    sample_id = sample_ids %||% vapply(basename(fcs_files), derive_sample_id, "", kw = NULL),
    patient_id = "", timepoint = "", is_control = "FALSE", panel = "",
    stringsAsFactors = FALSE)
  write.csv(tmpl, path, row.names = FALSE, na = "")
  log_msg("wrote sample map template: ", path,
          ", fill in patient_id/timepoint/is_control and re-run with --sample-map")
  invisible(tmpl)
}

# =============================================================================
