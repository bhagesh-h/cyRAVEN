# SECTION 1 -- FCS READING AND MARKER RESOLUTION
# =============================================================================

#' Read one FCS file and resolve marker symbols
#'
#' WHY marker symbols, not channel names: detector/channel names are
#' instrument-specific and `alter.names = TRUE` mangles them (e.g. `[[RB613]]-A`
#' becomes `X..RB613...A`). The biological identity lives in $PnS, which on this
#' instrument reads "CD45 : SparkUV-387 - Area" -- the symbol is the text before
#' " : ". Building an embedding on channel names (as the original template did)
#' means no marker identity ever reaches the analysis.
#'
#' Note: FCS files larger than 99,999,999 bytes cannot express their data offset
#' in the 8-character HEADER field and write 0 there; the true offsets live in
#' the $BEGINDATA/$ENDDATA keywords. flowCore handles this fallback internally.
#'
#' @param path File path.
#' @param sample_id The sample id.
#' @param max_events The max events. Default `0L`.
#' @return list(exprs, marker_cols, all_cols, keywords, n_events, file, sample_id)
#' @export
read_fcs_resolved <- function(path, sample_id = NULL, max_events = 0L) {
  # Bound the events read per file when asked.
  #
  # WHY EVENLY SPACED rather than the first N: acquisition order is time order,
  # and a flow run drifts -- the first 200k events of an 890k-event file are the
  # first minutes of acquisition, with their own fluidics settling and their own
  # sheath temperature. Thresholds derived from that slice are not thresholds for
  # the file. Taking every k-th event spans the whole run instead. Read the
  # header only first, so the decision costs nothing on small files.
  # WHY read the whole file then subsample in memory, rather than passing the
  # indices to read.FCS: read.FCS's which.lines path seeks to each requested
  # event individually. For ~10^5 scattered indices in a 10^6-event file that is
  # minutes per file -- 14 minutes was observed on a 75 MB file -- while reading the
  # same file contiguously takes seconds. The transient cost is one full matrix in
  # memory (a 10^6 x 21 double matrix is ~170 MB), which is released as soon as the
  # subsample is taken, so peak memory is bounded by the subsample, not the file.
  ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE,
                 alter.names = TRUE)
  n_read_full <- nrow(flowCore::exprs(ff))
  subsampled <- FALSE
  if (max_events > 0L && n_read_full > max_events) {
    keep_rows <- unique(round(seq(1L, n_read_full, length.out = max_events)))
    ff <- ff[keep_rows, ]
    subsampled <- TRUE
  }
  pd <- flowCore::pData(flowCore::parameters(ff))
  desc <- as.character(pd$desc)
  nm   <- as.character(pd$name)
  # marker symbol = text before " : "; fall back to the channel name
  sym <- marker_symbol(desc, nm)

  # AREA channels only. Height/width are redundant duplicates of area and would
  # double-weight every marker in the embedding.
  is_area <- grepl("(^|[^A-Za-z])A$|\\.A$|-A$|Area$", nm) |
             grepl("- Area$", desc)
  # scatter channels are kept separately for gating, keyed by a canonical name
  scatter <- c(`FSC-A` = "FSC.*A(rea)?$", `FSC-H` = "FSC.*H(eight)?$",
               `SSC-A` = "SSC.*A(rea)?$", `SSC-H` = "SSC.*H(eight)?$")
  sc_cols <- integer(0)
  for (k in names(scatter)) {
    hit <- which(grepl(scatter[[k]], nm) | grepl(scatter[[k]], desc))
    if (length(hit)) sc_cols[k] <- hit[1]
  }

  # HEIGHT-ONLY ACQUISITIONS. The area rule above assumes the instrument recorded
  # both, which older cytometers and some clinical archives did not: they write
  # FSC-H, SSC-H and FL1-H with no area channel anywhere. Applying the rule
  # unchanged to such a file resolves ZERO markers, and the run then fails several
  # stages later with an empty population table that blames the specification.
  #
  # The rule exists to avoid double-weighting a marker that appears twice. Where
  # nothing appears twice there is nothing to avoid, so height is used and the
  # substitution is announced. It is announced rather than silent because pulse
  # height and pulse area are not the same measurement: height understates a
  # bright wide event, so thresholds derived here are not interchangeable with
  # those from an area acquisition of the same panel.
  fluor <- !grepl("^FSC|^SSC|^Time|Time.Stamp", nm, ignore.case = TRUE) &
           !grepl("^FSC|^SSC|^Time", sym, ignore.case = TRUE)
  height_only <- !any(is_area & fluor) && any(fluor)
  if (height_only) {
    is_h <- grepl("(^|[^A-Za-z])H$|\\.H$|-H$|Height$", nm) |
            grepl("- Height$", desc)
    if (any(is_h & fluor)) {
      is_area <- is_area | is_h
      log_msg("  no area channel in this file; using pulse HEIGHT for ",
              sum(is_h & fluor), " marker(s). Height and area are different ",
              "measurements, so thresholds from this file are not comparable ",
              "with an area acquisition of the same panel")
    }
  }
  # The scatter gate is defined on area. Where the file has only height, the
  # height channel stands in for it, under the same caveat.
  for (k in c("FSC", "SSC"))
    if (!paste0(k, "-A") %in% names(sc_cols) && paste0(k, "-H") %in% names(sc_cols))
      sc_cols[paste0(k, "-A")] <- sc_cols[[paste0(k, "-H")]]

  # fluorescence markers = area channels that are not scatter/time
  keep <- which(is_area & fluor)
  marker_cols <- setNames(keep, sym[keep])
  # collapse duplicate symbols to the first occurrence
  marker_cols <- marker_cols[!duplicated(names(marker_cols))]

  ex <- flowCore::exprs(ff)
  # n_events is what was READ (the denominator of every downstream percentage);
  # n_events_file is what the file claims. They differ only under --max-events-per-file.
  n_file <- suppressWarnings(as.integer(flowCore::keyword(ff)$`$TOT`))
  list(exprs = ex, marker_cols = marker_cols, scatter_cols = sc_cols,
       channel_names = nm, descs = desc, keywords = flowCore::keyword(ff),
       n_events = nrow(ex),
       n_events_file = if (is.finite(n_file)) n_file else n_read_full,
       subsampled = subsampled, file = basename(path),
       sample_id = sample_id %||% derive_sample_id(basename(path), flowCore::keyword(ff)))
}

#' Derive a sample label from the FCS keywords or filename
#' WHY: wells rarely carry patient IDs; this gives a stable plot label when no
#' sample map is supplied. $WELLID is preferred, then $SMNO, then the filename.
#' @param fname The fname.
#' @param kw The kw.
#' @keywords internal
derive_sample_id <- function(fname, kw = NULL) {
  for (k in c("$WELLID", "$SMNO", "WELL ID", "$SRC")) {
    v <- kw[[k]]
    if (!is.null(v) && nzchar(trimws(v))) return(trimws(v))
  }
  # fall back to a leading well-like token, else the filename stem
  m <- regmatches(fname, regexpr("^[A-H][0-9]{1,2}", fname))
  if (length(m) && nzchar(m)) return(m)
  sub("\\.fcs$", "", fname, ignore.case = TRUE)
}

#' Group files into panels by their exact marker set
#'
#' WHY: markers define the feature space. Two files with different marker sets
#' have no common space, so binding them and embedding gives meaningless
#' coordinates -- the original template's blind row-bind produced NA columns for
#' every non-shared channel and then aborted on its own finiteness check.
#'
#' @param reads Named list of objects returned by [read_fcs_resolved()].
#' @param labels Vector of labels, one per row of coords.
#' @return list(assignment = named character (sample_id -> panel), panels = list)
#' @export
fingerprint_panels <- function(reads, labels = NULL) {
  fp <- vapply(reads, function(r) paste(sort(names(r$marker_cols)), collapse = "|"), "")
  uf <- unique(fp)
  pname <- if (!is.null(labels) && length(labels) == length(uf)) labels else
    paste0("panel_", seq_along(uf))
  assignment <- setNames(pname[match(fp, uf)], vapply(reads, `[[`, "", "sample_id"))
  panels <- lapply(seq_along(uf), function(i) {
    idx <- which(fp == uf[i])
    list(name = pname[i], markers = sort(names(reads[[idx[1]]]$marker_cols)),
         samples = vapply(reads[idx], `[[`, "", "sample_id"), n_files = length(idx))
  })
  names(panels) <- pname
  if (length(panels) > 1L) {
    log_msg(length(panels), " distinct panels detected, each gets its own embedding")
    all_m <- unique(unlist(lapply(panels, `[[`, "markers")))
    for (p in panels) {
      miss <- setdiff(all_m, p$markers)
      log_msg("  ", p$name, ": ", p$n_files, " file(s), ", length(p$markers),
              " markers", if (length(miss)) paste0("; lacks: ", paste(miss, collapse = ", ")) else "")
    }
  }
  list(assignment = assignment, panels = panels)
}

# =============================================================================

#' The marker name cyRAVEN resolves from an FCS description
#'
#' Spectral instruments write `$PnS` as `"CD45 : SparkUV-387 - Area"` — the
#' antibody, the detector, and the pulse statistic in one string. The name a
#' population specification has to use is the antibody alone, so everything from
#' `" : "` onward is stripped. Where `$PnS` is empty the channel name `$PnN`
#' stands in.
#'
#' EXTRACTED SO IT CANNOT DIVERGE. [read_fcs_resolved()] and
#' [report_input_check()] both need it, and for a while only the first had it:
#' `--check` then reported markers as the full `$PnS` string and told the user
#' that every population named a marker no file contained, on data where the run
#' itself would have resolved all of them. A false alarm from the one command
#' whose whole purpose is to catch real ones.
#'
#' @param desc `$PnS` values.
#' @param nm `$PnN` values, used where `desc` is empty.
#' @return Character vector of resolved marker names.
#' @keywords internal
marker_symbol <- function(desc, nm = NULL) {
  s <- trimws(sub(" : .*$", "", as.character(desc)))
  if (!is.null(nm)) {
    bad <- is.na(s) | !nzchar(s)
    s[bad] <- as.character(nm)[bad]
  }
  s
}
