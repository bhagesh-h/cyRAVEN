# SECTION 1 -- FCS READING AND MARKER RESOLUTION
# =============================================================================

#' Do these marker symbols match any ignore pattern
#'
#' Each pattern is matched as a whole name, case-insensitively. A pattern
#' containing `*` is treated as a glob, so `[AF color*` catches all three
#' autofluorescence channels while `CD16` matches CD16 and not CD161. Whole-name
#' matching is the point: substring matching would silently drop markers nobody
#' named.
#'
#' @param sym Character vector of resolved marker symbols.
#' @param patterns Character vector of patterns, or NULL.
#' @return Logical vector the length of `sym`.
#' @keywords internal
ignore_channel_match <- function(sym, patterns) {
  hit <- rep(FALSE, length(sym))
  if (!length(patterns)) return(hit)
  for (p in patterns) {
    p <- trimws(p)
    if (!nzchar(p)) next
    if (grepl("*", p, fixed = TRUE)) {
      esc <- gsub("([][{}()+.^$|\\\\?])", "\\\\\\1", p)
      hit <- hit | grepl(paste0("^", gsub("\\*", ".*", esc), "$"), sym,
                         ignore.case = TRUE)
    } else {
      hit <- hit | (tolower(sym) == tolower(p))
    }
  }
  hit
}

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
#' LOG-AMPLIFIED DATA ARE LINEARISED. An analog instrument with log amplifiers
#' writes integer channel numbers on a log scale and says so in `$PnE`. Every
#' transform this package applies assumes linear intensities, so those channels
#' are converted first, exactly as flowCore's own default does: see
#' [linearise_channels()].
#'
#' MARGIN EVENTS ARE REMOVED. An event at or beyond the top of a channel's
#' range has no true value in that channel, so it is dropped, as flowAI and
#' PeacoQC's `RemoveMargins()` do. See [margin_events()] for which channels are
#' judged, and `margins = "cap"` for the earlier behaviour.
#'
#' @param path File path.
#' @param sample_id The sample id.
#' @param max_events The max events. Default `0L`.
#' @param ignore_channels Character vector of marker names to drop before the
#'   panel fingerprint is taken, or NULL. Names match whole and
#'   case-insensitively; a name containing `*` is treated as a glob.
#' @param margins `"remove"` (default) drops every event at the top of the range
#'   of a marker or scatter channel; `"cap"` keeps the event and sets the value
#'   to the top of the range instead, which is what earlier versions did.
#' @return list(exprs, marker_cols, scatter_cols, channel_names, descs,
#'   keywords, channel_range, margin, n_events, n_events_file, subsampled,
#'   file, sample_id)
#' @export
read_fcs_resolved <- function(path, sample_id = NULL, max_events = 0L,
                              ignore_channels = NULL,
                              margins = c("remove", "cap")) {
  margins <- match.arg(margins)
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
  # truncate_max_range = FALSE, and the clamping done afterwards by
  # clamp_to_declared_range(), which can tell a stray value from a wrong $PnR.
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
  kw   <- flowCore::keyword(ff)
  feat_kw <- vapply(seq_along(nm), function(i)
    as.character(kw[[paste0("$P", i, "FEATURE")]] %||% NA_character_)[1], "")
  rc <- resolve_channels(nm, desc, features = feat_kw)
  marker_cols <- rc$markers; sc_cols <- rc$scatter
  # Pulse height and pulse area are not the same measurement: height
  # understates a bright wide event, so a cut on one is not interchangeable
  # with a cut on the other. Said whenever a marker had to come from height.
  if (length(rc$height_markers))
    log_msg("  pulse HEIGHT used for ", length(rc$height_markers), " marker(s) ",
            "whose detector records no area: ",
            paste(utils::head(rc$height_markers, 8L), collapse = ", "),
            if (length(rc$height_markers) > 8L) ", ..." else "",
            ". Their cuts are not comparable with an area acquisition of the ",
            "same panel")
  stop_if_mass_cytometry(nm, desc, kw, basename(path))

  # WHY DROPPING A CHANNEL HERE MATTERS MORE THAN IT LOOKS. The panel
  # fingerprint is the sorted set of marker names, so a channel present in some
  # files and absent in others splits one cohort into several panels, each with
  # its own cofactor, its own embedding and its own thresholds. Autofluorescence
  # channels extracted during spectral unmixing do exactly this: they are not
  # stains, they vary in number between acquisitions, and left in they can turn
  # twelve comparable files into seven incomparable panels. Removing them before
  # the fingerprint is taken is what makes the cohort one cohort again.
  if (length(ignore_channels)) {
    drop <- ignore_channel_match(names(marker_cols), ignore_channels)
    if (any(drop)) marker_cols <- marker_cols[!drop]
  }

  lin <- linearise_channels(flowCore::exprs(ff), kw, label = basename(path))
  ex <- lin$exprs
  # Judged on the channels the analysis reads, and nowhere else: an event
  # saturating a width or autofluorescence channel nobody uses is a cell.
  used <- unique(c(unname(marker_cols), unname(sc_cols)))
  mg <- margin_events(ex, lin$scale, used, label = basename(path))
  if (identical(margins, "remove")) {
    if (mg$n) ex <- ex[!mg$event, , drop = FALSE]
  } else {
    ex <- clamp_to_declared_range(ex, lin$scale, used, label = basename(path))
  }
  # n_events is what was READ and kept (the denominator of every downstream
  # percentage); n_events_file is what the file claims. They differ under
  # --max-events-per-file and by the margin events removed.
  n_file <- suppressWarnings(as.integer(kw$`$TOT`))
  list(exprs = ex, marker_cols = marker_cols, scatter_cols = sc_cols,
       channel_names = nm, descs = desc, keywords = kw,
       # The top of each channel's range in the units `exprs` now holds: the
       # declared $PnR for a linear channel, f2 * 10^f1 for a linearised log
       # channel, NA where the file's keywords cannot be trusted for it. The
       # logicle transform takes its top of scale from here.
       channel_range = setNames(lin$scale$range, nm),
       channel_scale = lin$scale,
       margin = list(n = if (identical(margins, "remove")) mg$n else 0L,
                     n_at_margin = mg$n, by_channel = mg$by_channel,
                     handling = margins),
       n_events = nrow(ex),
       n_events_file = if (is.finite(n_file)) n_file else n_read_full,
       subsampled = subsampled, file = basename(path),
       sample_id = sample_id %||% derive_sample_id(basename(path), kw))
}

#' Which pulse feature each channel records
#'
#' FCS 3.2 states it in `$PnFEATURE` (Area, Height or Width). Older files say
#' it only through the channel name, by the `-A`/`-H`/`-W` suffix every
#' acquisition package writes, or through a `$PnS` ending in "- Area" as some
#' spectral instruments write it. The keyword wins where it is present.
#'
#' @param nm `$PnN` values, as read (`alter.names` applied).
#' @param desc `$PnS` values.
#' @param features `$PnFEATURE` values, `NA` where absent. Default none.
#' @return character vector of `"A"`, `"H"`, `"W"` or `NA`.
#' @keywords internal
pulse_feature <- function(nm, desc, features = NULL) {
  desc <- as.character(desc); desc[is.na(desc)] <- ""
  f <- rep(NA_character_, length(nm))
  f[grepl("(^|[^A-Za-z])W$|Width$", nm) | grepl("- Width$", desc)] <- "W"
  f[grepl("(^|[^A-Za-z])H$|Height$", nm) | grepl("- Height$", desc)] <- "H"
  f[grepl("(^|[^A-Za-z])A$|Area$", nm) | grepl("- Area$", desc)] <- "A"
  if (!is.null(features)) {
    k <- toupper(substr(trimws(as.character(features)), 1L, 1L))
    use <- !is.na(k) & k %in% c("A", "H", "W")
    f[use] <- k[use]
  }
  f
}

#' The detector a channel belongs to, with its pulse-feature suffix removed
#'
#' `FL1.A` and `FL1.H` are two measurements of one detector, and so is the pair
#' `X..RB613...A` / `X..RB613...H` that `alter.names` makes of a spectral
#' channel. Stripping the suffix gives the key the two share.
#'
#' @param nm `$PnN` values, as read.
#' @return character vector.
#' @keywords internal
detector_of <- function(nm)
  sub("([^A-Za-z])?(A|H|W|Area|Height|Width)$", "", as.character(nm))

#' Decide which channels are markers and which are scatter
#'
#' ONE CHANNEL PER DETECTOR: area where the detector records it, height where
#' it records only height, width never. Area is preferred because it is the
#' integrated signal, and one channel per detector because a marker recorded as
#' both area and height would otherwise be counted twice in the embedding.
#'
#' The choice is made per detector, not per file. A file with area for one
#' detector and height for the others, which is how flowCore's FACSCalibur
#' example and many clinical archives look, used to keep the one area channel
#' and silently drop every height-only marker. FCS 3.2 records the feature per
#' parameter in `$PnFEATURE` for this reason (Spidlen et al. 2021).
#'
#' @param nm `$PnN` values, as read.
#' @param desc `$PnS` values.
#' @param features `$PnFEATURE` values, or NULL.
#' @return list(markers = named integer, marker symbol -> column; scatter =
#'   named integer keyed `FSC-A`, `FSC-H`, `SSC-A`, `SSC-H`; height_markers =
#'   markers taken from pulse height; height_only = TRUE when every marker was)
#' @keywords internal
resolve_channels <- function(nm, desc, features = NULL) {
  nm <- as.character(nm); desc <- as.character(desc)
  sym <- marker_symbol(desc, nm)
  feat <- pulse_feature(nm, desc, features)
  # scatter channels are kept separately for gating, keyed by a canonical name
  scatter <- c(`FSC-A` = "FSC.*A(rea)?$", `FSC-H` = "FSC.*H(eight)?$",
               `SSC-A` = "SSC.*A(rea)?$", `SSC-H` = "SSC.*H(eight)?$")
  sc_cols <- integer(0)
  for (k in names(scatter)) {
    hit <- which(grepl(scatter[[k]], nm) | grepl(scatter[[k]], desc))
    if (length(hit)) sc_cols[k] <- hit[1]
  }
  # The scatter gate is defined on area. Where the file has only height, the
  # height channel stands in for it, under the same caveat as a height marker.
  for (k in c("FSC", "SSC"))
    if (!paste0(k, "-A") %in% names(sc_cols) && paste0(k, "-H") %in% names(sc_cols))
      sc_cols[paste0(k, "-A")] <- sc_cols[[paste0(k, "-H")]]

  fluor <- !grepl("^FSC|^SSC|^Time|Time.Stamp", nm, ignore.case = TRUE) &
           !grepl("^FSC|^SSC|^Time", sym, ignore.case = TRUE) &
           feat %in% c("A", "H")
  det <- detector_of(nm)
  named <- !is.na(desc) & nzchar(trimws(sub(" : .*$", "", desc)))
  keep <- integer(0)
  for (d in unique(det[fluor])) {
    i <- which(fluor & det == d)
    pick <- i[feat[i] == "A"]
    if (!length(pick)) pick <- i[feat[i] == "H"]
    pick <- pick[1L]
    # Both channels of a detector measure the same stain. Where only one of
    # them carries it in $PnS, the chosen channel takes that name rather than
    # falling back to its detector name.
    if (!named[pick] && any(named[i])) sym[pick] <- sym[i[named[i]][1L]]
    keep <- c(keep, pick)
  }
  keep <- sort(keep)
  marker_cols <- setNames(keep, sym[keep])
  # collapse duplicate symbols to the first occurrence
  dup <- duplicated(names(marker_cols))
  marker_cols <- marker_cols[!dup]
  hm <- names(marker_cols)[feat[marker_cols] == "H"]
  list(markers = marker_cols, scatter = sc_cols, height_markers = hm,
       height_only = length(marker_cols) > 0L && length(hm) == length(marker_cols))
}

#' Stop on mass cytometry data, which this package does not process
#'
#' A CyTOF file has no scatter, carries ion counts rather than fluorescence, and
#' needs bead normalisation, a DNA-intercalator singlet gate and an arcsinh
#' cofactor of 5 before analysis (Finck et al. 2013; Nowicka et al. 2017). None
#' of that is implemented here, and the fluorescence path would run on it far
#' enough to produce numbers before failing at the scatter gate. Recognised by
#' FCS 3.2 `$PnTYPE` "Mass", by `$CYT`, or by channel names written as a metal
#' isotope with the `Di`/`Dd` suffix the acquisition software uses.
#'
#' @param nm,desc `$PnN` and `$PnS` values.
#' @param kw keyword list.
#' @param label file name for the message.
#' @return `invisible(FALSE)`, or an error.
#' @keywords internal
stop_if_mass_cytometry <- function(nm, desc, kw, label = "") {
  types <- vapply(seq_along(nm), function(i)
    as.character(kw[[paste0("$P", i, "TYPE")]] %||% "")[1], "")
  cyt <- as.character(kw[["$CYT"]] %||% "")[1]
  # "Yb176Di", and "(Yb176)Dd", which alter.names turns into ".Yb176.Dd"
  iso <- grepl("[A-Z][a-z]?[0-9]{2,3}[^A-Za-z0-9]?(Di|Dd)$", nm)
  if (any(tolower(types) == "mass") ||
      grepl("cytof|helios|hyperion|fluidigm|standard biotools", cyt, ignore.case = TRUE) ||
      sum(iso) >= 3L)
    stop(label, ": this looks like MASS CYTOMETRY data (",
         if (nzchar(cyt)) paste0("$CYT '", cyt, "', ") else "",
         sum(iso), " isotope channel(s)). cyRAVEN analyses fluorescence ",
         "cytometry only: mass cytometry has no scatter and needs bead ",
         "normalisation, a DNA-based singlet gate and an arcsinh cofactor of 5, ",
         "none of which this package applies. See vignette('limitations').",
         call. = FALSE)
  invisible(FALSE)
}

#' Put every channel on a linear scale and find the top of its range
#'
#' WHY. An analog instrument with log amplifiers (a FACSCalibur, for one) stores
#' integer channel numbers 0 to `$PnR` - 1 on a log scale, and `$PnE f1,f2`
#' declares it: f1 decades, offset f2. Every transform this package fits assumes
#' linear intensity, so an arcsinh or logicle laid over channel numbers puts a
#' log on top of a log and moves every valley. FCS 3.1 allows log storage only
#' for integer data (Spidlen et al. 2010).
#'
#' WHAT. Integer channels with f1 > 0 are converted to f2 * 10^(f1 * x / `$PnR`),
#' with f2 = 0 read as 1, which is exactly what flowCore's
#' `read.FCS(transformation = "linearize")` computes. Their range becomes
#' f2 * 10^f1.
#'
#' FLOATING-POINT DATA WITH A NON-ZERO `$PnE` breaks the standard, and the values
#' are not converted, because flowCore does not convert them either and nothing
#' says what scale they are on. flowCore's own GvHD example, written back out
#' with `write.FCS()`, is like this: linear floats up to 10,000 under the
#' FACSCalibur keywords `$PnE 4,0` and `$PnR 1024`. The `$PnR` of such a channel
#' describes the original channel numbers rather than the stored values, so it
#' is not used as the channel's range. Capping those values at 1024 is what put
#' a pile of identical values in the middle of the demo cohort's CD45 data.
#'
#' @param ex Expression matrix, events x channels, in file order.
#' @param kw FCS keywords.
#' @param label File name for the log line.
#' @return list(exprs, scale), `scale` a data.frame with one row per channel:
#'   `f1`, `f2`, `stored` ("linear", "log, linearised" or "log declared on
#'   floats, used as stored"), `range` (top of scale in the units now held, or
#'   NA) and `margin_at` (the value at and above which an event is at the
#'   margin, or NA).
#' @keywords internal
linearise_channels <- function(ex, kw, label = "") {
  n <- ncol(ex)
  dt <- toupper(trimws(as.character(kw[["$DATATYPE"]] %||% "")[1]))
  sc <- data.frame(channel = colnames(ex), f1 = 0, f2 = 0, stored = "linear",
                   range = NA_real_, margin_at = NA_real_,
                   stringsAsFactors = FALSE)
  for (i in seq_len(n)) {
    pnr <- suppressWarnings(as.numeric(kw[[paste0("$P", i, "R")]] %||% NA)[1])
    pne <- suppressWarnings(as.numeric(strsplit(
      as.character(kw[[paste0("$P", i, "E")]] %||% "0,0")[1], ",", fixed = TRUE)[[1]]))
    f1 <- if (length(pne) >= 1L && is.finite(pne[1])) pne[1] else 0
    f2 <- if (length(pne) >= 2L && is.finite(pne[2])) pne[2] else 0
    sc$f1[i] <- f1; sc$f2[i] <- f2
    ok_r <- is.finite(pnr) && pnr > 0
    if (f1 > 0) {
      if (identical(dt, "I") && ok_r) {
        g <- if (f2 == 0) 1 else f2
        ex[, i] <- g * 10^(f1 * ex[, i] / pnr)
        sc$stored[i] <- "log, linearised"
        sc$range[i] <- g * 10^f1
        # The highest channel number, pnr - 1, is the saturated one.
        sc$margin_at[i] <- g * 10^(f1 * (pnr - 1) / pnr)
      } else {
        sc$stored[i] <- "log declared on floats, used as stored"
      }
    } else if (ok_r) {
      sc$range[i] <- pnr
      # Integer data top out at the highest channel number, pnr - 1, and a
      # digitising instrument writing floats saturates at the same value
      # (262,143 under a declared 262,144).
      sc$margin_at[i] <- pnr - 1
    }
  }
  bad <- sc$stored == "log declared on floats, used as stored"
  if (any(bad))
    log_msg("  ", label, ": $PnE declares log amplification on ", sum(bad),
            " floating-point channel(s) (", paste(utils::head(sc$channel[bad], 6L),
                                                  collapse = ", "),
            "), which FCS 3.1 does not allow. The values are used as stored, ",
            "as linear, and their $PnR is not taken as their range. If they are ",
            "log-scaled channel numbers, convert the files before analysis")
  lc <- sc$stored == "log, linearised"
  if (any(lc))
    log_msg("  ", label, ": ", sum(lc), " log-amplified channel(s) linearised ",
            "from $PnE (", paste(utils::head(sc$channel[lc], 6L), collapse = ", "), ")")
  list(exprs = ex, scale = sc)
}

#' Find the events sitting at the top of a channel's range
#'
#' WHAT: an event is at the margin when, in any of the channels `used`, it is at
#' or above that channel's highest value (`margin_at` from
#' [linearise_channels()]).
#'
#' WHY THEY ARE REMOVED RATHER THAN KEPT. A saturated event has no true value in
#' the saturated channel: it is somewhere above the top, by an unknown amount.
#' flowAI removes "margin events in the upper limit of the dynamic range"
#' (Monaco et al. 2016) and PeacoQC's `RemoveMargins()` does the same (Emmaneel
#' et al. 2022). Capping them instead, which flowCore's `read.FCS()` and earlier
#' versions of this package did, keeps the event and invents a value for it: a
#' pile of identical values at the top of the scale. Floating-point files can
#' also carry values far above the declared range, up to 6.8e8 under a declared
#' 262,144 on one cohort, which are not measurements at all.
#'
#' On the 30 spectral files of one sepsis cohort the margin events were 0.01% to
#' 0.7% of each file. A file losing more than 1% is named in a warning, because
#' at that rate the voltage was set too high or the declared range is wrong.
#'
#' @param ex Expression matrix, linear units.
#' @param scale The `scale` table from [linearise_channels()].
#' @param used Column indices judged: the marker and scatter channels.
#' @param warn_frac Warn when more than this fraction of events is at a margin.
#'   Default `0.01`.
#' @param label File name for the log line.
#' @return list(event = logical per event, n = count, by_channel = named
#'   integer count per judged channel with any)
#' @keywords internal
margin_events <- function(ex, scale, used, warn_frac = 0.01, label = "") {
  hit <- rep(FALSE, nrow(ex)); by <- integer(0)
  for (i in intersect(used, seq_len(ncol(ex)))) {
    top <- scale$margin_at[i]
    if (!is.finite(top)) next
    h <- ex[, i] >= top
    h[is.na(h)] <- FALSE
    if (any(h)) {
      by[colnames(ex)[i]] <- sum(h)
      hit <- hit | h
    }
  }
  n <- sum(hit)
  if (n && nrow(ex)) {
    frac <- n / nrow(ex)
    log_msg("  ", label, ": ", n, " margin event(s) (", signif(100 * frac, 2),
            "%) at the top of the range of ", length(by), " channel(s)")
    if (frac > warn_frac)
      log_msg("  WARNING ", label, ": more than ", 100 * warn_frac, "% of events ",
              "are at a channel's margin (", paste(sprintf("%s %d", names(by), by),
                                                    collapse = ", "),
              "). Check the detector voltages and the declared ranges ($PnR)")
  }
  list(event = hit, n = n, by_channel = by)
}

#' Count the events below a channel's lower limit, by flowAI's rule
#'
#' flowAI removes "outliers in the lower limit" of the dynamic range (Monaco et
#' al. 2016): among a channel's negative values, those at or below their median
#' minus 3.5 * `mad()` / 0.6745, which with R's `mad()`, already scaled by
#' 1.4826, is 5.2 robust SDs. PeacoQC's `RemoveMargins()` removes events at the
#' channel's minimum instead (Emmaneel et al. 2022).
#'
#' COUNTED, NOT REMOVED. In compensated data a channel's most negative values
#' are mostly spreading: cells bright in a channel that spills into this one,
#' correctly centred and widely spread. Simulated with 10% of cells bright in a
#' source that widens the receiver five-fold, the rule removed 2% of events, every
#' one of them bright in the source, and took that population from 10.0% to
#' 8.1%. Removing them would bias the frequencies this package reports; the
#' count shows a lower tail the analyst should look at, an electronic baseline
#' fault for instance, without deciding for them.
#'
#' @param ex Expression matrix, after compensation.
#' @param cols Column indices judged: the marker channels. Scatter is not, as
#'   flowAI excludes it by default.
#' @return named integer, events at or below the limit per channel with any,
#'   with attribute `events`, the number of distinct events (one event can be
#'   below the limit in several channels, as flowAI's total_UNIQUE counts).
#' @keywords internal
lower_limit_outliers <- function(ex, cols) {
  by <- integer(0); hit <- rep(FALSE, nrow(ex))
  for (i in intersect(cols, seq_len(ncol(ex)))) {
    x <- ex[, i]
    neg <- x[is.finite(x) & x < 0]
    if (length(neg) < 2L) next
    lim <- stats::median(neg) - 3.5 * stats::mad(neg) / 0.6745
    h <- !is.na(x) & x <= lim
    if (any(h)) {
      by[colnames(ex)[i]] <- sum(h)
      hit <- hit | h
    }
  }
  structure(by, events = sum(hit))
}

#' Cap values at the top of each channel's range
#'
#' The behaviour of earlier versions, kept for `--keep-margin-events`: every
#' value above a channel's range is set to that range, the clamping flowCore's
#' `read.FCS()` does by default. Only channels whose range is known are capped,
#' so a floating-point channel whose `$PnR` describes some other scale (see
#' [linearise_channels()]) is left alone rather than cut in the middle of its
#' data.
#'
#' @param ex Expression matrix, linear units.
#' @param scale The `scale` table from [linearise_channels()].
#' @param used Column indices to cap. Default every column.
#' @param warn_frac Name a channel in a warning when more than this fraction of
#'   its events was capped. Default `0.01`.
#' @param label File name for the log line.
#' @return `ex` with the out-of-range values capped.
#' @keywords internal
clamp_to_declared_range <- function(ex, scale, used = seq_len(ncol(ex)),
                                    warn_frac = 0.01, label = "") {
  n_clamped <- 0L; many <- character(0)
  for (i in intersect(used, seq_len(ncol(ex)))) {
    r <- scale$range[i]
    if (!is.finite(r) || r <= 0) next
    over <- which(ex[, i] > r)
    if (!length(over)) next
    if (length(over) > warn_frac * nrow(ex))
      many <- c(many, sprintf("%s (%.1f%%)", colnames(ex)[i], 100 * length(over) / nrow(ex)))
    ex[over, i] <- r
    n_clamped <- n_clamped + length(over)
  }
  if (n_clamped > 0L)
    log_msg("  ", label, ": ", n_clamped, " value(s) above their channel's range ",
            "capped at it (--keep-margin-events)")
  if (length(many))
    log_msg("  WARNING ", label, ": more than ", 100 * warn_frac, "% of events were above ",
            "the range in ", paste(many, collapse = ", "), ". Capped, but check ",
            "that the declared range is right for these channels")
  ex
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
#' @param optional Marker names that do NOT take part in the fingerprint, so a
#'   reagent present in only some files does not split the cohort into two
#'   panels. Scored where present, UNAVAILABLE where absent. Default none.
#' @return list(assignment = named character (sample_id -> panel), panels = list)
#' @export
fingerprint_panels <- function(reads, labels = NULL, optional = character(0)) {
  # MARKERS THAT DO NOT DEFINE THE PANEL. A reagent stained in some files and
  # not others splits the cohort in two, and each half then gets its own
  # embedding, its own figures and its own file names. That is the right answer
  # when the two halves really are different panels. It is the wrong answer for
  # one added antibody: CD16 on this cohort is in 12 files of 20, and letting it
  # split the run buys three monocyte subsets at the cost of every UMAP position
  # becoming incomparable between the halves, plus a doubled figure set.
  #
  # Naming a marker here keeps it out of the fingerprint, so the files stay one
  # panel and one embedding, while the marker itself stays in `marker_cols` for
  # the files that carry it. A population needing it scores on those files and
  # is reported UNAVAILABLE on the rest, which is the same answer the split gave
  # without the cost.
  fp <- vapply(reads, function(r)
    paste(setdiff(sort(names(r$marker_cols)), optional), collapse = "|"), "")
  uf <- unique(fp)
  pname <- if (!is.null(labels) && length(labels) == length(uf)) labels else
    paste0("panel_", seq_along(uf))
  assignment <- setNames(pname[match(fp, uf)], vapply(reads, `[[`, "", "sample_id"))
  panels <- lapply(seq_along(uf), function(i) {
    idx <- which(fp == uf[i])
    # The panel's marker list is the INTERSECTION over its files, not the first
    # file's. With an optional marker the files in one panel no longer carry
    # identical sets, and taking the first file's would put a marker into the
    # embedding that 8 of 20 files cannot supply.
    .common <- Reduce(intersect, lapply(reads[idx], function(r) names(r$marker_cols)))
    list(name = pname[i], markers = sort(.common),
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
#' Spectral instruments write `$PnS` as `"CD45 : SparkUV-387 - Area"`, the
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
