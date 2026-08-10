# SECTION 3 -- DENSITY VALLEY AND THRESHOLD RESOLUTION
# =============================================================================

#' Find the deepest density valley between two modes
#'
#' WHAT: smooths a histogram, locates peaks above `peak_frac` of the maximum,
#' and returns the valley between the peak pair maximising depth * min(height).
#' Returns NA when the distribution is unimodal -- the honest answer, which then
#' routes threshold resolution to the unstained-control fallback.
#' WHY: for a bimodal marker the valley is the least-arbitrary cutoff available
#' and adapts to each sample's own staining intensity.
#' @param min_rel_depth minimum valley depth AS A FRACTION of the smaller of the
#'   two flanking peaks. This is the parameter that makes minority positive
#'   populations findable -- see the note on scoring below.
#' @param min_upper_frac the upper mode must hold at least this fraction of events,
#'   so a cut cannot be placed out in the extreme tail where a handful of events
#'   form a bump.
#' @param x A vector of values.
#' @param bins The bins. Default `220L`.
#' @param smooth The smooth. Default `4`.
#' @param peak_frac The peak frac. Default `0.02`.
#' @param min_gap_frac The min gap frac. Default `0.06`.
#' @param range_q The range q. Default `c(0.001, 0.999)`.
#' @param details when TRUE return a list carrying the cut plus the relative
#'   depth of the valley and why no cut was found, instead of the bare cut.
#'   The default returns exactly what it always has: one number or `NA_real_`.
#'   [threshold_uncertainty()] uses the list form, because how deep a valley is
#'   determines how far the cut can move under resampling, and recomputing the
#'   smoothed histogram elsewhere would let the two drift apart.
#' @export
density_valley <- function(x, bins = 220L, smooth = 4, peak_frac = 0.02,
                           min_gap_frac = 0.06, range_q = c(0.001, 0.999),
                           min_rel_depth = 0.30, min_upper_frac = 0.002,
                           details = FALSE) {
  none <- function(reason)
    if (details) list(cut = NA_real_, rel_depth = NA_real_, n_peaks = NA_integer_,
                      reason = reason) else NA_real_
  x <- x[is.finite(x)]
  if (length(x) < 500L) return(none("fewer than 500 finite values"))
  rg <- quantile(x, range_q, na.rm = TRUE)
  if (!all(is.finite(rg)) || diff(rg) <= 0) return(none("degenerate range"))
  h <- graphics::hist(x[x >= rg[1] & x <= rg[2]], breaks = seq(rg[1], rg[2], length.out = bins + 1L),
            plot = FALSE)
  y <- h$counts
  # moving-average smoothing
  k <- max(1L, as.integer(smooth))
  ys <- stats::filter(y, rep(1 / (2 * k + 1), 2 * k + 1), sides = 2)
  ys[is.na(ys)] <- y[is.na(ys)]
  ctr <- h$mids
  mx <- max(ys)
  pk <- which(c(FALSE, ys[-1] >= ys[-length(ys)]) &
              c(ys[-length(ys)] >= ys[-1], FALSE) & ys > peak_frac * mx)
  if (length(pk) < 2L) return(none("unimodal: fewer than two peaks"))
  min_gap <- min_gap_frac * diff(range(ctr))
  tot <- sum(ys[is.finite(ys)])
  best <- NULL
  for (i in seq_len(length(pk) - 1L)) for (j in (i + 1L):length(pk)) {
    a <- pk[i]; b <- pk[j]
    if (ctr[b] - ctr[a] < min_gap) next
    vi <- a + which.min(ys[a:b]) - 1L
    lo <- min(ys[a], ys[b])
    depth <- lo - ys[vi]
    if (depth <= 0) next
    # RELATIVE depth: how completely the two modes separate, independent of how
    # tall they are. The previous scoring used depth * min(height), which is an
    # ABSOLUTE quantity and therefore preferred whichever peak pair contained the
    # most events. That is wrong whenever the positive population is a minority:
    # in a sample where CD45+ cells are ~5% of events, a shallow ripple inside the
    # enormous negative peak scores higher than the clean, deep separation between
    # negative and positive, so the threshold lands INSIDE the negative
    # distribution and calls background positive. Because CD45 sits at the root of
    # the gate hierarchy, that one misplaced cut then corrupts every population
    # beneath it, and it does so silently -- the valley is real, just the wrong one.
    rel <- depth / lo
    if (rel < min_rel_depth) next
    # The upper mode must be a population, not a tail artefact.
    if (sum(ys[vi:length(ys)], na.rm = TRUE) < min_upper_frac * tot) next
    # Among candidates that genuinely separate, prefer the deepest separation and
    # break ties toward the wider peak spacing (a better-resolved pair).
    score <- rel + 0.05 * (ctr[b] - ctr[a]) / diff(range(ctr))
    if (is.null(best) || score > best$score)
      best <- list(score = score, cut = ctr[vi], rel_depth = rel)
  }
  if (is.null(best)) return(none("no peak pair separates enough"))
  if (details) list(cut = best$cut, rel_depth = best$rel_depth,
                    n_peaks = length(pk), reason = NA_character_) else best$cut
}

#' Resolve a threshold for one marker, recording where the value came from
#'
#' Priority (per the agreed contract):
#'   0. per-sample manual override     -- wins over everything, and is recorded
#'   1. explicit config value          -- applies to every sample
#'   2. density valley within parent   -- bimodal markers
#'   3. unstained-control quantile     -- unimodal markers with a control present
#'   4. quantile fallback              -- flagged needs_review
#' A valley falling BELOW the control-derived value is rejected: it lies inside
#' the background distribution and would call noise positive.
#'
#' WHY LEVEL 0 EXISTS AND WHY IT IS DISTINCT FROM LEVEL 1. When this package
#' flags a threshold for review, the only previous responses were to accept it or
#' to pin that marker across the whole run through `thresholds:`. Pinning is the
#' worse of the two: it applies one number to every sample and so reintroduces
#' exactly the fixed-coordinate bias described in README section 1, in order to
#' correct one tube. A per-sample override corrects the tube.
#'
#' It is a separate `source` value rather than reusing `config` because the two
#' are different claims. `config` says the assay declares this cut; `manual` says
#' a named person moved this one sample's cut for a stated reason on a stated
#' run. An override that is recorded and attributed is more defensible than an
#' automated cut nobody was permitted to correct, and less defensible than one
#' nobody needed to.
#'
#' @param marker Marker name.
#' @param x_parent The x parent.
#' @param cfg_value The cfg value.
#' @param control_x The control x.
#' @param control_q The control q. Default `0.995`.
#' @param fallback_q The fallback q. Default `0.9`.
#' @param override Optional list for THIS sample and marker carrying
#'   `threshold`, and optionally `reason` and `set_by`. See
#'   [sample_override()].
#' @param control_kind What `control_x` is, which decides the `source` string
#'   recorded. `"control_q995"` for an unstained tube, `"fmo_q995"` for a
#'   fluorescence-minus-one control. The arithmetic is identical; the two are
#'   named apart because they are different experiments and support different
#'   claims. See [parse_fmo_map()].
#' @export
resolve_threshold <- function(marker, x_parent, cfg_value = NULL,
                              control_x = NULL, control_q = 0.995,
                              fallback_q = 0.90, override = NULL,
                              control_kind = c("control_q995", "fmo_q995")) {
  control_kind <- match.arg(control_kind)
  out <- function(threshold, source, needs_review,
                  reason = NA_character_, by = NA_character_)
    list(threshold = threshold, source = source, needs_review = needs_review,
         override_reason = reason, override_by = by)

  ov <- suppressWarnings(as.numeric(override$threshold %||% NA_real_))
  if (is.finite(ov))
    return(out(ov, "manual", FALSE,
               reason = as.character(override$reason %||% NA_character_)[1],
               by = as.character(override$set_by %||% NA_character_)[1]))

  if (!is.null(cfg_value) && is.finite(cfg_value))
    return(out(cfg_value, "config", FALSE))
  ctrl <- if (!is.null(control_x) && length(control_x) > 100L)
    as.numeric(quantile(control_x, control_q, na.rm = TRUE)) else NA_real_
  v <- density_valley(x_parent)
  if (is.finite(v)) {
    if (is.finite(ctrl) && v < ctrl)
      return(out(ctrl, paste0(control_kind, "_valley_rejected"), FALSE))
    return(out(v, "valley", FALSE))
  }
  if (is.finite(ctrl))
    return(out(ctrl, control_kind, FALSE))
  out(as.numeric(quantile(x_parent, fallback_q, na.rm = TRUE)),
      "quantile_fallback", TRUE)
}

#' Look up a per-sample, per-marker threshold override
#'
#' Reads the `sample_overrides:` block of the config, which is keyed on sample
#' identifier and then on marker:
#'
#' ```yaml
#' sample_overrides:
#'   D07:
#'     CCR7:
#'       threshold: 2.15
#'       reason: "valley found inside the negative mode, see gating_qc.png"
#'       set_by: "initials or name"
#' ```
#'
#' Sample and marker names are matched exactly, as they appear in
#' `thresholds_used.csv`. A block naming a sample or marker that does not exist
#' is inert rather than an error, because a config is routinely reused across
#' cohorts that do not all contain the same tubes; unmatched entries are reported
#' once by [report_unused_overrides()] rather than failing the run.
#'
#' @param overrides the `sample_overrides` list, or NULL
#' @param sample_id sample identifier
#' @param marker marker name
#' @return a list carrying at least `threshold`, or NULL when none applies
#' @export
sample_override <- function(overrides, sample_id, marker) {
  if (is.null(overrides) || !length(overrides)) return(NULL)
  s <- overrides[[as.character(sample_id)]]
  if (is.null(s) || !length(s)) return(NULL)
  e <- s[[as.character(marker)]]
  if (is.null(e)) return(NULL)
  # A bare number is accepted as shorthand, but it carries no reason and no
  # attribution, which is most of the point of the feature.
  if (!is.list(e)) e <- list(threshold = e)
  if (!is.finite(suppressWarnings(as.numeric(e$threshold %||% NA_real_)))) return(NULL)
  e
}

#' Report override entries that matched nothing
#'
#' An override silently applying to no sample is the failure mode this guards:
#' the analyst believes a cut was corrected, the run says nothing, and the
#' uncorrected number is published.
#'
#' @param overrides the `sample_overrides` list, or NULL
#' @param applied character vector of "sample\\rmarker" keys that were used
#' @return invisibly, the character vector of unused keys
#' @keywords internal
report_unused_overrides <- function(overrides, applied) {
  if (is.null(overrides) || !length(overrides)) return(invisible(character(0)))
  want <- unlist(lapply(names(overrides), function(s)
    paste0(s, "\r", names(overrides[[s]]))), use.names = FALSE)
  unused <- setdiff(want, applied)
  if (length(unused)) {
    log_msg("WARNING ", length(unused), " sample_overrides entr(y/ies) matched no ",
            "sample and marker in this run and were NOT applied:")
    for (u in utils::head(unused, 8L))
      log_msg("    ", sub("\r", " / ", u))
    if (length(unused) > 8L) log_msg("    ... and ", length(unused) - 8L, " more")
    log_msg("  names are matched exactly as they appear in thresholds_used.csv")
  }
  invisible(unused)
}

# =============================================================================
