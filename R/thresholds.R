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
#' @param floor no cut is accepted at or below this value. Callers pass the
#'   transformed value of zero signal, `transform$fn(0, marker)`: a stained
#'   population cannot sit below zero. Without it, the small bump that
#'   compensated events form below zero, next to the negative peak, scores as a
#'   deep separation, and the cut lands under the negative peak, so every cell
#'   is called positive. On one cohort that put 87% of CD45+ cells "Vd1+" and
#'   emptied every population defined as Vd1-. Default `-Inf`, no floor.
#' @param ceiling no cut is accepted at or above this value. Default `Inf`; the
#'   viability step passes the dye's dominant mode to look for a minimum below
#'   it (see [apply_gate_hierarchy()]).
#' @param prefer which qualifying valley to return. `"deepest"`, the default,
#'   takes the best-separated pair of modes. `"lowest"` takes the qualifying
#'   valley nearest the bottom of the scale, which is what a gate meant to remove
#'   only NEGATIVE events needs: in whole blood the CD45 distribution has a
#'   negative mode, a dim mode (granulocytes) and a bright one (lymphocytes), the
#'   dim/bright valley is often the deepest, and taking it removed the
#'   granulocytes from the CD45+ parent. Used for the CD45 gate. Under
#'   `"lowest"` a valley also needs `min_lower_frac` of the events below it, so a
#'   ripple in the tail of the negative mode does not count as a separation.
#' @param min_lower_frac see `prefer`. Default `0.01`.
#' @param noise_alpha the family-wise rate at which counting noise alone may
#'   pass for a slope. A valley counts only where the density falls
#'   significantly from the first peak to it and rises significantly from it to
#'   the second, each judged against Poisson noise with a two-sided Bonferroni
#'   quantile over the non-overlapping smoothing windows (a SiZer-style test,
#'   Chaudhuri and Marron 1999). Relative depth alone does not see sample size:
#'   on single modes of 3,000 simulated events it found a valley in 81 of 150
#'   samples, most of them below 95% of the events, and on a flat shoulder of a
#'   few events per bin it split the shoulder. Default `0.05`; `NA` turns the
#'   test off.
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
                           floor = -Inf, prefer = c("deepest", "lowest"),
                           min_lower_frac = 0.01, details = FALSE, noise_alpha = 0.05,
                           ceiling = Inf) {
  prefer <- match.arg(prefer)
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
  # Moving-average smoothing, the window truncated at the two ends rather than
  # the end bins left raw. A raw end bin holding a few events beside an empty
  # one is a "peak" and a "valley" of relative depth 1, which outscores every
  # real separation: on one sample of 2,956 cells it cut Vd2 below every cell
  # and emptied every population defined as Vd2-negative.
  k <- max(1L, as.integer(smooth))
  nb <- length(y); cs <- c(0, cumsum(y))
  lo_b <- pmax(1L, seq_len(nb) - k); hi_b <- pmin(nb, seq_len(nb) + k)
  wb <- hi_b - lo_b + 1L
  ys <- (cs[hi_b + 1L] - cs[lo_b]) / wb
  ctr <- h$mids
  mx <- max(ys)
  pk <- which(c(FALSE, ys[-1] >= ys[-length(ys)]) &
              c(ys[-length(ys)] >= ys[-1], FALSE) & ys > peak_frac * mx)
  if (length(pk) < 2L) return(none("unimodal: fewer than two peaks"))
  min_gap <- min_gap_frac * diff(range(ctr))
  tot <- sum(ys[is.finite(ys)])
  # SIGNIFICANT SLOPE, SiZer-style (Chaudhuri and Marron 1999). At each bin the
  # events in the w bins to its right minus those in the w bins to its left
  # estimate the slope; the counts are Poisson, so that difference has variance
  # equal to their sum. A slope counts as real only beyond a two-sided
  # Bonferroni quantile over the nb / w windows that do not overlap. A dip is
  # then real only where the density falls significantly between the first peak
  # and the valley and rises significantly between the valley and the second:
  # judging the peak and valley heights themselves would not do, since they are
  # the largest and smallest of many noisy bins and differ by several standard
  # errors on a flat stretch with no dip at all.
  # As SiZer does, the slope is judged at several scales, the smoothing window
  # and 2, 4 and 8 times it, since a few thousand events leave too few per
  # window at the finest scale to show even a clear second population; the
  # Bonferroni count covers every window at every scale.
  sig_dn <- sig_up <- rep(FALSE, nb)
  if (is.finite(noise_alpha)) {
    ws <- (2L * k + 1L) * c(1L, 2L, 4L, 8L)
    ws <- ws[nb > 2L * ws]
    zc <- stats::qnorm(1 - noise_alpha / (2 * max(1, sum(nb %/% ws))))
    for (w in ws) {
      js <- w:(nb - w)
      rs <- cs[js + w + 1L] - cs[js + 1L]; ls <- cs[js + 1L] - cs[js - w + 1L]
      zs <- (rs - ls) / sqrt(pmax(rs + ls, 1))
      sig_dn[js] <- sig_dn[js] | zs < -zc; sig_up[js] <- sig_up[js] | zs > zc
    }
  }
  best <- NULL; cands <- list()
  for (i in seq_len(length(pk) - 1L)) for (j in (i + 1L):length(pk)) {
    a <- pk[i]; b <- pk[j]
    if (ctr[b] - ctr[a] < min_gap) next
    vi <- a + which.min(ys[a:b]) - 1L
    if (ctr[vi] <= floor || ctr[vi] >= ceiling) next
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
    # And, when the lowest valley is wanted, so must the lower one.
    if (prefer == "lowest" && sum(ys[seq_len(vi)], na.rm = TRUE) < min_lower_frac * tot) next
    # A dip counting noise could have made; see noise_alpha.
    if (is.finite(noise_alpha) && !(any(sig_dn[a:vi]) && any(sig_up[vi:b]))) next
    # Among candidates that genuinely separate, prefer the deepest separation and
    # break ties toward the wider peak spacing (a better-resolved pair).
    cands[[length(cands) + 1L]] <- list(
      score = rel + 0.05 * (ctr[b] - ctr[a]) / diff(range(ctr)),
      cut = ctr[vi], rel_depth = rel)
  }
  for (cd in cands) {
    better <- if (prefer == "lowest")
      is.null(best) || cd$cut < best$cut || (cd$cut == best$cut && cd$score > best$score)
    else is.null(best) || cd$score > best$score
    if (better) best <- cd
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
#'   3. control quantile               -- unimodal markers with a control present
#'   4. tail gate                      -- unimodal markers without a control;
#'                                        flagged needs_review
#'   5. quantile fallback              -- only where the tail gate cannot be
#'                                        placed; flagged needs_review
#' A valley falling BELOW the control-derived value is rejected: it lies inside
#' the background distribution and would call noise positive.
#'
#' WHY THE TAIL GATE COMES BEFORE THE QUANTILE. A marker with one mode and no
#' control used to be cut at its 90th percentile, which reports exactly 10%
#' positive whatever the sample holds: a constant, not a measurement, and no
#' source supports it. The tail gate places the cut a fixed number of robust
#' standard deviations above the dominant mode, so a sample with no positive
#' cells reports almost none and a sample with a bright tail reports the tail.
#' It answers the case openCyto's `gate_tail` exists for, a single mode with a
#' positive tail (Finak et al. 2014), by a different rule: see
#' [tail_threshold()]. It assumes the dominant mode is the negative one, which a
#' single mode cannot confirm (see [tail_gate_cut()] for why the side is not
#' guessed), so the cut stays flagged for review. The quantile is kept only for
#' the case the tail rule cannot handle, a distribution with no spread to
#' measure.
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
#' @param adaptive Sweep the kernel bandwidth and try the tail and Otsu rules
#'   before falling back to a quantile. See [best_threshold()]. `FALSE` by
#'   default, because it moves thresholds and a threshold that moves moves
#'   every frequency beneath it.
#' @param control_kind What `control_x` is, which decides the `source` string
#'   recorded. `"control_q995"` for an unstained tube, `"fmo_q995"` for a
#'   fluorescence-minus-one control, `"fmm_q995"` for a tube missing several
#'   reagents at once, whose negative is narrower than a true FMO's, so the
#'   cut it anchors is permissive. The arithmetic is identical in all three
#'   cases; they are named apart because they are different experiments and
#'   support different claims. See [parse_fmo_map()].
#' @param prefer Passed to [density_valley()] and [best_threshold()]: `"lowest"`
#'   for a gate that should remove only the negative mode, such as CD45.
#'   Default `"deepest"`.
#' @param floor No derived cut is accepted at or below this value; see
#'   [density_valley()]. Pass `transform$fn(0, marker)` for a fluorescence
#'   marker. Default `-Inf`. A manual, config or control cut is not floored:
#'   each is a decision someone made, not one the density made.
#' @param control_threshold A cut already placed on the control by the caller,
#'   used in place of `control_q` of `control_x`. [control_cut()] produces one
#'   and refuses a control that is brighter than the sample, which a minus
#'   control cannot be. `NULL` by default, in which case the quantile is taken
#'   here as before.
#' @param fmo_sets_cut When `TRUE` and the control is a true FMO
#'   (`control_kind = "fmo_q995"`), the FMO's cut is used even where the
#'   sample's own valley sits higher, and the source reads
#'   `fmo_q995_valley_lowered`. The FMO is the published reference for where
#'   positivity starts (Maecker & Trotter 2006; Cossarizza et al. 2019), so
#'   events between its edge and a higher valley are dim positives, not
#'   negatives. `FALSE` by default: an FMO argues a cut down validly only
#'   where its negative is the sample's own, the same cells stained the same
#'   day, and one FMO tube serves every sample of its panel or `control_group`,
#'   which this package cannot check against each sample; a valley is evidence
#'   from the sample itself. Raising a valley that sits below the FMO's edge is
#'   applied by default, since that errs towards fewer false positives. Never
#'   applied to an unstained or minus-multiple tube: their negatives lack the
#'   panel's spreading and sit too low to argue a cut down.
#' @param tail which side of the dominant mode a tail gate goes, where one is
#'   needed: `"upper"` (default) takes the mode to be the negative population;
#'   `"lower"` is the analyst's declaration that it is the positive one, read
#'   from `thresholds: <marker>: {tail: lower}`. See [tail_gate_cut()].
#' @export
resolve_threshold <- function(marker, x_parent, cfg_value = NULL,
                              control_x = NULL, control_q = 0.995,
                              fallback_q = 0.90, override = NULL,
                              control_kind = c("control_q995", "fmo_q995",
                                               "fmm_q995"),
                              control_threshold = NULL,
                              adaptive = FALSE, floor = -Inf,
                              prefer = c("deepest", "lowest"),
                              fmo_sets_cut = FALSE, tail = c("upper", "lower")) {
  prefer <- match.arg(prefer)
  tail <- match.arg(tail)
  control_kind <- match.arg(control_kind)
  out <- function(threshold, source, needs_review,
                  reason = NA_character_, by = NA_character_,
                  derived = NA_real_, derived_source = NA_character_)
    list(threshold = if (source %in% c("manual", "config")) threshold
                     else untie_cut(x_parent, threshold,
                                    if (identical(source, "quantile_fallback"))
                                      1 - fallback_q else NA_real_),
         source = source, needs_review = needs_review,
         override_reason = reason, override_by = by,
         # Carried so a control-anchored cut can still be compared against the
         # cut the sample's own density implied. Without this the agreement
         # diagnostic compares the control against itself and reports zero.
         derived_threshold = derived, derived_source = derived_source,
         # Kept so threshold_uncertainty() re-derives the cut under the same rule.
         floor = floor, prefer = prefer, tail = tail)

  ov <- suppressWarnings(as.numeric(override$threshold %||% NA_real_))
  if (is.finite(ov))
    return(out(ov, "manual", FALSE,
               reason = as.character(override$reason %||% NA_character_)[1],
               by = as.character(override$set_by %||% NA_character_)[1]))

  if (!is.null(cfg_value) && is.finite(cfg_value))
    return(out(cfg_value, "config", FALSE))
  # A caller that has already placed and validated the cut on the control passes
  # it here; see control_cut(). Everything downstream is identical either way.
  ctrl <- suppressWarnings(as.numeric(control_threshold %||% NA_real_))[1]
  if (!is.finite(ctrl))
    ctrl <- if (!is.null(control_x) && length(control_x) > 100L)
      untie_cut(control_x, as.numeric(quantile(control_x, control_q, na.rm = TRUE)),
                1 - control_q) else NA_real_
  # ADAPTIVE. One fixed smoothing is one hypothesis about how wide the kernel
  # should be, and on a spectral panel it is usually the wrong one: two thirds of
  # the cuts on the cohort this was written for came back as quantile fallbacks
  # because the histogram was fitted to noise. `best_threshold()` sweeps the
  # bandwidth, adds the tail and Otsu rules for the shapes a valley cannot
  # describe, and scores every candidate the same way. Off by default because it
  # moves thresholds, and a threshold that moves moves every frequency under it.
  # WHERE A CONTROL WINS, AND WHERE IT MUST NOT. A control-anchored cut replaces
  # a data-derived one in two cases, and the two are not the same case.
  #
  #   The derived cut has no support. `quantile_fallback` is a constant fraction
  #   of the parent, chosen because nothing in the distribution indicated a cut.
  #   An independent experiment beats that outright, so the control is taken
  #   whether it is higher or lower, and the review flag is cleared because the
  #   cut is no longer blind. This is what a bare `control_q995`, `fmo_q995` or
  #   `fmm_q995` has always meant, and it keeps that name; what changed is that
  #   the adaptive path used to refuse the control here and keep the blind
  #   quantile, which left the weakest cuts in the run unsupported.
  #
  #   The derived cut sits below the control. Here the data did indicate a cut,
  #   and the control says the negative population reaches further than the
  #   sample's own density suggested. Taking the control is the conservative
  #   reading, and the earlier cut is kept on the row so the two can be compared
  #   rather than one of them disappearing.
  #
  # A derived cut ABOVE the control is left alone. The control cannot argue a cut
  # downwards: its negative is at most as wide as the sample's, so a sample whose
  # own density puts the boundary higher has seen something the control cannot.
  # That holds for an unstained or minus-multiple tube, whose negative lacks the
  # panel's spreading. A true FMO's negative does not, and under fmo_sets_cut
  # it is taken as the definition of positivity in both directions.
  take_control <- function(b) {
    if (!is.finite(ctrl)) return(NULL)
    # No density minimum: the quantile and the tail gate both stand in for a
    # cut the data did not show, and an independent control beats either.
    blind <- b$source %in% c("quantile_fallback", "tailgate")
    lower <- isTRUE(fmo_sets_cut) && identical(control_kind, "fmo_q995") &&
      b$threshold > ctrl
    if (!blind && !(b$threshold < ctrl) && !lower) return(NULL)
    out(ctrl, if (blind) control_kind
              else if (lower) paste0(control_kind, "_valley_lowered")
              else paste0(control_kind, "_valley_rejected"),
        FALSE, derived = b$threshold, derived_source = b$source)
  }

  if (isTRUE(adaptive)) {
    b <- best_threshold(x_parent, fallback_q = fallback_q, floor = floor, prefer = prefer,
                        tail = tail)
    if (is.finite(b$threshold)) {
      tc <- take_control(b)
      if (!is.null(tc)) return(tc)
      return(out(b$threshold, b$source,
                 identical(b$source, "quantile_fallback")))
    }
  }

  v <- density_valley(x_parent, floor = floor, prefer = prefer)
  if (is.finite(v)) {
    tc <- take_control(list(threshold = v, source = "valley"))
    if (!is.null(tc)) return(tc)
    return(out(v, "valley", FALSE))
  }
  # No minimum. The tail gate above the dominant mode, where it can be placed
  # sensibly; see tail_gate_cut(). The quantile it used to fall through to
  # reported 10% positive whatever the sample held.
  tg <- if (identical(prefer, "lowest")) NA_real_ else tail_gate_cut(x_parent, floor, tail)
  if (is.finite(tg)) {
    tc <- take_control(list(threshold = tg, source = "tailgate"))
    if (!is.null(tc)) return(tc)
    return(out(tg, "tailgate", TRUE))
  }
  # Floored like every derived cut: a quantile below zero signal would call
  # events positive that no stain put there.
  qf <- untie_cut(x_parent, max(as.numeric(quantile(x_parent, fallback_q, na.rm = TRUE)),
                                floor), 1 - fallback_q)
  tc <- take_control(list(threshold = qf, source = "quantile_fallback"))
  if (!is.null(tc)) return(tc)
  out(qf, "quantile_fallback", TRUE)
}

#' The fixed cut a `thresholds:` configuration entry declares, if any
#'
#' An entry is a bare number, or a list as `write_config()` emits it
#' (`threshold`, `source`, `needs_review`, and optionally `tail`).
#' @param v the entry, or NULL.
#' @return a number, or NULL.
#' @keywords internal
cfg_threshold <- function(v) if (is.list(v)) v$threshold else v

#' The tail-gate side a `thresholds:` configuration entry declares
#'
#' `{tail: lower}` declares that the marker's dominant mode in its parent is
#' the positive population; see [tail_gate_cut()].
#' @param v the entry, or NULL.
#' @return `"lower"` where declared, `"upper"` otherwise.
#' @keywords internal
cfg_tail <- function(v)
  if (is.list(v) && identical(tolower(as.character(v$tail %||% "")), "lower")) "lower" else "upper"

#' Whether a cut was placed without a density minimum
#'
#' `quantile_fallback` and the default path's `tailgate` both mean the sample's
#' density showed no minimum for this marker. Counted together wherever the
#' run reports how often that happened, so the count does not fall just
#' because the tail gate now stands in for the quantile.
#' @param source `source` strings from [resolve_threshold()].
#' @return logical vector.
#' @keywords internal
no_minimum <- function(source) as.character(source) %in% c("quantile_fallback", "tailgate")

#' Move a cut off a block of identical values
#'
#' Populations are scored with strict inequalities, `>` for above and `<` for
#' below, so an event whose value EQUALS the cut is in neither. On continuous
#' data that never happens; on a pile of identical values it decides the
#' result. A quantile of data whose top 15% sit on one value IS that value, and
#' `> cut` then selects nothing: on the demo cohort that emptied every
#' population of 11 samples. Integer channel data, and the capped values of
#' `--keep-margin-events`, are piles of exactly this kind.
#'
#' The cut is moved to the midpoint between the tied value and its nearest
#' neighbour on one side, so the block falls wholly on the other. With a
#' `target_above` (a quantile's intended positive fraction) the side giving the
#' fraction nearer the target is taken; without one the block is kept below the
#' cut, which leaves a density-derived cut where it was to within one step.
#'
#' @param x the values the cut is applied to.
#' @param cut the proposed cut.
#' @param target_above intended fraction of `x` above the cut, or NA.
#' @return the cut, unchanged unless some value of `x` equals it.
#' @keywords internal
untie_cut <- function(x, cut, target_above = NA_real_) {
  if (length(cut) != 1L || !is.finite(cut)) return(cut)
  x <- x[is.finite(x)]
  if (!any(x == cut)) return(cut)
  below <- x[x < cut]; above <- x[x > cut]
  step <- max(abs(cut), 1) * 1e-9
  lo <- if (length(below)) (max(below) + cut) / 2 else cut - step   # block counted above
  hi <- if (length(above)) (cut + min(above)) / 2 else cut + step   # block counted below
  if (!is.finite(target_above)) return(hi)
  if (abs(mean(x > lo) - target_above) <= abs(mean(x > hi) - target_above)) lo else hi
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
