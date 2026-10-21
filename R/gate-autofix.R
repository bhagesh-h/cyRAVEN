# OBSERVATION-DRIVEN CORRECTION OF THE GATE HIERARCHY (--auto-fix-gates)
# =============================================================================
#
# THE OBSERVATION. resolve_threshold() places a cut at a density minimum. Where
# no minimum exists and no control is supplied, it falls back to a fixed
# quantile of the parent, `fallback_q = 0.90`. That number is not a measurement;
# it is an assumption that 10% of the parent is positive. Applied to the two
# gates in the hierarchy it produces a retention fixed by construction, and the
# direction of the comparison decides which artefact you get:
#
#   live_cells   live <- parent & x <  threshold     keeps exactly 90%
#   cd45_pos     cd45 <- live   & x >  threshold     keeps exactly 10%
#
# One constant, two opposite meanings, neither derived from the data. On the
# cohort that surfaced this, `live_cells` was 90.00% of its parent in 12 of 12
# samples -- range 90.00 to 90.00 -- and `cd45_pos` ranged from 10.0% to 98.6%
# depending on whether that sample happened to find a valley.
#
# WHY SKIPPING IS THE FIX, AND NOT A BETTER GUESS. A density minimum is absent
# when the data show one population, not two. In PBMC every cell is a leukocyte,
# so CD45 has no negative mode to separate; a viability dye on a healthy
# preparation has no bright dead mode. Choosing any quantile in that situation
# fabricates a boundary the data do not contain, and the fabricated boundary
# then defines the parent for every population downstream.
#
# Keeping the parent unchanged is the conservative direction for a hierarchy
# gate: it retains cells rather than discarding them, and the populations
# themselves are still discriminated by their own markers, which on this panel
# are valley-derived (CD3, CD4 and CD8 find a minimum in every sample). A gate
# that cannot discriminate should not pretend to.
#
# WHY IT IS OPT-IN. It changes every count in an affected file. That is the same
# rule --drop-unstable-events and --correct-batch follow: an analysis that
# changes a number the previous run reported is asked for, never assumed.
#
# POPULATION MARKERS ARE NOT SKIPPED, THEY BORROW. For a population marker a
# quantile fallback is a poor threshold, but skipping is not the alternative --
# "no threshold" for CD19 does not mean "every cell is a B cell". What the data
# do offer is the same marker in the other samples: every sample of a panel is on
# one pooled transform, so where most of them found a density minimum, the median
# of those minima is a measured cut on this sample's scale too. Under
# --auto-fix-gates peer_cut_replacements() below replaces two kinds of cut with
# it: a quantile fallback, and a tail-rule cut far out of line with the minima.
# On 30 sepsis PBMC files, a side-scatter cut at the 90th percentile counted only
# the top tenth of the granulocytes in two samples, and a CD14 tail cut at twice
# the others' let monocytes in in a third. The replaced cut keeps its
# needs_review flag, its original value stays in derived_threshold, and each
# replacement is a row of gate_adjustments.csv.

#' Decide whether a hierarchy gate rests on evidence or on the fallback constant
#'
#' Separated from `apply_gate_hierarchy()` so the decision is testable on its
#' own and so the reason travels with it as text, rather than being implied by
#' a threshold silently becoming NA.
#'
#' @param source The `source` string `resolve_threshold()` returned.
#' @param enabled Whether `--auto-fix-gates` is set.
#' @return `list(skip, reason)`. `skip = TRUE` means apply no cut at this gate.
#' @keywords internal
gate_needs_autofix <- function(source, enabled = FALSE, kept = NA_real_,
                               parent = NA_real_, min_retention = 0.05) {
  # OBSERVATION 1: the threshold is the fallback constant, not a measurement.
  fabricated <- identical(source, "quantile_fallback")

  # OBSERVATION 2: the gate retains an implausible fraction of its parent,
  # WHATEVER its source said.
  #
  # A density minimum can be found in the wrong place. On the cohort that
  # surfaced this, one acquisition of 20 was the only one where a viability
  # valley was found at all, and that valley kept 140 of 79,000 cells -- 0.18%.
  # A viability gate discarding 99.8% of a PBMC preparation is not measuring
  # viability; it is a misplaced cut, and "valley" as a source does not make it
  # right. Trusting the source alone let it through.
  #
  # The damage is not confined to that sample either. The embedding equalises
  # cells per sample so panel density is comparable between groups, and that
  # equalisation is bounded by the SMALLEST sample -- so one gate collapsing to
  # 140 cells throttled all 20 samples to 140 and cost ~98% of the embedding.
  # A single bad gate is therefore a cohort-wide failure, which is the argument
  # for catching it on its result rather than only on its provenance.
  #
  # 5% is deliberately far below any plausible value rather than near it. PBMC
  # viability runs 80-99%, and a genuinely poor preparation might reach 20-30%;
  # nothing legitimate reaches 5%. The floor is there to catch a broken gate,
  # not to enforce an expectation.
  collapsed <- is.finite(kept) && is.finite(parent) && parent > 0 &&
    (kept / parent) < min_retention

  if (!fabricated && !collapsed)
    return(list(skip = FALSE, reason = NA_character_))

  why <- if (fabricated && collapsed)
    paste0("threshold is a quantile fallback AND the gate retains ",
           sprintf("%.2f%%", 100 * kept / parent), " of its parent")
  else if (fabricated)
    "threshold is a quantile fallback, so retention is fixed by the constant"
  else
    paste0("gate retains ", sprintf("%.2f%%", 100 * kept / parent),
           " of its parent, below the ", sprintf("%.0f%%", 100 * min_retention),
           " floor: the cut is misplaced whatever its source says")

  if (!enabled)
    return(list(skip = FALSE,
                reason = paste0(why, "; --auto-fix-gates would skip it")))
  list(skip = TRUE,
       reason = paste0(why, ". No cut is applied and the parent is carried ",
                       "through unchanged."))
}

#' Record one gate adjustment for gate_adjustments.csv
#' @param sample_id Sample identifier.
#' @param gate Gate name.
#' @param marker Marker the gate is placed on.
#' @param source Threshold source.
#' @param threshold The threshold that would have been applied.
#' @param action "skipped" or "kept".
#' @param reason Free text.
#' @param n_parent Events entering the gate.
#' @param n_kept Events the gate retained.
#' @keywords internal
gate_adjustment_row <- function(sample_id, gate, marker, source, threshold,
                                action, reason, n_parent, n_kept) {
  data.frame(sample_id = sample_id, gate = gate, marker = marker,
             threshold_source = source,
             threshold_not_applied = threshold,
             action = action,
             pct_of_parent_kept = if (is.finite(n_parent) && n_parent > 0)
               round(100 * n_kept / n_parent, 4) else NA_real_,
             reason = reason,
             stringsAsFactors = FALSE)
}

#' Population cuts to take from the other samples' density minima
#'
#' For each panel and marker, the cuts placed at a density minimum (`valley_*`,
#' `otsu`) are the peers. Where at least `min_peers` and at least
#' `min_peer_frac` of the panel's samples have one, two kinds of cut are
#' replaced by the peers' median: a `quantile_fallback`, which is not a
#' measurement, and a `tailgate*` cut more than `z_max` robust SDs from the
#' peers (MAD floored at 0.25, as for `scale_outlier`). Cuts from a control, a
#' config entry or a manual override are never touched.
#'
#' @param thr Data frame with `sample_id`, `panel`, `marker`, `threshold`,
#'   `source`, one row per sample and marker.
#' @param z_max Robust z beyond which a tail-rule cut is replaced.
#' @param min_peers Minimum number of peer minima.
#' @param min_peer_frac Minimum share of the panel's samples with a minimum.
#' @return A data frame of replacements (`sample_id`, `panel`, `marker`,
#'   `old_threshold`, `old_source`, `threshold`, `n_peers`, `robust_z`,
#'   `reason`), or `NULL` when nothing is replaced.
#' @keywords internal
peer_cut_replacements <- function(thr, z_max = 3.5, min_peers = 3L,
                                  min_peer_frac = 0.5) {
  if (is.null(thr) || !nrow(thr)) return(NULL)
  measured <- (grepl("^valley", thr$source) | thr$source == "otsu") &
    is.finite(thr$threshold)
  key <- paste(thr$panel, thr$marker, sep = "\r")
  out <- list()
  for (k in unique(key)) {
    i <- which(key == k)
    pv <- thr$threshold[i][measured[i]]
    if (length(pv) < max(min_peers, ceiling(min_peer_frac * length(i)))) next
    med <- stats::median(pv)
    sc <- max(stats::mad(pv), 0.25)
    for (j in i) {
      src <- thr$source[j]; t0 <- thr$threshold[j]
      z <- if (is.finite(t0)) abs(t0 - med) / sc else NA_real_
      fb <- identical(src, "quantile_fallback")
      far_tail <- startsWith(src %||% "", "tailgate") && is.finite(z) && z > z_max
      if (!fb && !far_tail) next
      out[[length(out) + 1L]] <- data.frame(
        sample_id = thr$sample_id[j], panel = thr$panel[j], marker = thr$marker[j],
        old_threshold = t0, old_source = src, threshold = med,
        n_peers = length(pv), robust_z = round(z, 2),
        reason = if (fb)
          sprintf(paste("no density minimum in this sample, so the cut was the",
                        "fallback quantile; replaced by the median of the %d",
                        "minima found in the other samples"), length(pv))
        else
          sprintf(paste("tail-rule cut %.2f robust SDs from the %d minima found in",
                        "the other samples; replaced by their median"), z, length(pv)),
        stringsAsFactors = FALSE)
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}
