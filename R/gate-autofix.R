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
# of those minima is a measured cut on this sample's scale too. Under --peer-cuts
# (which --auto-fix-gates implies) peer_cut_replacements() below gives it to a
# sample whose own density showed no minimum; see that function for the rule
# when most samples show none. On 30 sepsis PBMC files, a side-scatter cut at the
# 90th percentile counted only the top tenth of the granulocytes in two samples.
# The replaced cut keeps its needs_review flag, its original value stays in
# derived_threshold, and each replacement is a row of gate_adjustments.csv.

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

#' Population cuts to take from the other samples
#'
#' Two cases, decided per panel and marker by how many samples placed their cut
#' at a density minimum (`valley_*`, `otsu`):
#'
#' * At least `min_peers` and at least `min_peer_frac` of the samples did. The
#'   marker is then bimodal in this assay, and a sample without a minimum -- a
#'   `quantile_fallback` or a `tailgate*` cut -- takes the median of the minima.
#'   On the KMT2C/KMT2D cohort one sample's CD19 tail cut at 3.24, against minima
#'   between 2.28 and 2.50 elsewhere, counted 0.9% B cells where manual gating
#'   counted 5.7%.
#' * Fewer did. The other samples' cuts, fallbacks excluded, are the reference:
#'   a fallback takes their median, and a tail cut more than `z_max` robust SDs
#'   from them (MAD floored at 0.25, as for `scale_outlier`) does too.
#'
#' Cuts from a control, a config entry or a manual override are never touched.
#'
#' @param thr Data frame with `sample_id`, `panel`, `marker`, `threshold`,
#'   `source`, one row per sample and marker.
#' @param z_max Robust z beyond which a tail cut is replaced when the marker has
#'   no majority of density minima.
#' @param min_peers Minimum number of reference cuts.
#' @param min_peer_frac Minimum share of the panel's samples the reference must
#'   cover.
#' @return A data frame of replacements (`sample_id`, `panel`, `marker`,
#'   `old_threshold`, `old_source`, `threshold`, `n_peers`, `robust_z`,
#'   `reason`), or `NULL` when nothing is replaced.
#' @keywords internal
peer_cut_replacements <- function(thr, z_max = 3.5, min_peers = 3L,
                                  min_peer_frac = 0.5) {
  if (is.null(thr) || !nrow(thr)) return(NULL)
  ok <- is.finite(thr$threshold)
  measured <- (grepl("^valley", thr$source) | thr$source == "otsu") & ok
  fallback <- thr$source == "quantile_fallback"
  tail <- startsWith(thr$source, "tailgate")
  key <- paste(thr$panel, thr$marker, sep = "\r")
  out <- list()
  for (k in unique(key)) {
    i <- which(key == k)
    need <- max(min_peers, ceiling(min_peer_frac * length(i)))
    bimodal <- sum(measured[i]) >= need
    for (j in i) {
      if (!fallback[j] && !tail[j]) next
      others <- setdiff(i, j)
      ref <- if (bimodal) thr$threshold[i][measured[i]]
             else thr$threshold[others][ok[others] & !fallback[others]]
      if (length(ref) < max(min_peers, ceiling(min_peer_frac * length(others)))) next
      med <- stats::median(ref)
      z <- if (ok[j]) abs(thr$threshold[j] - med) / max(stats::mad(ref), 0.25) else NA_real_
      if (!bimodal && tail[j] && !(is.finite(z) && z > z_max)) next
      out[[length(out) + 1L]] <- data.frame(
        sample_id = thr$sample_id[j], panel = thr$panel[j], marker = thr$marker[j],
        old_threshold = thr$threshold[j], old_source = thr$source[j], threshold = med,
        n_peers = length(ref), robust_z = round(z, 2),
        reason = if (bimodal)
          sprintf(paste("%d of %d samples found a density minimum for this marker and",
                        "this one did not (%s); replaced by the median of those minima"),
                  sum(measured[i]), length(i), thr$source[j])
        else if (fallback[j])
          sprintf(paste("no density minimum, so the cut was the fallback quantile;",
                        "replaced by the median of the other %d samples' cuts"), length(ref))
        else
          sprintf(paste("tail-rule cut %.2f robust SDs from the other %d samples' cuts;",
                        "replaced by their median"), z, length(ref)),
        stringsAsFactors = FALSE)
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}
