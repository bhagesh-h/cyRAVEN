# SECTION 3b -- GATE PLACEMENT UNCERTAINTY
# =============================================================================
#
# WHY THIS FILE EXISTS. Every frequency this package reports is a count of the
# cells on one side of a cut, and that cut was estimated from the same cells it
# divides. Reported on its own, the frequency claims a precision it does not
# have. A population separated by a deep, clean valley barely moves when the cut
# is displaced; one sitting on a shallow shoulder moves a great deal. Both are
# printed to the same number of decimal places, and nothing downstream can tell
# them apart.
#
# The consequence is not academic. Operator studies that follow the GUM
# convention report expanded uncertainty rising from about 12% on a three-gate
# strategy to about 16% on a five-gate one, with most of it entering at the first
# gate (Cadwell et al. 2021; Whitmore et al. 2021, Methods Protoc 4:24). Those
# figures describe humans placing gates. This file measures the same quantity for
# the placement rule the package uses instead, which is the only version of the
# number a reader of these outputs can act on.
#
# TWO COMPONENTS, COMBINED IN QUADRATURE.
#
#   sampling -- resample the parent-gate events with replacement and re-derive
#     the cut. This is the uncertainty from having counted a finite number of
#     cells, and it shrinks as events accumulate.
#
#   method -- re-derive the cut over a grid of the settings density_valley()
#     itself takes: histogram resolution, smoothing width, and how completely two
#     modes must separate before the gap between them counts. These are analyst
#     choices that this package makes on the analyst's behalf. Publishing the
#     spread they produce is the honest form of the claim that automation removed
#     the subjectivity; it did not, it fixed it, and this is how much was fixed.
#
# WHAT IT DELIBERATELY DOES NOT DO. It never changes the reported threshold. The
# cut written to thresholds_used.csv is still density_valley() at its defaults on
# the real events. Everything here perturbs a copy. A run with uncertainty
# enabled and a run without produce identical frequencies, MFIs and p-values.
#
# RNG. Every entry point saves and restores .Random.seed, for the reason spelled
# out at run_unsupervised_clusters(): run_cyraven() seeds once and the UMAP cell
# selection draws from that one stream, so a step that consumes draws without
# putting the stream back changes which cells are embedded and silently redraws
# every UMAP figure in the run.

#' The settings sweep that produces the method component
#'
#' Three factors, each at the package default and one step either side. Bins and
#' smoothing set how much structure the histogram retains; `min_rel_depth` sets
#' how completely two modes must separate before the gap between them is accepted
#' as a threshold, which is the setting that decides whether a minority positive
#' population is found at all.
#'
#' @return a data.frame of setting combinations, one per row
#' @export
valley_setting_grid <- function() {
  expand.grid(bins = c(160L, 220L, 300L), smooth = c(3, 4, 6),
              min_rel_depth = c(0.22, 0.30, 0.38),
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

#' Uncertainty in one per-sample threshold
#'
#' Returns the standard uncertainty of the cut, in the units of `x`, split into
#' the component from finite sampling and the component from the placement
#' settings, plus their quadrature sum.
#'
#' HOW EACH SOURCE IS TREATED. `source` is the string
#' [resolve_threshold()] recorded for this threshold, and it determines what
#' "uncertain" even means:
#'
#'   - `config`: the value was declared, not estimated. Uncertainty is zero by
#'     construction, and saying so is different from saying it is small.
#'   - `valley`: both components are computed as described at the top of this
#'     file.
#'   - `quantile_fallback`: no valley was found, so there is no valley to
#'     resample. Sampling uncertainty is the bootstrap spread of the quantile.
#'     The method component sweeps which quantile, because the choice of 0.90 is
#'     arbitrary and the spread it produces is normally large. That is the
#'     correct answer: a fallback threshold is barely determined by the data, and
#'     this is the number that says so.
#'   - `control_q995` and `control_q995_valley_rejected`: the cut came from a
#'     separate control tube not passed here, so neither component can be
#'     computed from `x`. Both are NA with a reason, rather than a small number
#'     that would read as confidence.
#'
#' @param x parent-gate values for one marker in one sample, on the analysis
#'   scale
#' @param source derivation string from [resolve_threshold()]
#' @param B bootstrap replicates
#' @param seed seed for the local RNG stream
#' @param max_events cap on events per bootstrap replicate; the cut is a
#'   histogram feature and stops moving long before the full parent gate is used
#' @param fallback_q the quantile [resolve_threshold()] falls back to
#' @param grid setting combinations for the method component
#' @return list(u_sampling, u_method, u_combined, rel_depth, valley_rate,
#'   n_events, basis)
#' @export
threshold_uncertainty <- function(x, source = "valley", B = 100L, seed = 42L,
                                  max_events = 20000L, fallback_q = 0.90,
                                  grid = valley_setting_grid()) {
  out <- list(u_sampling = NA_real_, u_method = NA_real_, u_combined = NA_real_,
              rel_depth = NA_real_, valley_rate = NA_real_,
              n_events = 0L, basis = NA_character_)
  x <- x[is.finite(x)]
  out$n_events <- length(x)

  if (identical(source, "config")) {
    out$u_sampling <- 0; out$u_method <- 0; out$u_combined <- 0
    out$basis <- "declared"
    return(out)
  }
  # Both control paths are outside this function's reach: the cut came from a
  # separate tube whose events were not passed in, so neither component can be
  # estimated from `x`. NA with a reason, rather than a small number that would
  # read as confidence.
  if (grepl("^control_q995|^fmo_q995", source %||% "")) {
    out$basis <- if (grepl("^fmo", source)) "FMO control not available here"
                 else "control tube not available here"
    return(out)
  }
  if (length(x) < 500L) {
    out$basis <- "fewer than 500 parent events"
    return(out)
  }

  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  n_draw <- min(length(x), as.integer(max_events))
  is_fb  <- identical(source, "quantile_fallback")

  boot <- vapply(seq_len(as.integer(B)), function(i) {
    xb <- x[sample.int(length(x), n_draw, replace = TRUE)]
    if (is_fb) as.numeric(quantile(xb, fallback_q, na.rm = TRUE))
    else density_valley(xb)
  }, numeric(1))
  ok <- is.finite(boot)
  out$valley_rate <- mean(ok)
  # Fewer than a third of replicates finding a cut is not a small uncertainty,
  # it is an unstable threshold, and an SD over the minority that did find one
  # would understate it by conditioning on success.
  out$u_sampling <- if (sum(ok) >= max(10L, 0.33 * B)) stats::sd(boot[ok]) else NA_real_

  xs <- if (length(x) > max_events) x[sample.int(length(x), max_events)] else x
  if (is_fb) {
    qs <- c(0.85, 0.90, 0.95)
    meth <- vapply(qs, function(q) as.numeric(quantile(xs, q, na.rm = TRUE)),
                   numeric(1))
    out$basis <- "quantile bootstrap; sweep over the fallback quantile"
  } else {
    meth <- vapply(seq_len(nrow(grid)), function(i)
      density_valley(xs, bins = grid$bins[i], smooth = grid$smooth[i],
                     min_rel_depth = grid$min_rel_depth[i]), numeric(1))
    d <- density_valley(xs, details = TRUE)
    out$rel_depth <- d$rel_depth
    out$basis <- "valley bootstrap; sweep over bins, smoothing and depth"
  }
  mok <- is.finite(meth)
  out$u_method <- if (sum(mok) >= 3L) stats::sd(meth[mok]) else NA_real_

  cmp <- c(out$u_sampling, out$u_method)
  out$u_combined <- if (any(is.finite(cmp)))
    sqrt(sum(cmp[is.finite(cmp)]^2)) else NA_real_
  out
}

# =============================================================================
# SECTION 3c -- COUNTING UNCERTAINTY AND DETECTION LIMITS
# =============================================================================
#
# WHY THIS IS SEPARATE FROM EVERYTHING ABOVE. The section above answers how far a
# frequency moves when the cut behind it moves. It does not answer how many cells
# there were to count, and only one of those two quantities was being reported. A
# population of twelve events and a population of twelve thousand, both sitting
# behind a deep clean valley, came out with the same small uncertainty, because
# displacing a well-separated cut moves neither of them much. One of those
# numbers is worth acting on and the other is not, and nothing in the output
# distinguished them.
#
# This is the oldest quantified source of imprecision in cytometry and the reason
# rare-event protocols specify event targets at all: about 4,000 events of a
# population gives a 2% coefficient of variation on its frequency, and 25 events
# gives 20%. Clinical practice fixes the two limits that follow. Roughly 20
# events is the fewest that can be called a detection and roughly 50 the fewest
# that can be called a measurement, against a stated denominator (Sommer et al.
# 2021, Cytometry B 100:42; the GIMEMA AML1310 post-hoc analysis, Haematologica
# 2022, 107:2823).
#
# WHY WILSON RATHER THAN THE TEXTBOOK FORMULA. The binomial standard error
# sqrt(p(1-p)/n) tends to zero as p does. At k = 0 it reports perfect certainty
# about a population nobody observed, which is the reverse of the truth, and at
# single-digit k it is still badly optimistic -- exactly the regime where the
# number is being consulted. The Wilson score interval has no such failure at the
# boundary. Its half-width at z = 1 is taken here as the standard uncertainty,
# and it agrees with the textbook value to three figures once a population has a
# few hundred events, so nothing is lost in the common case.
#
# THE DENOMINATOR IS THIS RUN'S, NOT THE ACQUISITION'S. n is the number of
# parent-gate events this analysis saw, so --max-events-per-file lowers it and
# raises every limit reported here in proportion. That is the honest reading: a
# population under the limit of detection of a subsample may be perfectly well
# resolved in the whole file, and the answer is to raise the cap rather than to
# believe the limit.

#' Counting uncertainty and detection limits for a population frequency
#'
#' The uncertainty a frequency carries from the number of events behind it,
#' independent of where the thresholds were placed. Returned in percentage points
#' so it combines directly with the gate placement uncertainty.
#'
#' @param k events in the population; may be a vector
#' @param n events in the parent gate the population is expressed against;
#'   scalar or the same length as `k`
#' @param z coverage factor for the Wilson half-width. The default of 1 gives a
#'   standard uncertainty, matching the convention used for the gate terms
#' @param lod_events events below which a population is not called detected
#' @param loq_events events below which a population is detected but not
#'   quantified
#' @return a data.frame with one row per element of `k`, carrying the counts, the
#'   standard uncertainty, both limits expressed as percentages of the parent
#'   gate, and a verdict
#' @export
counting_uncertainty <- function(k, n, z = 1, lod_events = 20L,
                                 loq_events = 50L) {
  k <- as.numeric(k)
  n <- rep_len(as.numeric(n), length(k))
  ok <- is.finite(k) & is.finite(n) & n > 0 & k >= 0 & k <= n

  u <- rep(NA_real_, length(k))
  u[ok] <- 100 * z * sqrt(k[ok] * (n[ok] - k[ok]) / n[ok] + z^2 / 4) /
    (n[ok] + z^2)

  lod <- ifelse(ok, 100 * lod_events / n, NA_real_)
  loq <- ifelse(ok, 100 * loq_events / n, NA_real_)

  data.frame(
    n_cells = k,
    n_parent_events = n,
    u_counting_pct_points = round(u, 4),
    lod_pct = round(lod, 4),
    loq_pct = round(loq, 4),
    detection = ifelse(!ok, NA_character_,
                ifelse(k < lod_events, "below LOD",
                ifelse(k < loq_events, "detected, below LOQ", "quantified"))),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Propagate threshold uncertainty into one sample's population frequencies
#'
#' For every marker a population's definition reads, the frequency is re-scored
#' with that marker's threshold displaced by plus and minus its standard
#' uncertainty, and half the resulting spread is taken as that marker's
#' contribution. Contributions are summed in quadrature across markers, which is
#' the GUM treatment of independent terms and is the same arithmetic the operator
#' studies apply to manual gates.
#'
#' THE PARENT TERM IS INCLUDED AND MATTERS MOST. Displacing the CD45 cut moves
#' cells into and out of the denominator every population is expressed against,
#' so it perturbs every population at once. The operator work finds the first
#' gate dominates the budget; reporting the per-marker terms without it would
#' leave out the largest one.
#'
#' Markers whose uncertainty is NA contribute nothing and are counted in
#' `n_terms_missing`, so a small total that is small only because most terms
#' could not be computed is distinguishable from a genuinely tight one.
#'
#' THE COUNTING TERM IS REPORTED BESIDE THE GATE TERMS, NOT MIXED INTO THEM.
#' `u_pct_points` keeps its meaning: gate placement only. `u_counting_pct_points`
#' is what the frequency carries from the number of events behind it, and
#' `u_total_pct_points` is their quadrature sum. Keeping the first column fixed
#' means every number this table published before is still the same number.
#'
#' The two are not strictly independent, since displacing a cut also changes the
#' count. Quadrature treats them as though they were, which is the same
#' approximation the GUM makes for the marker terms and is stated here rather
#' than left for the reader to find.
#'
#' @param tmat transformed marker matrix for one sample
#' @param thr named threshold vector for that sample
#' @param parent logical parent-gate mask
#' @param spec population specification
#' @param u named vector of standard uncertainties, one per marker
#' @param cd45_x CD45 values for the sample, or NULL when the gate was skipped
#' @param cd45_live logical mask of live cells, the parent of the CD45 gate
#' @param cd45_threshold the CD45 cut
#' @param cd45_u standard uncertainty of the CD45 cut
#' @param lod_events events below which a population is not called detected
#' @param loq_events events below which a population is detected but not
#'   quantified
#' @return list(per_population = data.frame, budget = data.frame)
#' @export
population_frequency_uncertainty <- function(tmat, thr, parent, spec, u,
                                             cd45_x = NULL, cd45_live = NULL,
                                             cd45_threshold = NA_real_,
                                             cd45_u = NA_real_,
                                             lod_events = 20L,
                                             loq_events = 50L) {
  pct_at <- function(thr2, parent2) {
    hi <- derive_intermediate_bounds(tmat, thr2, parent2, spec)
    sp <- score_populations(tmat, thr2, parent2, spec, hi_thr = hi)
    den <- max(1L, sum(parent2))
    vapply(sp$masks, function(m) 100 * sum(m) / den, numeric(1))
  }

  # The base case is scored once and its event counts kept, rather than being
  # recovered from the percentage afterwards. Same call, same arithmetic as
  # pct_at(thr, parent), so `base` is unchanged; the counts come out of it free.
  n_parent <- sum(parent)
  hi0  <- derive_intermediate_bounds(tmat, thr, parent, spec)
  sp0  <- score_populations(tmat, thr, parent, spec, hi_thr = hi0)
  k0   <- vapply(sp0$masks, sum, integer(1))
  base <- 100 * k0 / max(1L, n_parent)
  if (!length(base)) return(NULL)
  pops <- names(base)

  # Only markers a population actually reads can move it. Perturbing the rest
  # would cost time and add nothing but rounding noise to the budget.
  markers_of <- lapply(spec, function(d)
    unique(c(setdiff(names(d), "any_of"), names(d[["any_of"]]))))

  terms <- list()
  for (m in intersect(names(thr), names(u))) {
    um <- u[[m]]
    if (!is.finite(um) || um <= 0) next
    lo <- thr; lo[[m]] <- thr[[m]] - um
    hi <- thr; hi[[m]] <- thr[[m]] + um
    d <- abs(pct_at(hi, parent) - pct_at(lo, parent)) / 2
    terms[[m]] <- d[pops]
  }

  if (!is.null(cd45_x) && is.finite(cd45_threshold) && is.finite(cd45_u) &&
      cd45_u > 0 && !is.null(cd45_live)) {
    plo <- cd45_live & cd45_x > (cd45_threshold - cd45_u)
    phi <- cd45_live & cd45_x > (cd45_threshold + cd45_u)
    if (sum(plo) > 0L && sum(phi) > 0L)
      terms[["parent_CD45"]] <-
        abs(pct_at(thr, phi) - pct_at(thr, plo))[pops] / 2
  }

  budget <- list(); u_pop <- setNames(rep(NA_real_, length(pops)), pops)
  n_used <- setNames(integer(length(pops)), pops)
  n_miss <- setNames(integer(length(pops)), pops)
  for (p in pops) {
    rel <- c(markers_of[[p]] %||% character(0), "parent_CD45")
    v <- vapply(intersect(names(terms), rel), function(k) terms[[k]][[p]],
                numeric(1))
    v <- v[is.finite(v)]
    u_pop[[p]] <- if (length(v)) sqrt(sum(v^2)) else NA_real_
    n_used[[p]] <- length(v)
    n_miss[[p]] <- length(setdiff(rel, c(names(terms), "parent_CD45")))
    for (k in names(v))
      budget[[length(budget) + 1L]] <- data.frame(
        population = p, term = k, u_pct_points = round(unname(v[[k]]), 4),
        stringsAsFactors = FALSE)
  }

  # Quadrature over the finite components, matching how threshold_uncertainty()
  # forms u_combined. Where neither term could be computed the total is NA rather
  # than zero: an uncertainty nobody could estimate is not a small one.
  cu <- counting_uncertainty(k0, n_parent, lod_events = lod_events,
                             loq_events = loq_events)
  u_tot <- vapply(seq_along(pops), function(i) {
    v <- c(u_pop[[i]], cu$u_counting_pct_points[[i]])
    v <- v[is.finite(v)]
    if (length(v)) sqrt(sum(v^2)) else NA_real_
  }, numeric(1))

  list(per_population = data.frame(
         population = pops, pct_of_cd45_pos = unname(base),
         u_pct_points = round(unname(u_pop), 4),
         n_terms = unname(n_used), n_terms_missing = unname(n_miss),
         n_cells = unname(k0), n_parent_events = as.integer(n_parent),
         u_counting_pct_points = cu$u_counting_pct_points,
         u_total_pct_points = round(u_tot, 4),
         lod_pct = cu$lod_pct, loq_pct = cu$loq_pct,
         detection = cu$detection,
         row.names = NULL, stringsAsFactors = FALSE),
       budget = if (length(budget)) do.call(rbind, budget) else NULL)
}

#' Run the whole uncertainty analysis over a scored cohort
#'
#' Driver called by [run_cyraven()] after populations are scored. It touches
#' nothing the pipeline has already computed.
#'
#' @param pops the per-sample scoring results
#' @param gates the per-sample gate objects
#' @param verdicts the per-sample staining verdicts
#' @param panel_of named character vector mapping sample to panel
#' @param spec population specification
#' @param B bootstrap replicates
#' @param seed base seed
#' @param max_events cap on events per replicate
#' @param lod_events events below which a population is not called detected
#' @param loq_events events below which a population is detected but not
#'   quantified
#' @return list(thresholds, frequencies, budget), each a data.frame or NULL
#' @export
run_gate_uncertainty <- function(pops, gates, verdicts, panel_of, spec,
                                 B = 100L, seed = 42L, max_events = 20000L,
                                 lod_events = 20L, loq_events = 50L) {
  # Each (sample, marker) draws from its own stream, keyed on the two names, so a
  # threshold's uncertainty does not depend on how many markers happened to be
  # processed before it. Reordering a panel must not change a published number.
  seed_for <- function(key) {
    h <- 0L
    for (b in as.integer(charToRaw(key))) h <- (h * 31L + b) %% 100000L
    as.integer(seed) + h
  }
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  trows <- list(); frows <- list(); brows <- list()
  grid <- valley_setting_grid()

  for (s in names(pops)) {
    P <- pops[[s]]; g <- gates[[s]]
    if (is.null(P$tmat) || is.null(g)) next
    parent <- g$masks$cd45_pos
    thr <- P$thresholds
    u <- setNames(rep(NA_real_, length(thr)), names(thr))

    for (m in names(thr)) {
      src <- P$details[[m]]$source %||% "valley"
      tu <- threshold_uncertainty(P$tmat[parent, m], source = src, B = B,
                                  seed = seed_for(paste0(s, "\r", m)),
                                  max_events = max_events, grid = grid)
      u[[m]] <- tu$u_combined
      trows[[length(trows) + 1L]] <- data.frame(
        sample_id = s, panel = panel_of[[s]] %||% NA_character_, marker = m,
        threshold = unname(thr[[m]]), source = src,
        u_sampling = round(tu$u_sampling, 5), u_method = round(tu$u_method, 5),
        u_combined = round(tu$u_combined, 5),
        valley_rel_depth = round(tu$rel_depth, 4),
        bootstrap_valley_rate = round(tu$valley_rate, 3),
        n_parent_events = tu$n_events, basis = tu$basis,
        stringsAsFactors = FALSE)
    }

    cd45_u <- NA_real_
    if (!is.null(g$cd45_x) && is.finite(g$cd45_threshold %||% NA_real_)) {
      cu <- threshold_uncertainty(g$cd45_x[g$masks$live_cells],
                                  source = g$cd45_source %||% "valley", B = B,
                                  seed = seed_for(paste0(s, "\rCD45")),
                                  max_events = max_events, grid = grid)
      cd45_u <- cu$u_combined
      trows[[length(trows) + 1L]] <- data.frame(
        sample_id = s, panel = panel_of[[s]] %||% NA_character_,
        marker = "CD45 (parent gate)", threshold = g$cd45_threshold,
        source = g$cd45_source %||% NA_character_,
        u_sampling = round(cu$u_sampling, 5), u_method = round(cu$u_method, 5),
        u_combined = round(cu$u_combined, 5),
        valley_rel_depth = round(cu$rel_depth, 4),
        bootstrap_valley_rate = round(cu$valley_rate, 3),
        n_parent_events = cu$n_events, basis = cu$basis,
        stringsAsFactors = FALSE)
    }

    fu <- tryCatch(population_frequency_uncertainty(
      P$tmat, thr, parent, spec, u, cd45_x = g$cd45_x,
      cd45_live = g$masks$live_cells, cd45_threshold = g$cd45_threshold,
      cd45_u = cd45_u, lod_events = lod_events, loq_events = loq_events),
      error = function(e) NULL)
    if (is.null(fu)) next

    pp <- fu$per_population
    pp$sample_id <- s
    pp$panel <- panel_of[[s]] %||% NA_character_
    pp$qc_status <- verdicts[[s]]$qc_status %||% "pass"
    frows[[length(frows) + 1L]] <- pp
    if (!is.null(fu$budget)) {
      bb <- fu$budget; bb$sample_id <- s
      brows[[length(brows) + 1L]] <- bb
    }
  }

  list(thresholds  = if (length(trows)) do.call(rbind, trows) else NULL,
       frequencies = if (length(frows)) do.call(rbind, frows) else NULL,
       budget      = if (length(brows)) do.call(rbind, brows) else NULL)
}

#' Compare a between-group difference against the gate uncertainty behind it
#'
#' Adds two columns to a group comparison table: the typical within-sample
#' uncertainty on the quantity being compared, and the ratio of the observed
#' difference to it.
#'
#' READ IT AS A SCREEN, NOT A TEST. The gate uncertainty is partly common to
#' every sample in a run, since they share a panel, a transform and a placement
#' rule, so it cancels to some degree in a difference and the ratio is
#' conservative. A ratio below 1 says the groups differ by less than the typical
#' distance the cut itself can move, which is a reason to look at
#' threshold_uncertainty.csv before interpreting the result, not a p-value.
#'
#' TWO RATIOS, NOT ONE. `difference_over_gate_u` compares the difference against
#' gate placement alone and keeps exactly the value it has always had.
#' `difference_over_total_u` compares it against placement and counting together,
#' and is the stricter of the two. They separate for a rare population, where the
#' cut can be well placed and the frequency still be built on too few events: the
#' first ratio passes and the second does not, and the second is the one to
#' believe.
#'
#' @param gstats a table from [stats_group_comparison()]
#' @param ufreq the `frequencies` element of [run_gate_uncertainty()]
#' @return `gstats` with `gate_u_pct_points`, `difference_over_gate_u`,
#'   `total_u_pct_points` and `difference_over_total_u` appended, unchanged if
#'   the uncertainty table is missing
#' @export
annotate_gate_uncertainty <- function(gstats, ufreq) {
  if (is.null(gstats) || !nrow(gstats)) return(gstats)
  gstats$gate_u_pct_points   <- NA_real_
  gstats$difference_over_gate_u <- NA_real_
  gstats$total_u_pct_points  <- NA_real_
  gstats$difference_over_total_u <- NA_real_
  if (is.null(ufreq) || !nrow(ufreq)) return(gstats)
  pass <- (ufreq$qc_status %||% "pass") == "pass"
  d <- abs(gstats$median_comparison - gstats$median_reference)

  keep <- ufreq[pass & is.finite(ufreq$u_pct_points), , drop = FALSE]
  if (nrow(keep)) {
    um <- tapply(keep$u_pct_points, keep$population, median, na.rm = TRUE)
    u <- unname(um[gstats$population])
    gstats$gate_u_pct_points <- round(u, 4)
    gstats$difference_over_gate_u <-
      ifelse(is.finite(u) & u > 0, round(d / u, 2), NA_real_)
  }

  # Absent when the run predates the counting term, so guard on the column
  # rather than assuming the shape of the table handed in.
  if ("u_total_pct_points" %in% names(ufreq)) {
    k2 <- ufreq[pass & is.finite(ufreq$u_total_pct_points), , drop = FALSE]
    if (nrow(k2)) {
      tm <- tapply(k2$u_total_pct_points, k2$population, median, na.rm = TRUE)
      ut <- unname(tm[gstats$population])
      gstats$total_u_pct_points <- round(ut, 4)
      gstats$difference_over_total_u <-
        ifelse(is.finite(ut) & ut > 0, round(d / ut, 2), NA_real_)
    }
  }
  gstats
}
