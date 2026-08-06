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
  if (grepl("^control_q995", source %||% "")) {
    out$basis <- "control tube not available here"
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
#' @param tmat transformed marker matrix for one sample
#' @param thr named threshold vector for that sample
#' @param parent logical parent-gate mask
#' @param spec population specification
#' @param u named vector of standard uncertainties, one per marker
#' @param cd45_x CD45 values for the sample, or NULL when the gate was skipped
#' @param cd45_live logical mask of live cells, the parent of the CD45 gate
#' @param cd45_threshold the CD45 cut
#' @param cd45_u standard uncertainty of the CD45 cut
#' @return list(per_population = data.frame, budget = data.frame)
#' @export
population_frequency_uncertainty <- function(tmat, thr, parent, spec, u,
                                             cd45_x = NULL, cd45_live = NULL,
                                             cd45_threshold = NA_real_,
                                             cd45_u = NA_real_) {
  pct_at <- function(thr2, parent2) {
    hi <- derive_intermediate_bounds(tmat, thr2, parent2, spec)
    sp <- score_populations(tmat, thr2, parent2, spec, hi_thr = hi)
    den <- max(1L, sum(parent2))
    vapply(sp$masks, function(m) 100 * sum(m) / den, numeric(1))
  }

  base <- pct_at(thr, parent)
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

  list(per_population = data.frame(
         population = pops, pct_of_cd45_pos = unname(base),
         u_pct_points = round(unname(u_pop), 4),
         n_terms = unname(n_used), n_terms_missing = unname(n_miss),
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
#' @return list(thresholds, frequencies, budget), each a data.frame or NULL
#' @export
run_gate_uncertainty <- function(pops, gates, verdicts, panel_of, spec,
                                 B = 100L, seed = 42L, max_events = 20000L) {
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
      cd45_u = cd45_u), error = function(e) NULL)
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
#' @param gstats a table from [stats_group_comparison()]
#' @param ufreq the `frequencies` element of [run_gate_uncertainty()]
#' @return `gstats` with `gate_u_pct_points` and `difference_over_gate_u`
#'   appended, unchanged if the uncertainty table is missing
#' @export
annotate_gate_uncertainty <- function(gstats, ufreq) {
  if (is.null(gstats) || !nrow(gstats)) return(gstats)
  gstats$gate_u_pct_points   <- NA_real_
  gstats$difference_over_gate_u <- NA_real_
  if (is.null(ufreq) || !nrow(ufreq)) return(gstats)
  keep <- ufreq[(ufreq$qc_status %||% "pass") == "pass" &
                  is.finite(ufreq$u_pct_points), , drop = FALSE]
  if (!nrow(keep)) return(gstats)
  um <- tapply(keep$u_pct_points, keep$population, median, na.rm = TRUE)
  u <- unname(um[gstats$population])
  d <- abs(gstats$median_comparison - gstats$median_reference)
  gstats$gate_u_pct_points <- round(u, 4)
  gstats$difference_over_gate_u <-
    ifelse(is.finite(u) & u > 0, round(d / u, 2), NA_real_)
  gstats
}
