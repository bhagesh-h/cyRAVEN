# SECTION 9 -- SPECIFICATION CONFORMANCE ACROSS RUNS
# =============================================================================
#
# WHAT THIS ADDS THAT threshold_scale_qc.csv DOES NOT. That check compares each
# threshold against the other samples in the same run. It catches one bad tube
# among twenty good ones. It cannot catch the case where the whole run has moved:
# if every sample drifts together, the leave-one-out peer median drifts with it
# and every sample looks fine relative to its peers. A cohort acquired six months
# later, after a laser service or a reagent lot change, is exactly that case.
#
# So this file adds a second reference, and keeps the first one. A baseline
# written from an accepted run records where each marker's threshold sat, how
# variable it was, how often it needed the quantile fallback, and what the
# populations came out at. Later runs are then measured against that fixed point
# as well as against themselves.
#
# WHY IT IS A SEPARATE FILE AND NOT A NEW REFERENCE INSIDE THE OLD ONE.
# threshold_scale_qc.csv keeps its within-run meaning unchanged. Repointing that
# file at a baseline would have kept the filename and the column names while
# silently changing what the numbers mean, which is worse than a rename: anyone
# comparing an old copy against a new one would be comparing two different
# quantities without being told.
#
# WHAT A FAILURE MEANS. Not that the run is bad. It means the gating
# specification no longer places its cuts where it did when the baseline was
# accepted, so the frequencies from the two runs are not the same measurement and
# should not be pooled or compared until someone has looked. The verdict is
# advisory by default; --fail-on-drift makes it an exit code for scheduled runs.

#' Canonical text form of a population specification
#'
#' Used to detect that the specification itself changed between runs, which makes
#' a threshold comparison meaningless: a population defined differently is a
#' different population, however similar its frequency.
#'
#' The full text is stored rather than a hash, so a mismatch can name the
#' population that changed instead of only asserting that something did.
#'
#' @param spec population specification
#' @return named character vector, one canonical string per population
#' @export
spec_fingerprint <- function(spec) {
  if (!length(spec)) return(character(0))
  out <- vapply(spec, function(d) {
    if (!length(d)) return("")
    keys <- sort(setdiff(names(d), "any_of"))
    main <- paste(keys, vapply(keys, function(k) paste(as.character(d[[k]]),
                                                       collapse = "|"),
                               character(1)), sep = "=", collapse = ",")
    ao <- d[["any_of"]]
    if (!is.null(ao) && length(ao)) {
      ak <- sort(names(ao))
      main <- paste0(main, ";any_of(",
                     paste(ak, vapply(ak, function(k)
                       paste(as.character(ao[[k]]), collapse = "|"),
                       character(1)), sep = "=", collapse = ","), ")")
    }
    main
  }, character(1))
  out[order(names(out))]
}

#' Write a conformance baseline from an accepted run
#'
#' Summaries, not cells: the file holds per-marker threshold location and spread,
#' the fallback rate, per-population frequency location and spread, and the
#' specification text. It carries no event-level data and no patient data, so it
#' can be version-controlled next to the config it describes.
#'
#' @param path destination `.rds`
#' @param thr_all the thresholds table (`thresholds_used.csv` in memory)
#' @param freq the frequency table
#' @param spec population specification
#' @param opt the resolved option list
#' @param cofactors per-panel cofactors
#' @return `path`, invisibly
#' @export
write_spec_baseline <- function(path, thr_all, freq, spec, opt = list(),
                                cofactors = list()) {
  mk <- NULL
  if (!is.null(thr_all) && nrow(thr_all)) {
    key <- paste(thr_all$panel, thr_all$marker, sep = "\r")
    mk <- do.call(rbind, lapply(split(thr_all, key), function(d) data.frame(
      panel = d$panel[1], marker = d$marker[1], n = nrow(d),
      median_threshold = median(d$threshold, na.rm = TRUE),
      mad_threshold = mad(d$threshold, na.rm = TRUE),
      fallback_rate = mean(d$source == "quantile_fallback", na.rm = TRUE),
      stringsAsFactors = FALSE)))
    rownames(mk) <- NULL
  }
  pp <- NULL
  if (!is.null(freq) && nrow(freq)) {
    f <- freq[(freq$qc_status %||% "pass") == "pass" & !freq$is_control, ,
              drop = FALSE]
    if (nrow(f)) {
      pp <- do.call(rbind, lapply(split(f, f$population), function(d) data.frame(
        population = d$population[1], n = nrow(d),
        median_pct = median(d$pct_of_cd45_pos, na.rm = TRUE),
        mad_pct = mad(d$pct_of_cd45_pos, na.rm = TRUE),
        stringsAsFactors = FALSE)))
      rownames(pp) <- NULL
    }
  }
  b <- list(
    format = 1L,
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    cyraven_version = as.character(utils::packageVersion("cyRAVEN")),
    transform = opt$transform %||% "arcsinh",
    cofactors = cofactors,
    spec = spec_fingerprint(spec),
    markers = mk,
    populations = pp)
  saveRDS(b, path)
  invisible(path)
}

#' Read a conformance baseline
#'
#' @param path path written by [write_spec_baseline()]
#' @return the baseline list
#' @export
read_spec_baseline <- function(path) {
  if (!file.exists(path))
    stop("baseline not found: ", path, call. = FALSE)
  b <- readRDS(path)
  if (!is.list(b) || is.null(b$format))
    stop("not a cyRAVEN baseline: ", path, call. = FALSE)
  b
}

#' Test this run against a baseline
#'
#' Compares each marker's threshold location against where the baseline put it,
#' scaled by the baseline's own spread for that marker, so the tolerance adapts to
#' how variable the marker legitimately is. The same comparison is made for
#' population frequencies.
#'
#' THE SCALE IS FLOORED, for the reason it is floored in the within-run check: a
#' marker whose baseline samples happened to agree closely would otherwise divide
#' by near zero and fail on a difference of no consequence.
#'
#' A CHANGED TRANSFORM INVALIDATES EVERY THRESHOLD COMPARISON, because thresholds
#' are expressed on the analysis scale. That case is reported as its own verdict
#' rather than as dozens of drifted markers.
#'
#' @param thr_all this run's thresholds table
#' @param freq this run's frequency table
#' @param spec this run's population specification
#' @param baseline the list from [read_spec_baseline()]
#' @param transform this run's transform name
#' @param z_qualify robust z above which a marker is qualified
#' @param z_fail robust z above which a marker fails
#' @return list(markers, populations, spec_changes, summary)
#' @export
specification_conformance <- function(thr_all, freq, spec, baseline,
                                      transform = "arcsinh",
                                      z_qualify = 3.5, z_fail = 6) {
  verdict_of <- function(z) ifelse(!is.finite(z), "not comparable",
                            ifelse(z > z_fail, "fail",
                            ifelse(z > z_qualify, "qualify", "pass")))

  scale_changed <- !identical(as.character(transform),
                              as.character(baseline$transform %||% "arcsinh"))

  # ---- specification drift, checked before anything derived from it ----------
  now <- spec_fingerprint(spec); was <- baseline$spec %||% character(0)
  added   <- setdiff(names(now), names(was))
  removed <- setdiff(names(was), names(now))
  common  <- intersect(names(now), names(was))
  altered <- common[now[common] != was[common]]
  sc <- NULL
  if (length(c(added, removed, altered)))
    sc <- data.frame(
      population = c(added, removed, altered),
      change = c(rep("added", length(added)), rep("removed", length(removed)),
                 rep("redefined", length(altered))),
      baseline_definition = c(rep(NA_character_, length(added)),
                              unname(was[removed]), unname(was[altered])),
      current_definition = c(unname(now[added]), rep(NA_character_, length(removed)),
                             unname(now[altered])),
      row.names = NULL, stringsAsFactors = FALSE)

  # ---- per-marker threshold location ----------------------------------------
  mk <- NULL
  bm <- baseline$markers
  if (!is.null(thr_all) && nrow(thr_all) && !is.null(bm) && nrow(bm)) {
    key <- paste(thr_all$panel, thr_all$marker, sep = "\r")
    cur <- do.call(rbind, lapply(split(thr_all, key), function(d) data.frame(
      panel = d$panel[1], marker = d$marker[1], n = nrow(d),
      median_threshold = median(d$threshold, na.rm = TRUE),
      fallback_rate = mean(d$source == "quantile_fallback", na.rm = TRUE),
      stringsAsFactors = FALSE)))
    rownames(cur) <- NULL
    m <- merge(cur, bm, by = c("panel", "marker"), suffixes = c("", "_baseline"),
               all.x = TRUE)
    s <- pmax(m$mad_threshold, 0.25)
    z <- abs(m$median_threshold - m$median_threshold_baseline) / s
    mk <- data.frame(
      panel = m$panel, marker = m$marker,
      n_samples = m$n, n_samples_baseline = m$n_baseline,
      median_threshold = round(m$median_threshold, 4),
      baseline_threshold = round(m$median_threshold_baseline, 4),
      baseline_mad = round(m$mad_threshold, 4),
      robust_z = round(z, 2),
      fallback_rate = round(m$fallback_rate, 3),
      baseline_fallback_rate = round(m$fallback_rate_baseline, 3),
      verdict = if (scale_changed) "not comparable" else verdict_of(z),
      note = ifelse(is.na(m$median_threshold_baseline), "absent from baseline",
             ifelse(scale_changed, "transform differs from baseline",
             ifelse(m$fallback_rate > (m$fallback_rate_baseline + 0.25),
                    "more samples fell back to a quantile than in the baseline",
                    NA_character_))),
      row.names = NULL, stringsAsFactors = FALSE)
    mk <- mk[order(-ifelse(is.finite(mk$robust_z), mk$robust_z, -1)), ]
  }

  # ---- per-population frequency location ------------------------------------
  pp <- NULL
  bp <- baseline$populations
  if (!is.null(freq) && nrow(freq) && !is.null(bp) && nrow(bp)) {
    f <- freq[(freq$qc_status %||% "pass") == "pass" & !freq$is_control, ,
              drop = FALSE]
    if (nrow(f)) {
      cur <- do.call(rbind, lapply(split(f, f$population), function(d) data.frame(
        population = d$population[1], n = nrow(d),
        median_pct = median(d$pct_of_cd45_pos, na.rm = TRUE),
        stringsAsFactors = FALSE)))
      rownames(cur) <- NULL
      m <- merge(cur, bp, by = "population", suffixes = c("", "_baseline"),
                 all.x = TRUE)
      # Floored at 0.5 percentage points: a population whose baseline samples
      # agreed to within a tenth of a point is not thereby entitled to fail the
      # run over a quarter-point move.
      s <- pmax(m$mad_pct, 0.5)
      z <- abs(m$median_pct - m$median_pct_baseline) / s
      # A population whose definition changed is not drifting, it is a different
      # population, and the z would be a category error.
      redefined <- !is.null(sc) & m$population %in%
        (if (is.null(sc)) character(0) else sc$population[sc$change == "redefined"])
      pp <- data.frame(
        population = m$population, n_samples = m$n,
        n_samples_baseline = m$n_baseline,
        median_pct = round(m$median_pct, 3),
        baseline_pct = round(m$median_pct_baseline, 3),
        baseline_mad = round(m$mad_pct, 3),
        robust_z = round(ifelse(redefined, NA_real_, z), 2),
        verdict = ifelse(redefined, "not comparable", verdict_of(z)),
        note = ifelse(redefined, "definition changed since the baseline",
               ifelse(is.na(m$median_pct_baseline), "absent from baseline",
                      NA_character_)),
        row.names = NULL, stringsAsFactors = FALSE)
      pp <- pp[order(-ifelse(is.finite(pp$robust_z), pp$robust_z, -1)), ]
    }
  }

  n_fail <- sum(c(mk$verdict, pp$verdict) == "fail", na.rm = TRUE)
  n_qual <- sum(c(mk$verdict, pp$verdict) == "qualify", na.rm = TRUE)
  list(markers = mk, populations = pp, spec_changes = sc,
       summary = list(n_fail = n_fail, n_qualify = n_qual,
                      transform_changed = scale_changed,
                      spec_changed = !is.null(sc),
                      baseline_created = baseline$created %||% NA_character_))
}
