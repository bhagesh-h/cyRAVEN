# SECTION 7d -- TIMEPOINT FIGURES
# =============================================================================
#
# WHY THESE ARE NOT JUST "THE GROUP FIGURES WITH A DIFFERENT COLUMN". A
# timepoint column is not a grouping in the sense the between-group figures
# mean. d0, d3 and d7 come from the SAME patients, so the samples are paired,
# not independent, and a boxplot per timepoint throws away the pairing -- the
# only structure a longitudinal design has. Two cohorts can have identical
# boxplots while every individual rises in one and falls in the other.
#
# The visualisation the longitudinal immune-monitoring literature settles on is
# the per-subject trajectory: one line per patient across visits, with the group
# summary drawn over it rather than instead of it. Time-series flow cytometry is
# still recognised as methodologically unsettled at the modelling level
# (Fukushima et al., CYBERTRACK, on population tracking across time), but on the
# display question there is broad agreement: keep the subject visible, and read
# each timepoint against that subject's own baseline rather than against the
# cohort mean. Longitudinal sepsis work follows the same convention -- reporting
# d3 and d7 against each patient's d0 -- because the between-patient spread is
# larger than the within-patient change, so pooling hides the effect.
#
# WHAT IS DELIBERATELY NOT HERE. No p-values. On this cohort the timepoints are
# 9/6/5 samples from 9 patients, and the tests that would be honest about the
# pairing (mixed model, Friedman) need more complete triplets than exist. These
# figures are descriptive, and `--no-group-tests` is set on the runs that write
# them.

#' Per-patient trajectories across timepoints, one panel per population
#'
#' The paired view: each line is one patient followed through their visits, so
#' a consistent within-patient direction is visible even when the cohort spread
#' swamps it. The cohort median is drawn over the lines, not in place of them.
#'
#' @param freq population_frequencies table.
#' @param smap Sample map carrying `timepoint` and `patient_id`.
#' @param outfile Destination PNG.
#' @param value_col Column to plot.
#' @param log_scale Use a log10 y axis. Population sizes span orders of
#'   magnitude, so a fold change is the comparable quantity.
#' @param max_pops Facet cap, largest populations first.
#' @keywords internal
fig_timepoint_trajectories <- function(freq, smap, outfile,
                                       value_col = "pct_of_cd45_pos",
                                       log_scale = TRUE, max_pops = 12L) {
  if (is.null(freq) || !nrow(freq) || !value_col %in% names(freq))
    return(invisible(NULL))
  if (is.null(smap) || !all(c("timepoint", "patient_id") %in% names(smap)))
    return(invisible(NULL))

  key <- smap[, c("sample_id", "timepoint", "patient_id")]
  d <- merge(freq, key, by = "sample_id")
  d <- d[!is.na(d$timepoint) & nzchar(as.character(d$timepoint)), , drop = FALSE]
  d$value <- suppressWarnings(as.numeric(d[[value_col]]))
  d <- d[is.finite(d$value), , drop = FALSE]
  if (log_scale) d <- d[d$value > 0, , drop = FALSE]
  if (!nrow(d) || length(unique(d$timepoint)) < 2L) return(invisible(NULL))

  ord <- names(sort(tapply(d$value, d$population, stats::median, na.rm = TRUE),
                    decreasing = TRUE))
  keep <- utils::head(ord, max_pops)
  d <- d[d$population %in% keep, , drop = FALSE]
  d$population <- factor(d$population, levels = keep)
  # Timepoints in collection order, not alphabetical. d10 would otherwise sit
  # between d0 and d3.
  d$timepoint <- natural_cluster_factor(d$timepoint)

  med <- stats::aggregate(value ~ population + timepoint, d, stats::median)

  nc <- min(4L, length(keep)); nr <- ceiling(length(keep) / nc)
  p <- ggplot2::ggplot(d, ggplot2::aes(timepoint, value)) +
    ggplot2::geom_line(ggplot2::aes(group = patient_id, colour = patient_id),
                       alpha = 0.65, linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = patient_id), size = 1.4, alpha = 0.85) +
    # The cohort median is context, not the subject of the figure. Drawn in
    # black at 1.1 it read as the primary series and the per-patient lines --
    # the whole point of a paired view -- became its backdrop. Grey and thin
    # enough to follow, light enough to stay behind them.
    ggplot2::geom_line(data = med, ggplot2::aes(group = 1), colour = "grey35",
                       linewidth = 0.7) +
    ggplot2::geom_point(data = med, colour = "grey35", size = 1.6) +
    ggplot2::facet_wrap(~ population, ncol = nc, scales = "free_y") +
    ggplot2::labs(
      title = paste0("Per-patient trajectories across timepoints (", value_col, ")"),
      subtitle = paste("one line per patient; grey is the cohort median.",
                       "Samples are PAIRED -- the same patients at each",
                       "timepoint -- so read the direction of each line, not",
                       "the spread between them"),
      caption = paste("Descriptive only: no test is drawn. The timepoints are",
                      "the same patients, and there are too few complete",
                      "triplets for a paired test to be meaningful."),
      x = NULL, y = value_col, colour = "patient") +
    theme_cyto() + theme_panel_borders()
  if (log_scale)
    p <- p + ggplot2::scale_y_log10(labels = function(x)
      format(x, big.mark = ",", scientific = FALSE, trim = TRUE))

  safe_ggsave(outfile, plot = p, width = max(9, 2.9 * nc + 1.8),
              height = max(4.5, 2.9 * nr + 1.6), dpi = 200, limitsize = FALSE)
  invisible(outfile)
}

#' One population's two subsets against each other, per timepoint
#'
#' WHY A DEDICATED FIGURE. A pair of subsets that divide one compartment --
#' Vd1 and Vd2 inside the gamma-delta compartment, CD4 and CD8 inside T cells --
#' is read as a balance, not as two independent numbers, and the published
#' comparisons report it that way: the subsets shift against each other while
#' the compartment stays flat. Two separate panels make that shift hard to see,
#' because the reader has to hold one panel in mind while looking at the other.
#'
#' @param freq population_frequencies table.
#' @param smap Sample map carrying `timepoint` and `patient_id`.
#' @param outfile Destination PNG.
#' @param pairs List of two-element character vectors naming the subsets.
#' @param value_col Column to plot.
#' @keywords internal
fig_subset_balance <- function(freq, smap, outfile,
                               pairs = list(c("Vd1 T cells", "Vd2 T cells"),
                                            c("CD4 T cells", "CD8 T cells")),
                               value_col = "pct_of_cd45_pos") {
  if (is.null(freq) || !nrow(freq) || !value_col %in% names(freq))
    return(invisible(NULL))
  if (is.null(smap) || !all(c("timepoint", "patient_id") %in% names(smap)))
    return(invisible(NULL))
  pairs <- Filter(function(pr) all(pr %in% freq$population), pairs)
  if (!length(pairs)) return(invisible(NULL))

  key <- smap[, c("sample_id", "timepoint", "patient_id")]
  rows <- lapply(pairs, function(pr) {
    f <- freq[freq$population %in% pr, c("sample_id", "population", value_col)]
    names(f)[3] <- "value"
    w <- stats::reshape(f, idvar = "sample_id", timevar = "population",
                        direction = "wide")
    names(w) <- sub("^value[.]", "", names(w))
    if (!all(pr %in% names(w))) return(NULL)
    m <- merge(w, key, by = "sample_id")
    a <- suppressWarnings(as.numeric(m[[pr[1]]]))
    b <- suppressWarnings(as.numeric(m[[pr[2]]]))
    ok <- is.finite(a) & is.finite(b) & a > 0 & b > 0
    if (!any(ok)) return(NULL)
    data.frame(pair = paste(pr[1], "vs", pr[2]),
               timepoint = m$timepoint[ok], patient_id = m$patient_id[ok],
               ratio = a[ok] / b[ok], stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(d) || !nrow(d)) return(invisible(NULL))
  d$timepoint <- natural_cluster_factor(d$timepoint)

  p <- ggplot2::ggplot(d, ggplot2::aes(timepoint, ratio)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(group = patient_id, colour = patient_id),
                       alpha = 0.7, linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = patient_id), size = 1.8) +
    ggplot2::facet_wrap(~ pair, scales = "free_y") +
    # A ratio is symmetric in logs and badly asymmetric in linear space: 2 and
    # 0.5 are the same magnitude of shift in opposite directions, and only a log
    # axis places them the same distance from the dashed line of equality.
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Subset balance across timepoints",
                  subtitle = paste("ratio of the two subsets within one",
                                   "compartment; the dashed line is parity.",
                                   "Log axis, so a doubling and a halving are",
                                   "the same distance from it"),
                  x = NULL, y = "ratio (log scale)", colour = "patient") +
    theme_cyto() + theme_panel_borders()
  safe_ggsave(outfile, plot = p, width = max(8, 4.2 * length(pairs) + 1.5),
              height = 5, dpi = 200)
  invisible(outfile)
}

#' Every marker's median intensity per population, split by timepoint
#'
#' The state heatmap asks what a population expresses; this asks whether that
#' changed over the admission. Split rather than pooled, because pooling the
#' timepoints is what hides a change that happens in all patients at once.
#' @param mfi population_marker_mfi table (median_asinh per sample x population x marker).
#' @param smap Sample map carrying `timepoint`.
#' @param outfile Destination PNG.
#' @param max_pops Population cap, largest first.
#' @keywords internal
fig_marker_by_timepoint <- function(mfi, smap, outfile, max_pops = 8L) {
  if (is.null(mfi) || !nrow(mfi) || !"median_asinh" %in% names(mfi))
    return(invisible(NULL))
  if (is.null(smap) || !"timepoint" %in% names(smap)) return(invisible(NULL))
  d <- merge(mfi, smap[, c("sample_id", "timepoint")], by = "sample_id")
  d <- d[!is.na(d$timepoint) & nzchar(as.character(d$timepoint)) &
           is.finite(d$median_asinh), , drop = FALSE]
  if (!nrow(d) || length(unique(d$timepoint)) < 2L) return(invisible(NULL))

  big <- names(sort(table(d$population), decreasing = TRUE))[seq_len(
    min(max_pops, length(unique(d$population))))]
  d <- d[d$population %in% big, , drop = FALSE]
  d$timepoint <- natural_cluster_factor(d$timepoint)
  agg <- stats::aggregate(median_asinh ~ population + marker + timepoint, d,
                          stats::median)

  p <- ggplot2::ggplot(agg, ggplot2::aes(timepoint, marker, fill = median_asinh)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::facet_wrap(~ population, nrow = 1) +
    ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::labs(title = "Marker intensity per population, by timepoint",
                  subtitle = paste("cohort median of each sample's median",
                                   "intensity; split by timepoint rather than",
                                   "pooled, so a shift shared by every patient",
                                   "is visible"),
                  x = NULL, y = NULL, fill = "median\n(asinh)") +
    theme_cyto() + theme_panel_borders()
  safe_ggsave(outfile, plot = p,
              width = max(9, 1.7 * length(big) + 2.5),
              height = max(5, 0.22 * length(unique(agg$marker)) + 2.2),
              dpi = 200, limitsize = FALSE)
  invisible(outfile)
}
