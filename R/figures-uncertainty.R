# SECTION 9d -- GATE UNCERTAINTY FIGURES
# =============================================================================

#' Population frequencies with the uncertainty of the gate that produced them
#'
#' One row per population, one point per sample, and a bar through each point
#' spanning plus and minus the standard uncertainty of that sample's frequency
#' propagated from where its thresholds were placed.
#'
#' WHAT TO LOOK FOR. The question this figure answers is whether the spread
#' BETWEEN samples is larger than the bar WITHIN each of them. Where it is, the
#' variation is something the gate is measuring. Where the bars are as wide as
#' the scatter, the variation is the gate moving, and a group difference in that
#' population is not interpretable however small its p-value.
#'
#' Controls and QC-failed samples are excluded here rather than labelled, unlike
#' the QC figures: their percentages are fractions of an arbitrary parent slice,
#' so plotting them next to real samples would put two different quantities on
#' one axis.
#'
#' @param ufreq the `frequencies` element of [run_gate_uncertainty()]
#' @param outfile Path to write the figure to.
#' @param group_of Named character vector mapping sample_id to group label.
#' @param dpi Resolution in dots per inch. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. Default `fcs_colors()`.
#' @export
fig_frequency_uncertainty <- function(ufreq, outfile, group_of = NULL, dpi = 200,
                                      colors = fcs_colors()) {
  if (is.null(ufreq) || !nrow(ufreq)) {
    log_msg("[fig] no uncertainty table, frequency uncertainty figure skipped")
    return(invisible(NULL))
  }
  d <- ufreq[(ufreq$qc_status %||% "pass") == "pass" &
               is.finite(ufreq$pct_of_cd45_pos), , drop = FALSE]
  if (!nrow(d)) {
    log_msg("[fig] no QC-passing samples, frequency uncertainty figure skipped")
    return(invisible(NULL))
  }
  d$u <- ifelse(is.finite(d$u_pct_points), d$u_pct_points, 0)
  d$lo <- pmax(0, d$pct_of_cd45_pos - d$u)
  d$hi <- d$pct_of_cd45_pos + d$u

  ord <- stats::aggregate(pct_of_cd45_pos ~ population, d, median)
  d$population <- factor(d$population,
                         levels = ord$population[order(ord$pct_of_cd45_pos)])

  has_grp <- !is.null(group_of) && any(!is.na(group_of[d$sample_id]))
  if (has_grp) {
    d$grp <- unname(group_of[d$sample_id])
    d <- d[!is.na(d$grp), , drop = FALSE]
    fills <- population_colours(sort(unique(d$grp)), colors = colors)
  }

  # How many populations have bars wider than the between-sample spread. Stated
  # in the subtitle because it is the reading of the figure, and a reader who
  # takes only the title should take that with it.
  agg <- do.call(rbind, lapply(split(d, d$population), function(g) data.frame(
    population = g$population[1],
    spread = stats::IQR(g$pct_of_cd45_pos, na.rm = TRUE),
    u = median(g$u, na.rm = TRUE), stringsAsFactors = FALSE)))
  n_swamped <- sum(is.finite(agg$u) & agg$u > 0 & agg$u >= agg$spread, na.rm = TRUE)

  fig <- ggplot(d, aes(x = pct_of_cd45_pos, y = population)) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                   linewidth = 0.45, colour = colors$na_fill %||% "grey60") +
    (if (has_grp) geom_point(aes(fill = grp), shape = 21, size = 2, stroke = 0.3,
                             colour = colors$tile_border %||% "grey20")
     else geom_point(shape = 21, size = 2, stroke = 0.3,
                     fill = colors$reference_fill %||% "white",
                     colour = colors$tile_border %||% "grey20")) +
    (if (has_grp) scale_fill_manual(name = NULL, values = fills) else NULL) +
    labs(title = "Population abundance with gate placement uncertainty",
         subtitle = paste0(
           "bar = standard uncertainty of that sample's frequency, propagated ",
           "from its thresholds in quadrature\n",
           n_swamped, " of ", nrow(agg), " population(s) have an uncertainty as ",
           "wide as the between-sample spread"),
         x = "% of CD45+", y = NULL) +
    theme_cyto(9, colors = colors) +
    theme(panel.grid.major.y = element_blank())
  safe_ggsave(outfile, plot = fig, width = 8.5,
              height = max(3.5, 0.34 * nlevels(d$population) + 1.8), dpi = dpi,
              limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' How many samples each population is actually measurable in
#'
#' One bar per population, split by whether that population's event count in each
#' sample was enough to quantify it, enough to detect it, or neither.
#'
#' WHY A COUNT RATHER THAN A SCATTER. The limit of detection depends only on how
#' many parent-gate events a sample contributed, so it is a property of the
#' sample and identical across the populations within it. Plotting frequencies
#' against it would put every population on the same reference line and say
#' nothing extra. What varies, and what decides whether a population is worth
#' testing, is how many samples clear the limit at all.
#'
#' HOW TO READ IT. A population quantified in every sample is measurable in this
#' cohort. One that is mostly below the limit of quantification is not, and no
#' change to the gating strategy will fix it: the numerator is small because few
#' cells were acquired, so the answer is a longer acquisition or a higher
#' `--max-events-per-file`, not a different threshold. A population split between
#' the two is the dangerous case, because its group difference can be driven
#' entirely by which samples happened to clear the limit.
#'
#' @param ufreq the `frequencies` element of [run_gate_uncertainty()]
#' @param outfile Path to write the figure to.
#' @param dpi Resolution in dots per inch. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. Default `fcs_colors()`.
#' @export
fig_detection_limits <- function(ufreq, outfile, dpi = 200,
                                 colors = fcs_colors()) {
  if (is.null(ufreq) || !nrow(ufreq) || !"detection" %in% names(ufreq)) {
    log_msg("[fig] no detection limits in the uncertainty table, figure skipped")
    return(invisible(NULL))
  }
  d <- qc_pass_rows(ufreq)
  d <- d[!is.na(d$detection), , drop = FALSE]
  if (!nrow(d)) {
    log_msg("[fig] no QC-passing samples, detection limits figure skipped")
    return(invisible(NULL))
  }

  lv <- c("quantified", "detected, below LOQ", "below LOD")
  d$detection <- factor(d$detection, levels = lv)
  m <- as.data.frame(table(population = d$population, detection = d$detection),
                     stringsAsFactors = FALSE)
  names(m)[3] <- "n"
  m$detection <- factor(m$detection, levels = lv)

  # Order populations by how often they are quantifiable, worst at the top, so
  # the populations a reader must not test read first.
  q <- stats::aggregate(n ~ population, m[m$detection == "quantified", ], sum)
  tot <- stats::aggregate(n ~ population, m, sum)
  q <- merge(tot["population"], q, by = "population", all.x = TRUE)
  q$n[is.na(q$n)] <- 0
  m$population <- factor(m$population, levels = q$population[order(-q$n)])

  n_never <- sum(q$n == 0)
  n_pop <- nrow(q)
  n_samp <- length(unique(d$sample_id))

  fig <- ggplot(m, aes(x = n, y = population, fill = detection)) +
    geom_col(width = 0.7, colour = colors$tile_border %||% "grey20",
             linewidth = 0.2) +
    scale_fill_viridis_d(name = NULL, option = colors$count_viridis %||% "D",
                         direction = -1, drop = FALSE) +
    labs(title = "Populations measurable at this acquisition depth",
         subtitle = paste0(
           "events behind each frequency, against the limits of detection and ",
           "quantification set by that sample's parent gate\n",
           n_never, " of ", n_pop, " population(s) are quantifiable in no ",
           "sample; ", n_samp, " QC-passing sample(s)"),
         x = "samples", y = NULL) +
    theme_cyto(9, colors = colors) +
    theme(panel.grid.major.y = element_blank())
  safe_ggsave(outfile, plot = fig, width = 8.5,
              height = max(3.5, 0.34 * nlevels(m$population) + 1.8), dpi = dpi,
              limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' Where each population's uncertainty comes from
#'
#' Stacked contribution per term, median across samples, in percentage points.
#' The terms are the markers the population's definition reads, plus the CD45
#' parent gate, which enters every population at once because it sets the
#' denominator.
#'
#' WHY IT IS WORTH A FIGURE OF ITS OWN. It says which gate to fix. A population
#' whose budget is dominated by one marker has one threshold worth inspecting in
#' gating_qc.png; a cohort whose budgets are dominated by the parent term has a
#' CD45 gate problem that no amount of attention to the downstream markers will
#' improve. Operator studies of manual gating find the same concentration at the
#' first gate, which is the one place this figure and that literature can be read
#' against each other.
#'
#' @param budget the `budget` element of [run_gate_uncertainty()]
#' @param outfile Path to write the figure to.
#' @param dpi Resolution in dots per inch. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. Default `fcs_colors()`.
#' @export
fig_uncertainty_budget <- function(budget, outfile, dpi = 200,
                                   colors = fcs_colors()) {
  if (is.null(budget) || !nrow(budget)) {
    log_msg("[fig] no uncertainty budget, figure skipped")
    return(invisible(NULL))
  }
  d <- budget[is.finite(budget$u_pct_points), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  m <- stats::aggregate(u_pct_points ~ population + term, d, median)
  tot <- stats::aggregate(u_pct_points ~ population, m, sum)
  m$population <- factor(m$population,
                         levels = tot$population[order(tot$u_pct_points)])
  terms <- sort(unique(m$term))
  # The parent term last in the stack and in the key, so it reads as the shared
  # component it is rather than as one marker among the others.
  terms <- c(setdiff(terms, "parent_CD45"), intersect("parent_CD45", terms))
  m$term <- factor(m$term, levels = terms)

  fig <- ggplot(m, aes(x = u_pct_points, y = population, fill = term)) +
    geom_col(width = 0.7, colour = colors$tile_border %||% "grey20",
             linewidth = 0.2) +
    scale_fill_viridis_d(name = NULL, option = colors$count_viridis %||% "D",
                         direction = -1) +
    labs(title = "Uncertainty budget per population",
         subtitle = paste("median contribution of each threshold to the",
                          "population's standard uncertainty, in percentage points;",
                          "\nterms combine in quadrature, so the bar is longer",
                          "than the total it sums to"),
         x = "contribution (percentage points)", y = NULL) +
    theme_cyto(9, colors = colors) +
    theme(panel.grid.major.y = element_blank())
  safe_ggsave(outfile, plot = fig, width = 9,
              height = max(3.5, 0.34 * nlevels(m$population) + 1.8), dpi = dpi,
              limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}
