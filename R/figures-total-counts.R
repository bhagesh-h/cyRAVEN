# Figures for --total-counts.
#
# Two, in the order they must be read.
#
#   1. the QC figure, which is about the INPUT: are these yields plausible, and
#      does every acquisition in the run actually have one? Everything derived
#      from a total inherits that total's errors, so a mistyped yield is
#      invisible in the derived table and obvious here. This mirrors the rule
# --absolute-counts already follows with absolute_counts_qc.png: inspect
#      the external numbers before quoting anything computed from them.
#
#   2. the derived figure, which is about the RESULT: each cluster as a cell
#      number rather than a share, so that a compartment moving as a whole is
#      visible at all.

#' Run an expression with the global RNG stream left exactly as it was found
#'
#' WHY THIS EXISTS. `run_cyraven()` seeds once and every later draw comes from
#' that one stream: the embedding's cell selection, the clustering, the
#' bootstraps. A step that spends draws therefore does not merely get its own
#' random numbers, it shifts every draw after it, and the visible symptom is
#' that adding an unrelated flag silently moves the UMAP. `geom_jitter()` draws,
#' so a figure is enough to do it -- on a two-panel explore run, drawing panel
#' 1's figure would re-roll panel 2's embedding.
#'
#' The package's contract is that a run without `--total-counts` produces
#' byte-identical output to one with it, apart from the files the flag adds.
#' Restoring the stream is what makes that true. Same idiom as
#' `clin_boot_ci()` and `run_flowsom()`.
#' @param expr Expression to evaluate.
#' @keywords internal
without_spending_draws <- function(expr) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  force(expr)
}

#' Sanity figure for externally supplied total cell counts
#'
#' WHY A LOG SCALE. Yields across a clinical cohort routinely span an order of
#' magnitude, and the failure this figure exists to catch -- a value entered in
#' the wrong unit -- lands three or six decades away. On a linear axis one such
#' value flattens every real one into the same pixel row and hides exactly the
#' problem being looked for.
#' @param tc Matched total counts from `load_total_counts()`.
#' @param outfile Destination PNG.
#' @param group_of Optional named vector, sample_id -> group.
#' @param all_samples Optional character vector of every sample in the run, so
#'   the ones with no external total can be named rather than silently absent.
#' @param panel_label Marker-panel name added to the title; empty for none.
#'   A multi-panel run writes one of these per panel, and without the label in
#'   the image the two are distinguishable only by filename -- which fails the
#'   moment either is pasted into a slide or a report. Same convention as
#'   `fig_umap_overview()` and `fig_group_comparison()`.
#' @keywords internal
fig_total_counts_qc <- function(tc, outfile, group_of = NULL,
                                all_samples = NULL, panel_label = "") {
  if (is.null(tc) || !nrow(tc)) return(invisible(NULL))
  d <- tc[order(tc$total_cells), , drop = FALSE]
  d$sample_id <- factor(d$sample_id, levels = d$sample_id)
  d$group <- if (!is.null(group_of)) unname(group_of[as.character(d$sample_id)]) else NA
  if (all(is.na(d$group))) d$group <- "all samples"

  med <- stats::median(d$total_cells)
  missing_n <- if (!is.null(all_samples))
    length(setdiff(all_samples, tc$sample_id)) else 0L

  sub <- paste0(nrow(d), " acquisition(s) with an external total; median ",
                format(med, big.mark = ",", scientific = FALSE), " cells",
                if (missing_n) paste0("; ", missing_n, " acquisition(s) in this run have NONE",
                                      " and are absent from every absolute figure") else "")

  # Points on a stem, NOT bars. A bar encodes its value as the length from a
  # baseline, and a log axis has no meaningful baseline -- drawn as columns,
  # yields an order of magnitude apart differ by a few percent of bar height
  # and the unit error this figure exists to catch becomes invisible, which is
  # the exact failure it was drawn to prevent.
  # Computed outside aes(): referring to d$ inside it is deprecated in ggplot2
  # and warns on every draw, which is noise in a log a reader is meant to scan
  # for the NOTEs that matter.
  stem_base <- min(d$total_cells) * 0.75
  p <- ggplot2::ggplot(d, ggplot2::aes(sample_id, total_cells, colour = group)) +
    ggplot2::geom_segment(ggplot2::aes(xend = sample_id, yend = total_cells),
                          y = stem_base, linewidth = 0.4, alpha = 0.55) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_hline(yintercept = med, linetype = "dashed", linewidth = 0.4) +
    ggplot2::scale_y_log10(labels = function(x)
      format(x, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    ggplot2::labs(title = paste0("External total cell counts, as supplied",
                                 if (nzchar(panel_label))
                                   paste0(" - ", panel_label) else ""),
                  subtitle = sub,
                  caption = paste("Read this before any cells_absolute number.",
                                  "A yield in the wrong unit sits decades off the",
                                  "dashed median and is invisible downstream."),
                  x = NULL, y = "total cells (log scale)", fill = NULL) +
    theme_cyto() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  safe_ggsave(outfile, plot = p, width = max(7, 0.42 * nrow(d) + 3.5),
              height = 5.2, dpi = 200)
  invisible(outfile)
}

#' Cluster abundance as a cell number, beside the same cluster as a share
#'
#' WHY BOTH PANELS AND NOT JUST THE ABSOLUTE ONE. The pair is the argument. A
#' cluster that moves on the right panel and not the left is a compartment that
#' changed size, which no share can show; one that moves on the left and not the
#' right is a redistribution at constant size. Shown apart, either panel invites
#' the wrong reading of the other.
#' @param ab Abundance table carrying a share column and cells_absolute.
#' @param outfile Destination PNG.
#' @param group_of Named vector, sample_id -> group.
#' @param max_clusters Facet cap, largest populations first.
#' @param share_col The share column. Explore writes `pct_of_gated`, the
#'   declared run writes `pct_of_cd45_pos`; resolved rather than hard-coded so
#'   one figure serves both, and so a table carrying neither is skipped rather
#'   than drawn against a column of NULL.
#' @keywords internal
fig_absolute_vs_share <- function(ab, outfile, group_of = NULL, max_clusters = 12L,
                                  share_col = NULL, panel_label = "") {
  if (is.null(ab) || !"cells_absolute" %in% names(ab)) return(invisible(NULL))
  share_col <- share_col %||%
    intersect(c("pct_of_gated", "pct_of_cd45_pos"), names(ab))[1]
  if (is.na(share_col) || !share_col %in% names(ab)) return(invisible(NULL))
  d <- ab[is.finite(ab$cells_absolute), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  d$group <- if (!is.null(group_of)) unname(group_of[d$sample_id]) else "all samples"
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))

  # Largest clusters first: a facet grid of 40 rare clusters is unreadable and
  # the rare ones are where the counting limit bites hardest anyway.
  ord <- names(sort(tapply(d[[share_col]], d$population, stats::median,
                           na.rm = TRUE), decreasing = TRUE))
  keep <- utils::head(ord, max_clusters)
  d <- d[d$population %in% keep, , drop = FALSE]
  d$population <- factor(d$population, levels = keep)

  # The share row's label carries the column it came from, so a reader can tell
  # a declared run's pct_of_cd45_pos from explore's pct_of_gated at a glance.
  # The factor levels are built FROM that label rather than written out again --
  # spelling them a second time is how the share row silently became NA and
  # dropped out of the figure.
  share_lab <- paste0("share (", share_col, ", %)")
  # Explore's rows are numbered clusters (k1, k2, ...); the declared run's are
  # named populations. Decided from the data rather than passed in, so a caller
  # cannot label one as the other by forgetting an argument.
  unit_lab <- if (all(grepl("^k?[0-9]+$", unique(as.character(d$population)))))
    "Explore clusters" else "Populations"
  long <- rbind(
    data.frame(population = d$population, group = d$group, value = d[[share_col]],
               measure = share_lab, stringsAsFactors = FALSE),
    data.frame(population = d$population, group = d$group, value = d$cells_absolute,
               measure = "cells (absolute)", stringsAsFactors = FALSE))
  long$measure <- factor(long$measure, levels = c(share_lab, "cells (absolute)"))
  stopifnot(!anyNA(long$measure))
  n_zero <- sum(!is.finite(long$value) | long$value <= 0, na.rm = TRUE)
  long <- long[is.finite(long$value) & long$value > 0, , drop = FALSE]
  if (!nrow(long)) return(invisible(NULL))

  p <- ggplot2::ggplot(long, ggplot2::aes(group, value, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.55, linewidth = 0.3) +
    # position_jitter(seed=) rather than a bare geom_jitter(): it fixes the
    # offsets so the figure is reproducible, and ggplot2 evaluates it inside a
    # temporary RNG state so the draws do not reach the shared stream. The
    # without_spending_draws() wrapper around the call site is the belt to this
    # braces, because the guarantee is version-dependent and the cost of it
    # being wrong is a silently moved embedding.
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.16, height = 0,
                                                            seed = 42L),
                        size = 0.9, alpha = 0.8) +
    ggplot2::facet_grid(measure ~ population, scales = "free_y") +
    # LOG SCALE ON BOTH ROWS. Population sizes in blood span orders of
    # magnitude -- T cells and regulatory T cells differ by a factor of ~20 in
    # share and rather more in cell number -- so on a linear axis the small
    # populations collapse onto the baseline and any change in them is
    # invisible, which is exactly where a change matters most. A log axis gives
    # equal vertical space to an equal FOLD change, which is the quantity being
    # compared between groups anyway.
    #
    # Zeroes cannot be drawn on a log axis, so they are dropped rather than
    # silently shifted by a pseudocount: a zero and a very small number are
    # different findings and a pseudocount would merge them. The subtitle says
    # how many were dropped.
    ggplot2::scale_y_log10(labels = function(x)
      format(x, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    ggplot2::labs(
      # "Clusters" only where they ARE clusters. This figure serves both the
      # declared run, whose rows are named populations from the specification,
      # and explore, whose rows are numbered clusters; calling the declared
      # populations clusters invites the reader to treat a curated gate as an
      # unsupervised guess, which is the one distinction the package is built on.
      title = paste0(unit_lab, ": share against absolute cell number",
                     if (nzchar(panel_label)) paste0(" - ", panel_label) else ""),
      # The sample count goes in the subtitle for the same reason the panel goes
      # in the title: on a multi-panel run these figures are otherwise identical
      # furniture, and "12 samples" versus "8 samples" is the quickest way to
      # see which cohort you are looking at.
      subtitle = paste0("top ", length(keep), " ", tolower(unit_lab),
                        " by median share, ", length(unique(d$sample_id)),
                        " sample(s), log scale. Absolute = share x that",
                        " acquisition's external total",
                        if (n_zero) paste0(". ", n_zero,
                          " zero value(s) omitted: a log axis cannot draw them,",
                          " and a pseudocount would make a zero and a very",
                          " small number look alike") else ""),
      caption = paste("Moves on the lower row only: the compartment changed size,",
                      "which a frequency cannot show.",
                      "Upper row only: redistribution at constant size.",
                      "Absolute values are dual-platform and carry the external",
                      "instrument's error as well as this run's."),
      x = NULL, y = NULL, fill = NULL) +
    theme_cyto() + theme_panel_borders() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none")
  safe_ggsave(outfile, plot = p, width = max(8, 1.5 * length(keep) + 2),
              height = 7.5, dpi = 200)
  invisible(outfile)
}
