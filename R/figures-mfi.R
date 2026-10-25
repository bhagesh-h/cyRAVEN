# SECTION 7m -- MARKER INTENSITY (MFI) AND COMPOSITION ON THE EMBEDDING
# =============================================================================
#
# Two summaries a cytometry report is expected to carry, which this package
# produced only in pieces.
#
# MFI. population_marker_mfi.csv has always held one median per sample,
# population and marker, on the analysis scale. That is the right input for the
# differential-state tests, and the wrong thing to hand a reader: a median of
# 2.31 on a logicle scale is not a number anyone compares with a published
# value. The summary below gives the cohort-level median fluorescence intensity
# in the instrument's own units, with its spread across donors, per study group,
# and the two figures that make the table readable at a glance:
#
#   heatmap   the value itself printed in every tile, the tile coloured on the
#             analysis scale rescaled per marker between its 1st and 99th
#             percentile -- the CATALYST convention (plotExprHeatmap), so the
#             colour means "bright for this marker" and the number means "how
#             bright". A bar beside each row gives the population's share, as
#             CATALYST's frequency annotation does, because a bright population
#             of three cells is not the same finding as a bright population of
#             thirty percent.
#   dot plot  the Seurat DotPlot layout: dot size is the percent of the
#             population above its own sample's threshold, colour is the same
#             per-marker scaled intensity. Size and colour disagreeing is the
#             informative case -- a population that is uniformly dim reads as
#             a large pale dot, one with a bright minority as a small dark one.
#
# WHY MEDIANS OF MEDIANS, BACK-TRANSFORMED. Donors are the replicates (the
# aggregation strategy of Weber et al. 2019), so the cohort value is the median
# of the per-sample medians, not the median of pooled cells, which would let the
# donor with the most events set it. The median is computed on the analysis
# scale and then inverted through the run's own transform. Every transform here
# is monotone, and a median commutes with a monotone function, so this IS the
# median in instrument units -- not an approximation of it -- while avoiding a
# second pass over the raw events, which the pipeline frees after scoring.
#
# COMPOSITION. umap_populations.png already prints each population's share in
# its key. What it cannot show is what each colour MEANS: the key gives a name,
# and the definition behind the name lives in a YAML file the reader may never
# have seen. The composition figure puts both beside the map -- name, defining
# markers, share of the embedded cells, and the share within each study group --
# so the figure can be read on its own. The table under it carries the same
# numbers with the per-donor spread.

#' The defining markers of each declared population, as a short string
#'
#' "CD3+ CD4+ CD127- CD25+" from `list(CD3 = "above", CD4 = "above", ...)`.
#' Intermediate and bright requirements read `CD14int` and `CD56hi`, and an
#' `any_of` disjunction is written in brackets: `(CD16+ or CD56+)`.
#' @param spec Population specification, as [default_population_spec()].
#' @return Named character vector, one string per population.
#' @keywords internal
population_phenotype <- function(spec) {
  if (!length(spec)) return(character(0))
  one <- function(mk, d) switch(as.character(d),
    above = paste0(mk, "+"), below = paste0(mk, "-"),
    intermediate = paste0(mk, "int"), bright = paste0(mk, "hi"),
    paste0(mk, " ", d))
  terms <- function(req) {
    if (!length(req)) return(character(0))
    unname(unlist(Map(one, names(req),
                      vapply(req, function(v) as.character(v)[1], character(1)))))
  }
  vapply(spec, function(req) {
    s <- terms(req[setdiff(names(req), "any_of")])
    a <- req[["any_of"]]
    if (length(a)) s <- c(s, paste0("(", paste(terms(a), collapse = " or "), ")"))
    paste(s, collapse = " ")
  }, character(1))
}

#' Analysis-scale values back to instrument units
#'
#' Inverts the run's transform marker by marker. A scatter channel that explore
#' mode put on a log10 scale is inverted as such; a marker the transform has no
#' parameters for comes back NA rather than through a wrong inverse.
#' @param y Numeric vector on the analysis scale.
#' @param marker Character vector, same length, naming each value's marker.
#' @param tr A transform from [make_transform()].
#' @param scatter Markers that are log10 scatter rather than fluorescence.
#' @return Numeric vector in instrument units.
#' @keywords internal
mfi_back_transform <- function(y, marker, tr, scatter = character(0)) {
  out <- rep(NA_real_, length(y))
  if (is.null(tr) || !length(y)) return(out)
  marker <- as.character(marker)
  for (m in unique(marker)) {
    i <- which(marker == m & is.finite(y))
    if (!length(i)) next
    out[i] <- if (m %in% scatter) 10^y[i]
              else tryCatch(tr$inv(y[i], m), error = function(e) NA_real_)
  }
  out
}

#' Compact intensity labels for heatmap tiles: 840, 12.3k, 1.2M
#' @param x Numeric vector.
#' @keywords internal
format_mfi <- function(x) {
  out <- rep("", length(x))
  ok <- is.finite(x)
  a <- abs(x[ok])
  out[ok] <- ifelse(a < 1e3, sprintf("%.0f", x[ok]),
             ifelse(a < 1e5, sprintf("%.1fk", x[ok] / 1e3),
             ifelse(a < 1e6, sprintf("%.0fk", x[ok] / 1e3),
                    sprintf("%.1fM", x[ok] / 1e6))))
  out
}

#' Cohort-level MFI per population and marker, overall and per study group
#'
#' One row per panel, population (or cluster), marker and group, where the group
#' is "All samples" plus each study group when there are two or more. Controls
#' and staining-QC failures are dropped first, as for every group comparison.
#'
#' @param mfi Per-sample table: `sample_id`, the label column, `marker`, the
#'   analysis-scale median, and optionally `panel`, `pct_positive`, `n_cells`.
#' @param transforms Either one transform or a list of them named by panel.
#' @param group_of Optional named vector, sample_id -> group.
#' @param label_col Name of the population or cluster column.
#' @param value_col Name of the analysis-scale median column.
#' @param pos_col Name of the percent-positive column, if any.
#' @param count_col Name of the cell-count column, if any.
#' @param scatter Markers on a log10 scatter scale. See [mfi_back_transform()].
#' @return data.frame, or NULL when there is nothing to summarise.
#' @keywords internal
summarise_mfi <- function(mfi, transforms = NULL, group_of = NULL,
                          label_col = "population", value_col = "median_asinh",
                          pos_col = "pct_positive", count_col = "n_cells",
                          scatter = character(0)) {
  if (is.null(mfi) || !nrow(mfi)) return(NULL)
  d <- qc_pass_rows(mfi)
  d <- d[is.finite(d[[value_col]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  if (!"panel" %in% names(d)) d$panel <- ""
  lab_lv <- unique(as.character(d[[label_col]]))
  mk_lv <- unique(as.character(d$marker))
  d$group <- "All samples"
  parts <- list(d)
  grp_lv <- character(0)
  if (!is.null(group_of)) {
    g <- unname(group_of[as.character(d$sample_id)])
    ok <- !is.na(g) & nzchar(g)
    if (length(unique(g[ok])) > 1L) {
      dg <- d[ok, , drop = FALSE]
      dg$group <- g[ok]
      parts <- c(parts, list(dg))
      grp_lv <- sort(unique(g[ok]))
    }
  }
  d <- do.call(rbind, parts)
  key <- paste(d$panel, d[[label_col]], d$marker, d$group, sep = "\r")
  idx <- split(seq_len(nrow(d)), factor(key, levels = unique(key)))
  has_pos <- pos_col %in% names(d)
  has_n <- count_col %in% names(d)
  out <- do.call(rbind, lapply(idx, function(i) {
    q <- stats::quantile(d[[value_col]][i], c(0.25, 0.5, 0.75), names = FALSE)
    data.frame(panel = d$panel[i[1]], label = d[[label_col]][i[1]],
               marker = d$marker[i[1]], group = d$group[i[1]],
               n_samples = length(unique(d$sample_id[i])),
               median_transformed = q[2], q1_transformed = q[1],
               q3_transformed = q[3],
               pct_positive_median = if (has_pos)
                 stats::median(d[[pos_col]][i], na.rm = TRUE) else NA_real_,
               cells_median = if (has_n)
                 stats::median(d[[count_col]][i], na.rm = TRUE) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out$mfi <- out$mfi_q1 <- out$mfi_q3 <- NA_real_
  out$transform <- NA_character_
  for (pn in unique(out$panel)) {
    tr <- if (is.function(transforms$inv)) transforms else transforms[[pn]]
    if (is.null(tr)) next
    i <- which(out$panel == pn)
    out$mfi[i]    <- mfi_back_transform(out$median_transformed[i], out$marker[i], tr, scatter)
    out$mfi_q1[i] <- mfi_back_transform(out$q1_transformed[i], out$marker[i], tr, scatter)
    out$mfi_q3[i] <- mfi_back_transform(out$q3_transformed[i], out$marker[i], tr, scatter)
    out$transform[i] <- ifelse(out$marker[i] %in% scatter, "log10 scatter",
                               tr$method %||% NA_character_)
  }
  out <- out[order(out$panel, match(out$group, c("All samples", grp_lv)),
                   match(out$label, lab_lv), match(out$marker, mk_lv)), , drop = FALSE]
  out <- data.frame(
    panel = out$panel, label = out$label, marker = out$marker, group = out$group,
    n_samples = out$n_samples,
    mfi = signif(out$mfi, 5), mfi_q1 = signif(out$mfi_q1, 5),
    mfi_q3 = signif(out$mfi_q3, 5),
    median_transformed = round(out$median_transformed, 4),
    q1_transformed = round(out$q1_transformed, 4),
    q3_transformed = round(out$q3_transformed, 4),
    pct_positive_median = round(out$pct_positive_median, 2),
    cells_median = out$cells_median, transform = out$transform,
    stringsAsFactors = FALSE)
  names(out)[names(out) == "label"] <- label_col
  out
}

#' Wrap a subtitle to the width of the figure it sits on
#'
#' A narrow heatmap of five markers is six inches wide, and a one-line subtitle
#' that explains the colour scale is longer than that; unwrapped, ggplot2 runs
#' it off the edge of the image.
#' @param x Subtitle text, or NULL.
#' @param width_in Figure width in inches.
#' @keywords internal
wrap_to <- function(x, width_in) {
  if (is.null(x) || !nzchar(x)) return(x)
  paste(strwrap(x, width = max(40L, floor(width_in * 13))), collapse = "\n")
}

#' Per-marker 0-1 rescaling of analysis-scale medians, for colour
#'
#' Between the 1st and 99th percentile of the values being drawn, per marker,
#' as CATALYST scales its expression heatmaps. Scaling over every group at once
#' keeps one colour meaning one intensity across the faceted figure.
#' @param d Rows of [summarise_mfi()] output for one panel.
#' @keywords internal
mfi_scale <- function(d) {
  s <- rep(NA_real_, nrow(d))
  for (m in unique(d$marker)) {
    i <- which(d$marker == m)
    v <- d$median_transformed[i]
    r <- stats::quantile(v, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    s[i] <- if (is.finite(r[2] - r[1]) && r[2] > r[1])
      pmin(1, pmax(0, (v - r[1]) / (r[2] - r[1]))) else 0.5
  }
  s
}

#' Row and column order for the MFI figures
#' @param d Rows for the "All samples" group of one panel.
#' @param label_col Label column.
#' @param cluster_rows Order rows by hierarchical clustering; otherwise keep
#'   the order they arrive in, which for declared populations is the order the
#'   specification lists them.
#' @keywords internal
mfi_axis_order <- function(d, label_col, cluster_rows = FALSE) {
  labs <- unique(as.character(d[[label_col]]))
  mks <- unique(as.character(d$marker))
  M <- matrix(NA_real_, length(labs), length(mks), dimnames = list(labs, mks))
  M[cbind(match(d[[label_col]], labs), match(d$marker, mks))] <- d$scaled
  M[!is.finite(M)] <- 0
  ord <- function(X) if (nrow(X) >= 3L)
    stats::hclust(stats::dist(X), method = "average")$order else seq_len(nrow(X))
  list(rows = if (cluster_rows) labs[ord(M)] else labs,
       cols = mks[ord(t(M))])
}

#' Heatmap of median fluorescence intensity per population and marker
#'
#' @param summ Output of [summarise_mfi()], restricted to one panel.
#' @param outfile Destination PNG.
#' @param label_col Population or cluster column in `summ`.
#' @param by_group Facet by study group instead of drawing "All samples".
#' @param freq Optional named numeric, label -> percent, drawn as a bar beside
#'   each row. Ignored when `by_group`.
#' @param freq_label Axis title for that bar.
#' @param cluster_rows See [mfi_axis_order()].
#' @param row_labels Optional named character, label -> text shown on the axis.
#' @param title,subtitle Figure title and subtitle.
#' @param dpi Resolution.
#' @return The plot, invisibly, or NULL when there was nothing to draw.
#' @keywords internal
fig_mfi_heatmap <- function(summ, outfile, label_col = "population",
                            by_group = FALSE, freq = NULL,
                            freq_label = "% of cells", cluster_rows = FALSE,
                            row_labels = NULL,
                            title = "Median fluorescence intensity",
                            subtitle = NULL, dpi = 200) {
  if (is.null(summ) || !nrow(summ)) return(invisible(NULL))
  summ$scaled <- mfi_scale(summ)
  base <- summ[summ$group == "All samples", , drop = FALSE]
  if (!nrow(base)) return(invisible(NULL))
  ax <- mfi_axis_order(base, label_col, cluster_rows)
  d <- if (by_group) summ[summ$group != "All samples", , drop = FALSE] else base
  if (!nrow(d)) return(invisible(NULL))
  ng <- length(unique(d$group))
  ylab <- function(x) if (is.null(row_labels)) x else
    ifelse(is.na(row_labels[x]), x, row_labels[x])
  d$lab <- factor(as.character(d[[label_col]]), levels = rev(ax$rows))
  d$marker <- factor(d$marker, levels = ax$cols)
  d$txt <- format_mfi(d$mfi)
  d$ink <- ifelse(d$scaled > 0.6, "dark", "light")
  if (by_group) d$group <- factor(d$group, levels = sort(unique(d$group)))
  p <- ggplot(d, aes(marker, lab, fill = scaled)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = txt, colour = ink), size = 2.5) +
    scale_colour_manual(values = c(dark = "grey10", light = "white"), guide = "none") +
    scale_fill_gradientn(colours = viridisLite::viridis(64), limits = c(0, 1),
                         breaks = c(0, 0.5, 1), labels = c("low", "", "high"),
                         name = "per marker") +
    scale_y_discrete(labels = ylab) +
    labs(x = NULL, y = NULL) +
    theme_cyto() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(), legend.position = "bottom",
          legend.key.height = grid::unit(0.35, "lines"),
          legend.key.width = grid::unit(2.2, "lines"),
          legend.title = element_text(size = 8.5))
  p <- p + legend_title_on_top()
  if (by_group) p <- p + facet_wrap(~ group, nrow = 1) + theme_panel_borders()
  nr <- length(ax$rows); nc <- length(ax$cols)
  lab_in <- 0.075 * max(nchar(ylab(ax$rows)), 6)
  w <- lab_in + 0.5 * nc * (if (by_group) ng else 1) + 1.4
  fig <- p
  if (!by_group && length(freq)) {
    fd <- data.frame(lab = factor(ax$rows, levels = rev(ax$rows)),
                     pct = as.numeric(freq[ax$rows]))
    fd$pct[!is.finite(fd$pct)] <- 0
    pb <- ggplot(fd, aes(pct, lab)) +
      geom_col(fill = "grey62", width = 0.7) +
      geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.12, size = 2.5) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
      labs(x = freq_label, y = NULL) + theme_cyto() +
      # No y ticks: theme_cyto() draws them inward, and beside a bar of zero
      # length the tick reads as a minus sign on its label ("-0.0%").
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.text.x = element_text(size = 7),
            axis.title.x = element_text(size = 8.5))
    fig <- patchwork::wrap_plots(p, pb, widths = c(nc, 3))
    w <- w + 1.9
  }
  w <- max(6.5, w)
  fig <- fig + patchwork::plot_annotation(
    title = title, subtitle = wrap_to(subtitle, w),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(colour = fcs_colors()$subtitle_text,
                                               size = 9)))
  safe_ggsave(outfile, plot = fig, width = w,
              height = max(4.2, 0.3 * nr + 2.8), dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' A horizontal colour bar's title above the bar rather than beside it
#'
#' Beside it, the title runs into the bar's first label. ggplot2 3.5 moved the
#' setting from the guide to the theme and warns on the old argument, so the
#' version decides which one is used.
#' @keywords internal
legend_title_on_top <- function() {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.title.position = "top")
  else guides(fill = guide_colourbar(title.position = "top"))
}

#' Dot plot: percent positive as size, scaled intensity as colour
#'
#' @inheritParams fig_mfi_heatmap
#' @return The plot, invisibly, or NULL when there was nothing to draw.
#' @keywords internal
fig_mfi_dotplot <- function(summ, outfile, label_col = "population",
                            cluster_rows = FALSE, row_labels = NULL,
                            title = "Marker positivity and intensity",
                            subtitle = NULL, dpi = 200) {
  if (is.null(summ) || !nrow(summ)) return(invisible(NULL))
  summ$scaled <- mfi_scale(summ)
  d <- summ[summ$group == "All samples", , drop = FALSE]
  if (!nrow(d) || all(!is.finite(d$pct_positive_median))) return(invisible(NULL))
  ax <- mfi_axis_order(d, label_col, cluster_rows)
  ylab <- function(x) if (is.null(row_labels)) x else
    ifelse(is.na(row_labels[x]), x, row_labels[x])
  d$lab <- factor(as.character(d[[label_col]]), levels = rev(ax$rows))
  d$marker <- factor(d$marker, levels = ax$cols)
  nr <- length(ax$rows); nc <- length(ax$cols)
  w <- max(6.5, 0.075 * max(nchar(ylab(ax$rows)), 6) + 0.42 * nc + 2.6)
  subtitle <- wrap_to(subtitle, w)
  p <- ggplot(d, aes(marker, lab)) +
    geom_point(aes(size = pct_positive_median, colour = scaled)) +
    scale_size_area(max_size = 6.5, limits = c(0, 100),
                    breaks = c(10, 25, 50, 75, 100), name = "% positive") +
    scale_colour_gradientn(colours = viridisLite::viridis(64), limits = c(0, 1),
                           breaks = c(0, 0.5, 1), labels = c("low", "", "high"),
                           name = "MFI, per marker") +
    scale_y_discrete(labels = ylab) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_cyto() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25))
  safe_ggsave(outfile, plot = p, width = w,
              height = max(4.2, 0.3 * nr + 2.4 + 0.15 * lengths(gregexpr("\n", subtitle))),
              dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(p)
}

#' Share of the embedded cells per label, with its spread across donors
#'
#' Each embedded cell carries ONE label, its most specific population (or its
#' cluster), so the shares sum to 100. That is why the numbers differ from
#' population_frequencies.csv, where the gates overlap -- "CD3 pos" there
#' includes the CD4 and CD8 T cells, and here it holds only the CD3+ cells no
#' deeper population claimed. Under `--subsample rare` the cells carry a
#' `sampling_weight`, and the shares are weighted back to the true composition.
#'
#' @param cells Embedded cell table with `sample_id` and the label column.
#' @param label_col Label column.
#' @param phenotype Optional named character, label -> defining markers.
#' @param group_col Optional study-group column in `cells`.
#' @param colours Optional named character, label -> colour, carried into the
#'   table so a row can be matched to the figure.
#' @param name_of Optional named character, label -> a second name (explore
#'   mode's matched subset).
#' @return data.frame, largest share first and the catch-all last.
#' @keywords internal
composition_table <- function(cells, label_col = "population_label",
                              phenotype = NULL, group_col = NULL,
                              colours = NULL, name_of = NULL) {
  if (is.null(cells) || !nrow(cells) || !label_col %in% names(cells)) return(NULL)
  d <- cells[!is.na(cells[[label_col]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  lab <- as.character(d[[label_col]])
  w <- if ("sampling_weight" %in% names(d) && all(is.finite(d$sampling_weight)))
    d$sampling_weight else rep(1, nrow(d))
  tot <- tapply(w, lab, sum)
  lv <- names(tot)
  m <- tapply(w, list(as.character(d$sample_id), factor(lab, levels = lv)), sum)
  m[is.na(m)] <- 0
  sh <- 100 * m / rowSums(m)
  qs <- function(p) apply(sh, 2, stats::quantile, p, names = FALSE)
  out <- data.frame(label = lv,
                    name = if (is.null(name_of)) NA_character_
                           else unname(name_of[lv]),
                    phenotype = if (is.null(phenotype)) NA_character_
                                else unname(phenotype[lv]),
                    colour = if (is.null(colours)) NA_character_
                             else unname(colours[lv]),
                    cells = as.integer(table(factor(lab, levels = lv))),
                    pct_of_embedded = round(100 * as.numeric(tot) / sum(w), 3),
                    pct_median_per_sample = round(qs(0.5), 3),
                    pct_q1_per_sample = round(qs(0.25), 3),
                    pct_q3_per_sample = round(qs(0.75), 3),
                    n_samples_present = as.integer(colSums(m > 0)),
                    n_samples = nrow(m),
                    stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(name_of)) out$name <- NULL
  if (!is.null(group_col) && group_col %in% names(d)) {
    g <- as.character(d[[group_col]])
    ok <- !is.na(g) & nzchar(g)
    gl <- sort(unique(g[ok]))
    if (length(gl) > 1L) for (G in gl) {
      i <- ok & g == G
      tg <- tapply(w[i], factor(lab[i], levels = lv), sum)
      tg[is.na(tg)] <- 0
      out[[paste0("pct_", G)]] <- round(100 * as.numeric(tg) / sum(w[i]), 3)
    }
  }
  out <- out[order(is_catch_all_label(out$label), -out$pct_of_embedded), , drop = FALSE]
  rownames(out) <- NULL
  names(out)[names(out) == "label"] <- label_col
  out
}

#' The shared embedding coloured by label, with a legend of composition
#'
#' Beside the map, one row per label: its colour, its name, its defining
#' markers, a bar of its share of the embedded cells, and its share within each
#' study group (up to four; the table carries any number). The catch-all is
#' drawn in grey beneath everything else and listed last.
#'
#' @param cells Embedded cell table with `umap_1`, `umap_2` and the label column.
#' @param comp Output of [composition_table()], with a `colour` column.
#' @param outfile Destination PNG.
#' @param label_col Label column.
#' @param title,subtitle Figure title and subtitle.
#' @param header Heading over the legend's name column.
#' @param dpi Resolution.
#' @return The plot, invisibly, or NULL when there was nothing to draw.
#' @keywords internal
fig_umap_composition <- function(cells, comp, outfile,
                                 label_col = "population_label",
                                 title = "Composition of the shared embedding",
                                 subtitle = NULL,
                                 header = "Cell type and defining markers",
                                 dpi = 220) {
  if (is.null(comp) || !nrow(comp) || is.null(cells) ||
      !all(c("umap_1", "umap_2", label_col) %in% names(cells)))
    return(invisible(NULL))
  pal <- stats::setNames(comp$colour, comp[[label_col]])
  pal[is.na(pal)] <- "grey70"
  d <- cells[!is.na(cells[[label_col]]), c("umap_1", "umap_2", label_col), drop = FALSE]
  names(d)[3] <- "lab"
  d$lab <- as.character(d$lab)
  # Draw order without spending random draws: the catch-all first, then the
  # largest populations, so the smallest land on top and stay visible.
  rank <- match(d$lab, rev(comp[[label_col]]))
  d <- d[order(rank), , drop = FALSE]
  d$lab <- factor(d$lab, levels = comp[[label_col]])
  a <- auto_point_aes(nrow(d))
  pm <- ggplot(d, aes(umap_1, umap_2, colour = lab)) +
    geom_point(size = a$size * 1.3, alpha = min(1, a$alpha + 0.25), stroke = 0) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_equal() + theme_cyto()

  # ---- the legend, drawn as a table ----------------------------------------
  # In character units on x, so the column positions follow the longest name
  # rather than being guessed; the panel is then given that many characters'
  # width in the layout.
  nr <- nrow(comp)
  name2 <- if ("name" %in% names(comp)) comp$name else rep(NA_character_, nr)
  line1 <- ifelse(is.na(name2) | !nzchar(name2), comp[[label_col]],
                  paste0(comp[[label_col]], "  ", name2))
  ph <- if ("phenotype" %in% names(comp)) comp$phenotype else rep("", nr)
  ph[is.na(ph)] <- ""
  ph <- ifelse(nchar(ph) > 56, paste0(substr(ph, 1, 55), "..."), ph)
  gcols <- grep("^pct_", setdiff(names(comp), c("pct_of_embedded",
    "pct_median_per_sample", "pct_q1_per_sample", "pct_q3_per_sample")), value = TRUE)
  if (length(gcols) > 4L) gcols <- character(0)
  tw <- max(nchar(line1) * 1.05, nchar(ph) * 0.85, nchar(header) * 0.95, 14) + 1.5
  x_sw <- 0.9; x_tx <- 2.7; x_bar <- x_tx + tw; bar_w <- 15
  x_pct <- x_bar + bar_w + 0.7
  # A group column is as wide as its name needs, wrapped onto two lines: study
  # groups are routinely called "GvHD grade 1" and "GvHD grade 3", and cut to a
  # fixed width the two headers read identically.
  gname <- vapply(sub("^pct_", "", gcols), function(g)
    paste(strwrap(g, width = max(8L, ceiling(nchar(g) / 2) + 2L)), collapse = "\n"),
    character(1))
  gw <- vapply(gname, function(g) max(7.5, max(nchar(strsplit(g, "\n")[[1]])) + 2),
               numeric(1))
  x_g <- x_pct + 7.5 + c(0, cumsum(gw))[seq_along(gcols)]
  W <- x_pct + 7.5 + sum(gw) + 0.5
  rows <- data.frame(y = rev(seq_len(nr)), line1 = line1, ph = ph,
                     pct = comp$pct_of_embedded, col = unname(pal[comp[[label_col]]]),
                     stringsAsFactors = FALSE)
  mx <- max(rows$pct, 1e-9)
  pl <- ggplot(rows) +
    geom_point(aes(x = x_sw, y = y), shape = 15, size = 5.2, colour = rows$col) +
    geom_text(aes(x = x_tx, y = y + 0.17, label = line1), hjust = 0, size = 3.2,
              fontface = "bold", colour = "grey10") +
    geom_text(aes(x = x_tx, y = y - 0.21, label = ph), hjust = 0, size = 2.6,
              colour = "grey35") +
    geom_rect(aes(xmin = x_bar, xmax = x_bar + bar_w * pct / mx,
                  ymin = y - 0.26, ymax = y + 0.26), fill = rows$col) +
    geom_text(aes(x = x_pct, y = y, label = sprintf("%.1f%%", pct)), hjust = 0,
              size = 3, colour = "grey10") +
    annotate("text", x = x_tx, y = nr + 0.7, label = header, hjust = 0,
             vjust = 0, size = 3, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_bar, y = nr + 0.7, label = "Share of embedded cells",
             hjust = 0, vjust = 0, size = 3, fontface = "bold", colour = "grey30")
  for (k in seq_along(gcols)) {
    gv <- comp[[gcols[k]]]
    pl <- pl +
      annotate("text", x = x_g[k], y = rows$y, label = sprintf("%.1f%%", gv),
               hjust = 0, size = 2.8, colour = "grey20") +
      annotate("text", x = x_g[k], y = nr + 0.7, label = gname[[k]], hjust = 0,
               vjust = 0, size = 3, fontface = "bold", colour = "grey30",
               lineheight = 0.9)
  }
  h <- max(6.8, 0.42 * nr + 1.7)
  # Rows stack from the top at a fixed pitch rather than spreading over the
  # whole height, so a legend of five rows beside a tall map stays a compact
  # block instead of five lines an inch apart.
  y_hi <- nr + 1.75
  y_lo <- min(0.45, y_hi - (h - 1.1) / 0.42)
  pl <- pl +
    scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_lo, y_hi), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA))

  side <- max(6.2, min(10, h - 1.2))
  leg_in <- W * 0.078
  fig <- patchwork::wrap_plots(pm, pl, widths = c(side, leg_in)) +
    patchwork::plot_annotation(
      title = title, subtitle = wrap_to(subtitle, side + leg_in + 0.9),
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(colour = fcs_colors()$subtitle_text,
                                                 size = 9.5),
                    plot.background = element_rect(fill = "white", colour = NA)))
  safe_ggsave(outfile, plot = fig, width = side + leg_in + 0.9, height = h,
              dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}
