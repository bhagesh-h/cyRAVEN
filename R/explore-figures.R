# Figures for explore mode. All written into <outdir>/explore/.
#
# Deliberately a small set. An unsupervised run can generate a figure per marker
# per grouping without limit, and most of them are never looked at. These are the
# four that carry a decision:
#
#   clusters      what was found, and whether the clustering is cutting through
#                 continuous structure rather than separating islands
#   heatmap       what each cluster IS -- the one figure that names them
#   by group      whether the cohorts occupy different regions
#   markers       whether the map is organised by biology at all

#' Figures for explore mode
#'
#' @param cells data.frame with umap_1, umap_2, cluster, sample_id and the
#'   feature columns.
#' @param feats Feature names.
#' @param prof Cluster profile from the orchestrator, carrying `phenotype`.
#' @param ex_dir Output directory.
#' @param tag Panel suffix.
#' @param group_col Name of the group column, or NULL.
#' @return Character vector of file names written.
#' @keywords internal
explore_figures <- function(cells, feats, prof, ex_dir, tag = "",
                            group_col = NULL) {
  written <- character(0)
  keep <- function(f) written <<- c(written, basename(f))

  lab <- setNames(paste0(prof$cluster, ": ", prof$phenotype), prof$cluster)
  cells$cluster_lab <- unname(lab[cells$cluster])
  # k1, k2, ... k30 rather than k1, k10, k11, ... k2. See natural_cluster_factor().
  cells$cluster <- natural_cluster_factor(cells$cluster)
  nk <- length(levels(cells$cluster))

  # ---- 1. the embedding, coloured by cluster -------------------------------
  # THE LEGEND IS SIZED TO THE FIGURE, not left to run off it. At --explore-k 30
  # a single-column legend is 30 keys tall, which is taller than the 6.5in plot,
  # so ggplot silently clipped it and the last third of the clusters -- k7
  # onwards -- had no visible key at all. The number of columns is chosen so
  # the legend fits the panel height, and the canvas widens to hold them.
  leg_ncol <- max(1L, ceiling(nk / 22L))
  f1 <- file.path(ex_dir, sprintf("explore_umap_clusters%s.png", tag))
  p1 <- ggplot2::ggplot(cells, ggplot2::aes(umap_1, umap_2, colour = cluster)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.55, show.legend = TRUE) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 2.5, alpha = 1), ncol = leg_ncol)) +
    ggplot2::labs(title = "Explore: unsupervised clusters",
                  subtitle = paste("every eligible channel, no population",
                                   "specification used"),
                  x = "UMAP1", y = "UMAP2", colour = NULL) +
    theme_cyto() + ggplot2::theme(aspect.ratio = 1)
  # Square panel plus room for the legend beside it, rather than a wide canvas
  # the panel is stretched to fill.
  .h1 <- max(7, 0.28 * min(nk, 22L) + 2)
  safe_ggsave(f1, plot = p1, width = .h1 + 1.1 * leg_ncol + 0.8,
              height = .h1, dpi = 200, limitsize = FALSE)
  keep(f1)

  # ---- 2. what each cluster is ---------------------------------------------
  # Positivity fractions, not medians. A median is a number on a transformed
  # scale whose meaning depends on the colour limits; a fraction positive is
  # "this share of the cluster is above its own sample's cut", which is the
  # same quantity a person reads a gate for.
  fp <- prof[, grep("^frac_pos\\.", names(prof)), drop = FALSE]
  if (ncol(fp)) {
    names(fp) <- sub("^frac_pos\\.", "", names(fp))
    hm <- data.frame(
      cluster = rep(prof$cluster, ncol(fp)),
      marker = rep(names(fp), each = nrow(fp)),
      value = as.numeric(as.matrix(fp)), stringsAsFactors = FALSE)
    hm$cluster <- factor(hm$cluster, levels = rev(prof$cluster))
    f2 <- file.path(ex_dir, sprintf("explore_cluster_heatmap%s.png", tag))
    p2 <- ggplot2::ggplot(hm, ggplot2::aes(marker, cluster, fill = value)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(low = "#2c7fb8", mid = "#f7f7f7",
                                    high = "#d7301f", midpoint = 0.5,
                                    limits = c(0, 1), labels = scales::percent) +
      ggplot2::labs(title = "Explore: what each cluster is",
                    subtitle = paste("share of the cluster above that sample's",
                                     "own threshold for the marker"),
                    x = NULL, y = NULL, fill = "positive") +
      theme_cyto() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    safe_ggsave(f2, plot = p2, width = max(7, 0.45 * ncol(fp) + 4),
                height = max(4, 0.28 * nrow(prof) + 2), dpi = 200)
    keep(f2)
  }

  # ---- 2b. the same clusters as scaled median expression --------------------
  #
  # WHY BOTH THIS AND THE FRACTION-POSITIVE HEATMAP ABOVE. They answer different
  # questions and disagree usefully. A fraction positive is "what share of this
  # cluster is above its own sample's cut", which is the quantity a person reads
  # a gate for, and it is the honest one when thresholds exist. A scaled median
  # is "how bright is this cluster for this marker relative to the other
  # clusters", which needs no threshold at all -- so it still says something
  # when every threshold is a quantile fallback, and it is what an unsupervised
  # tool with no gating step (cyCONDOR's cluster_marker_heatmap.png, FlowSOM's
  # own star plots) shows. Reading one cohort through both is how you notice
  # that a cluster is uniformly dim rather than genuinely negative.
  #
  # Rows and columns are ordered by hierarchical clustering rather than left in
  # cluster-number order, because the point of the figure is which clusters
  # resemble each other, and k7 and k19 being adjacent on the axis is otherwise
  # invisible.
  md <- prof[, grep("^median\\.", names(prof)), drop = FALSE]
  if (ncol(md) >= 2L && nrow(prof) >= 2L) {
    names(md) <- sub("^median\\.", "", names(md))
    M <- as.matrix(md)
    rownames(M) <- prof$cluster
    # z-score per MARKER (column): the colour has to mean "bright for this
    # marker relative to other clusters". Scaling per cluster instead would
    # make every cluster look like it has one bright marker.
    M <- M[, apply(M, 2, function(v) is.finite(stats::sd(v)) && stats::sd(v) > 0),
           drop = FALSE]
    if (ncol(M) >= 2L) {
      Z <- scale(M)
      Z[!is.finite(Z)] <- 0
      f2b <- file.path(ex_dir, sprintf("explore_cluster_median_heatmap%s.png", tag))
      ok <- tryCatch({
        grDevices::png(f2b, width = max(1400, 46 * ncol(Z) + 700),
                       height = max(900, 34 * nrow(Z) + 500), res = 150)
        on.exit(grDevices::dev.off(), add = TRUE)
        stats::heatmap(t(Z), scale = "none",
                       col = grDevices::colorRampPalette(
                         c("#2166ac", "#f7f7f7", "#b2182b"))(64),
                       margins = c(7, 12),
                       xlab = "cluster", ylab = NULL,
                       main = "Explore: scaled median expression per cluster")
        TRUE
      }, error = function(e) FALSE)
      if (isTRUE(ok)) keep(f2b)
    }
  }

  # ---- 3. by group ----------------------------------------------------------
  if (!is.null(group_col) && group_col %in% names(cells) &&
      length(unique(stats::na.omit(cells[[group_col]]))) > 1L) {
    ng <- length(unique(stats::na.omit(cells[[group_col]])))
    f3 <- file.path(ex_dir, sprintf("explore_umap_by_group%s.png", tag))
    p3 <- ggplot2::ggplot(cells, ggplot2::aes(umap_1, umap_2, colour = cluster)) +
      ggplot2::geom_point(size = 0.2, alpha = 0.5, show.legend = FALSE) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", group_col))) +
      ggplot2::labs(title = "Explore: clusters by group",
                    subtitle = paste("one embedding, split by group -- equal",
                                     "cells per sample, so panel density is",
                                     "comparable"),
                    x = "UMAP1", y = "UMAP2") +
      theme_cyto() + theme_panel_borders() +
      ggplot2::theme(aspect.ratio = 1)
    # Height follows the rows the facets actually occupy, and aspect.ratio
    # above keeps each panel square, so the embedding has the same shape here
    # as in every other figure of it.
    .nc3 <- min(4L, ng); .nr3 <- ceiling(ng / .nc3)
    safe_ggsave(f3, plot = p3, width = 2.9 * .nc3 + 1.4,
                height = 2.9 * .nr3 + 1.4, dpi = 200, limitsize = FALSE)
    keep(f3)
  }

  # ---- 3b. each marker, split by group -------------------------------------
  # The cross of figures 3 and 4: where a marker sits differently between
  # groups, as opposed to how much of it there is. Readable only because one
  # embedding covers every sample, and honest about density only because the
  # cells were equalised per sample first.
  if (!is.null(group_col) && group_col %in% names(cells) &&
      length(unique(stats::na.omit(cells[[group_col]]))) > 1L) {
    ng <- length(unique(stats::na.omit(cells[[group_col]])))
    bg_dir <- file.path(ex_dir, "explore_marker_umaps_by_group")
    dir.create(bg_dir, showWarnings = FALSE, recursive = TRUE)
    for (m in intersect(feats, names(cells))) {
      d <- data.frame(umap_1 = cells$umap_1, umap_2 = cells$umap_2,
                      value = cells[[m]], grp = cells[[group_col]],
                      stringsAsFactors = FALSE)
      f <- file.path(bg_dir, sprintf("explore_umap_%s_by_group%s.png",
                                     gsub("[^A-Za-z0-9]+", "_", m), tag))
      .nc <- min(4L, ng); .nr <- ceiling(ng / .nc)
      p <- ggplot2::ggplot(d, ggplot2::aes(umap_1, umap_2, colour = value)) +
        ggplot2::geom_point(size = 0.18, alpha = 0.5) +
        ggplot2::facet_wrap(~ grp, ncol = .nc) +
        ggplot2::scale_colour_viridis_c(option = "C") +
        ggplot2::labs(title = paste(m, "by", group_col),
                      subtitle = paste("one embedding for every group;",
                                       "colour is expression, shared across panels"),
                      x = "UMAP1", y = "UMAP2", colour = m) +
        theme_cyto() + theme_panel_borders() +
        ggplot2::theme(aspect.ratio = 1)
      # Height was pinned at 3.9in whatever the panel count, so with four
      # groups each panel was ~1.7in wide and 2.6in tall of usable area and the
      # embedding came out letterboxed. aspect.ratio squares the panel and the
      # canvas is sized from the grid rather than fixed.
      safe_ggsave(f, plot = p, width = 2.9 * .nc + 1.8,
                  height = 2.9 * .nr + 1.4, dpi = 170, limitsize = FALSE)
      written <- c(written, file.path(basename(bg_dir), basename(f)))
    }
  }

  # ---- 4. marker expression over the embedding -----------------------------
  mf <- intersect(feats, names(cells))
  if (length(mf)) {
    long <- do.call(rbind, lapply(mf, function(m) data.frame(
      umap_1 = cells$umap_1, umap_2 = cells$umap_2, marker = m,
      value = cells[[m]], stringsAsFactors = FALSE)))
    # SQUARE-ISH GRID, not four columns. At 30 markers a 4-column grid is 8
    # rows deep, so the figure came out 2244 x 4284 -- a strip twice as tall as
    # it is wide. Scaled to fit any page or report pane that makes every panel
    # tiny, which is what "squished" looks like. ceiling(sqrt(n)) is what
    # fig_marker_grid() already uses for the declared equivalent, and it keeps
    # the whole figure near 1:1 whatever the marker count.
    ncol_grid <- max(1L, ceiling(sqrt(length(mf))))
    nrow_grid <- ceiling(length(mf) / ncol_grid)
    f4 <- file.path(ex_dir, sprintf("explore_umap_markers%s.png", tag))
    p4 <- ggplot2::ggplot(long, ggplot2::aes(umap_1, umap_2, colour = value)) +
      ggplot2::geom_point(size = 0.15, alpha = 0.5) +
      # scales = "fixed", not "free". Every panel is the SAME embedding, so a
      # free scale gave each marker its own axis range and the clouds no longer
      # lined up between panels -- the one thing this figure exists to let you
      # do. It also made each panel a slightly different shape, which is most of
      # why the grid looked squashed.
      ggplot2::facet_wrap(~ marker, ncol = ncol_grid) +
      ggplot2::scale_colour_viridis_c(option = "C") +
      ggplot2::labs(title = "Explore: marker expression over the embedding",
                    subtitle = paste("one shared embedding; every panel is the",
                                     "same cells, coloured by a different marker"),
                    x = "UMAP1", y = "UMAP2", colour = NULL) +
      theme_cyto() + theme_panel_borders() +
      # aspect.ratio = 1 pins each PANEL square regardless of how the canvas is
      # divided. Without it the panel takes whatever shape is left after the
      # strips, legend and margins, and a UMAP drawn into a wide-short box is
      # stretched horizontally -- the clusters change shape between figures of
      # the same embedding, which is exactly what must not happen.
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     aspect.ratio = 1)
    safe_ggsave(f4, plot = p4, width = 2.7 * ncol_grid + 1.6,
                height = 2.7 * nrow_grid + 1.4, dpi = 170, limitsize = FALSE)
    keep(f4)
  }

  written
}

#' Explore mode: MFI per cluster, and the composition of the embedding
#'
#' The explore counterparts of population_mfi_summary.csv and
#' umap_composition.png, built by the same functions so the two modes read
#' alike: `explore_cluster_mfi.csv` with its heatmap and dot plot, and
#' `explore_cluster_composition.csv` with `explore_umap_composition.png`.
#'
#' Per-cluster medians are taken per SAMPLE first and summarised over samples,
#' as for declared populations, so the donor with the most events does not set
#' a cluster's intensity. A sample contributes to a cluster only where it holds
#' at least ten of its cells. Percent positive is against each sample's own
#' threshold where the declared run supplied one; otherwise the "All samples"
#' row falls back to the cluster's pooled call, the basis explore_cluster_profile
#' already reports.
#'
#' @param cells Explore cell table: `sample_id`, `cluster`, `umap_1`, `umap_2`,
#'   the feature columns and optionally `group`.
#' @param feats Feature names.
#' @param prof Cluster profile, with `phenotype` and `frac_pos.*` columns.
#' @param subs Optional subset annotation from annotate_clusters_with_subsets().
#' @param tr The panel's transform.
#' @param thr_by_sample Named list: sample -> named thresholds on the analysis scale.
#' @param scatter Feature names that are log10 scatter.
#' @param ex_dir Output directory.
#' @param tag Panel suffix.
#' @return Character vector of file names written.
#' @keywords internal
explore_mfi_composition <- function(cells, feats, prof, subs = NULL, tr = NULL,
                                    thr_by_sample = list(),
                                    scatter = character(0), ex_dir, tag = "") {
  written <- character(0)
  feats <- intersect(feats, names(cells))
  if (!length(feats) || !nrow(cells)) return(written)
  plab <- if (nzchar(tag)) paste0(", ", sub("^_", "", tag)) else ""
  cells$cluster <- as.character(cells$cluster)
  ks <- levels(natural_cluster_factor(cells$cluster))
  name_of <- NULL
  if (!is.null(subs) && nrow(subs) && all(c("cluster", "subset_label") %in% names(subs)))
    name_of <- stats::setNames(as.character(subs$subset_label), subs$cluster)
  row_labels <- stats::setNames(ks, ks)
  if (!is.null(name_of)) {
    nm <- name_of[ks]
    row_labels[!is.na(nm)] <- paste0(ks[!is.na(nm)], "  ", nm[!is.na(nm)])
  }
  grp_of <- if ("group" %in% names(cells)) {
    g <- tapply(as.character(cells$group), cells$sample_id, function(x) x[1])
    g[!is.na(g)]
  } else NULL

  # ---- composition ----------------------------------------------------------
  cols <- stats::setNames(scales::hue_pal()(length(ks)), ks)
  # The marker grid comes from the profile's fraction positive, which is each
  # cluster's share above its own sample's threshold where the declared run
  # supplied thresholds and the pooled call otherwise. Scatter is left out.
  fpc <- grep("^frac_pos\\.", names(prof), value = TRUE)
  fpc <- fpc[!sub("^frac_pos\\.", "", fpc) %in% scatter]
  mpct <- if (length(fpc)) {
    M <- 100 * as.matrix(prof[, fpc, drop = FALSE])
    dimnames(M) <- list(prof$cluster, sub("^frac_pos\\.", "", fpc))
    M
  } else NULL
  comp <- composition_table(cells, "cluster",
                            phenotype = stats::setNames(prof$phenotype, prof$cluster),
                            group_col = if ("group" %in% names(cells)) "group" else NULL,
                            colours = cols, name_of = name_of, marker_pct = mpct)
  if (!is.null(comp) && nrow(comp)) {
    f <- sprintf("explore_cluster_composition%s.csv", tag)
    utils::write.csv(comp, file.path(ex_dir, f), row.names = FALSE)
    written <- c(written, f)
    fz <- file.path(ex_dir, sprintf("explore_umap_composition%s.png", tag))
    if (!is.null(fig_umap_composition(cells, comp, fz, label_col = "cluster",
          title = paste0("Explore: composition of the embedding", plab),
          subtitle = paste0(format(nrow(cells), big.mark = ","),
                            " cells after the cluster-level QC gate; markers are ",
                            "each cluster's measured phenotype, not a definition; ",
                            "on the right, the share of each cluster positive for ",
                            "every marker"),
          header = "Cluster and measured phenotype")))
      written <- c(written, basename(fz))
  }

  # ---- MFI per cluster ------------------------------------------------------
  key <- paste(cells$sample_id, cells$cluster, sep = "\r")
  idx <- split(seq_len(nrow(cells)), factor(key, levels = unique(key)))
  long <- do.call(rbind, lapply(idx, function(i) {
    if (length(i) < 10L) return(NULL)
    s <- cells$sample_id[i[1]]
    X <- as.matrix(cells[i, feats, drop = FALSE])
    thr <- thr_by_sample[[s]]
    pp <- vapply(feats, function(m) {
      t <- if (!is.null(thr) && m %in% names(thr)) as.numeric(thr[[m]]) else NA_real_
      if (is.finite(t)) 100 * mean(X[, m] > t, na.rm = TRUE) else NA_real_
    }, numeric(1))
    data.frame(sample_id = s, cluster = cells$cluster[i[1]], n_cells = length(i),
               marker = feats,
               median_transformed = apply(X, 2, stats::median, na.rm = TRUE),
               pct_positive = pp, stringsAsFactors = FALSE)
  }))
  if (is.null(long) || !nrow(long)) return(written)
  long <- long[order(match(long$cluster, ks)), , drop = FALSE]
  ms <- summarise_mfi(long, transforms = tr, group_of = grp_of,
                      label_col = "cluster", value_col = "median_transformed",
                      scatter = scatter)
  if (is.null(ms) || !nrow(ms)) return(written)
  ms$panel <- NULL
  fp <- prof[, c("cluster", grep("^frac_pos\\.", names(prof), value = TRUE)), drop = FALSE]
  miss <- which(ms$group == "All samples" & !is.finite(ms$pct_positive_median))
  for (r in miss) {
    col <- paste0("frac_pos.", ms$marker[r])
    v <- if (col %in% names(fp)) fp[[col]][match(ms$cluster[r], fp$cluster)] else NA
    if (length(v) && is.finite(v)) ms$pct_positive_median[r] <- round(100 * v, 2)
  }
  f <- sprintf("explore_cluster_mfi%s.csv", tag)
  utils::write.csv(ms, file.path(ex_dir, f), row.names = FALSE)
  written <- c(written, f)
  share <- if (!is.null(comp)) stats::setNames(comp$pct_of_embedded, comp$cluster) else NULL
  fh <- file.path(ex_dir, sprintf("explore_cluster_mfi_heatmap%s.png", tag))
  if (!is.null(fig_mfi_heatmap(ms, fh, label_col = "cluster", freq = share,
        freq_label = "% of embedded cells", cluster_rows = TRUE,
        row_labels = row_labels,
        title = paste0("Explore: median fluorescence intensity per cluster", plab),
        subtitle = paste0("Median over donors of each sample's cluster median, in ",
                          "instrument units; colour scaled per marker; rows ",
                          "ordered by similarity"))))
    written <- c(written, basename(fh))
  fd <- file.path(ex_dir, sprintf("explore_cluster_mfi_dotplot%s.png", tag))
  if (!is.null(fig_mfi_dotplot(ms, fd, label_col = "cluster", cluster_rows = TRUE,
        row_labels = row_labels,
        title = paste0("Explore: marker positivity and intensity per cluster", plab),
        subtitle = paste0("Dot size: median % positive; colour: median intensity ",
                          "scaled per marker"))))
    written <- c(written, basename(fd))
  if (!is.null(grp_of) && length(unique(grp_of)) > 1L) {
    fg <- file.path(ex_dir, sprintf("explore_cluster_mfi_heatmap_by_group%s.png", tag))
    if (!is.null(fig_mfi_heatmap(ms, fg, label_col = "cluster", by_group = TRUE,
          cluster_rows = TRUE, row_labels = row_labels,
          title = paste0("Explore: median fluorescence intensity by group", plab),
          subtitle = "One panel per group on one colour scale; values in instrument units")))
      written <- c(written, basename(fg))
  }
  written
}
