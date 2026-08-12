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

  # ---- 1. the embedding, coloured by cluster -------------------------------
  f1 <- file.path(ex_dir, sprintf("explore_umap_clusters%s.png", tag))
  p1 <- ggplot2::ggplot(cells, ggplot2::aes(umap_1, umap_2, colour = cluster)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.55, show.legend = TRUE) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
    ggplot2::labs(title = "Explore: unsupervised clusters",
                  subtitle = paste("every eligible channel, no population",
                                   "specification used"),
                  x = "UMAP1", y = "UMAP2", colour = NULL) +
    theme_cyto()
  safe_ggsave(f1, plot = p1, width = 9, height = 6.5, dpi = 200)
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
      theme_cyto()
    safe_ggsave(f3, plot = p3, width = max(7, 3.2 * ng + 1), height = 4.2,
                dpi = 200)
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
      p <- ggplot2::ggplot(d, ggplot2::aes(umap_1, umap_2, colour = value)) +
        ggplot2::geom_point(size = 0.18, alpha = 0.5) +
        ggplot2::facet_wrap(~ grp) +
        ggplot2::scale_colour_viridis_c(option = "C") +
        ggplot2::labs(title = paste(m, "by", group_col),
                      subtitle = paste("one embedding, equal cells per sample;",
                                       "colour is expression, shared across panels"),
                      x = "UMAP1", y = "UMAP2", colour = m) +
        theme_cyto()
      safe_ggsave(f, plot = p, width = max(7, 3.1 * ng + 1.4), height = 3.9,
                  dpi = 170)
      written <- c(written, file.path(basename(bg_dir), basename(f)))
    }
  }

  # ---- 4. marker expression over the embedding -----------------------------
  mf <- intersect(feats, names(cells))
  if (length(mf)) {
    long <- do.call(rbind, lapply(mf, function(m) data.frame(
      umap_1 = cells$umap_1, umap_2 = cells$umap_2, marker = m,
      value = cells[[m]], stringsAsFactors = FALSE)))
    ncol_grid <- min(4L, length(mf))
    f4 <- file.path(ex_dir, sprintf("explore_umap_markers%s.png", tag))
    p4 <- ggplot2::ggplot(long, ggplot2::aes(umap_1, umap_2, colour = value)) +
      ggplot2::geom_point(size = 0.15, alpha = 0.5) +
      ggplot2::facet_wrap(~ marker, ncol = ncol_grid, scales = "free") +
      ggplot2::scale_colour_viridis_c(option = "C") +
      ggplot2::labs(title = "Explore: marker expression over the embedding",
                    x = "UMAP1", y = "UMAP2", colour = NULL) +
      theme_cyto() +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank())
    safe_ggsave(f4, plot = p4, width = 3.1 * ncol_grid + 1,
                height = 2.8 * ceiling(length(mf) / ncol_grid) + 1, dpi = 170)
    keep(f4)
  }

  written
}
