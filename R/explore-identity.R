# =============================================================================
# CLUSTER IDENTITY: what does each unsupervised cluster correspond to?
# =============================================================================
#
# THE GAP THIS CLOSES. Explore mode already writes the cross-tabulation of
# clusters against declared population labels (explore_vs_populations.csv) and a
# marker heatmap that says what each cluster EXPRESSES. Neither answers the
# question a reader actually arrives with -- "cluster 12 is what?" -- because
# the cross-tab is a few hundred rows of long-format counts and the heatmap
# names markers, not cell types. This is the figure that maps one to the other.
#
# WHAT THE LITERATURE DOES. Three representations recur, and they answer
# different questions:
#
#   1. Marker heatmap per cluster (median or scaled intensity, clusters as rows).
#      The standard annotation figure, and the one cytometry pipelines ship by
#      default -- cyCONDOR's own annotation vignette works this way, and
#      cyRAVEN already writes it as explore_cluster_median_heatmap.png. It says
#      what a cluster expresses. It does not say what it corresponds to, and a
#      reader still has to do the naming themselves.
#
#   2. Correspondence (confusion) matrix of clusters against reference
#      populations. Weber & Robinson (Cytometry A 2016, 89:1084-1096), comparing
#      clustering methods, match clusters to manually gated reference
#      populations and score the match with the F1 statistic -- the harmonic
#      mean of precision and recall -- assigning one cluster to one population
#      with the Hungarian algorithm. This is the representation that answers the
#      correspondence question directly, and it is what this figure draws.
#
#   3. Interactive annotation trees. CyCadas (Bioinformatics 2024, 40:btae595)
#      builds cell types as a tree over clusters, with everything unassigned
#      sitting at the root until a definition claims it, plus a UMAP highlight
#      per node. The idea worth borrowing from it is not the tree -- that needs
#      an interactive session -- but the INSISTENCE ON THE UNASSIGNED REMAINDER
#      being visible rather than quietly absent, which this figure follows.
#
# WHY F1 AND NOT SIMPLY THE LARGEST OVERLAP. The plurality label is the obvious
# assignment rule and it is wrong in a way that matters here. A cluster of 2,000
# cells of which 600 are CD4 T cells has CD4 T cells as its plurality even when
# those 600 are 1% of all the CD4 T cells in the run -- the label describes the
# cluster, but the cluster does not describe the label. F1 is the harmonic mean
# of both directions, so it can only be high when the cluster is mostly that
# population AND that population is mostly in this cluster, which is what
# "corresponds to" means. Both components are written to the table beside it so
# a disagreement between them is legible rather than hidden inside one number.
#
# WHY THE CATCH-ALL IS RANKED SEPARATELY. "Other CD45+" is not a cell type, it
# is the complement of the specification. Letting it win the assignment would
# label most clusters "Other" on a specification with modest coverage -- true,
# and useless. It is excluded from the ranking and reported as its own column,
# so a cluster the specification cannot describe is visibly undescribed instead
# of being given a name that means "we did not look".

#' Assign each unsupervised cluster the declared population it corresponds to
#'
#' @param ct Long cross-tabulation with columns `cluster`, `population`,
#'   `cells`, and optionally `phenotype`: one row per cluster x population pair.
#' @return data.frame, one row per cluster, carrying the called identity, the
#'   precision/recall/F1 behind it, the runner-up, and the catch-all share.
#' @keywords internal
cluster_identity_table <- function(ct) {
  if (is.null(ct) || !nrow(ct) ||
      !all(c("cluster", "population", "cells") %in% names(ct)))
    return(NULL)
  ct <- ct[is.finite(ct$cells) & ct$cells > 0, , drop = FALSE]
  if (!nrow(ct)) return(NULL)

  n_cluster <- tapply(ct$cells, ct$cluster, sum)
  n_pop     <- tapply(ct$cells, ct$population, sum)
  is_other  <- is_catch_all_label(ct$population)

  rows <- lapply(unique(as.character(ct$cluster)), function(k) {
    sub <- ct[as.character(ct$cluster) == k, , drop = FALSE]
    tot <- sum(sub$cells)
    other_pct <- 100 * sum(sub$cells[is_catch_all_label(sub$population)]) / tot
    real <- sub[!is_catch_all_label(sub$population), , drop = FALSE]
    if (!nrow(real))
      return(data.frame(
        cluster = k, cells = tot,
        identity = "undescribed", precision = NA_real_, recall = NA_real_,
        f1 = NA_real_, runner_up = NA_character_, runner_up_f1 = NA_real_,
        pct_catch_all = round(other_pct, 1),
        stringsAsFactors = FALSE))
    # precision: of this cluster, how much is that population.
    # recall:    of that population, how much is in this cluster.
    prec <- real$cells / tot
    rec  <- real$cells / as.numeric(n_pop[as.character(real$population)])
    f1   <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
    o    <- order(-f1)
    data.frame(
      cluster = k, cells = tot,
      identity = as.character(real$population)[o[1]],
      precision = round(100 * prec[o[1]], 1),
      recall    = round(100 * rec[o[1]], 1),
      f1        = round(f1[o[1]], 3),
      runner_up = if (length(o) > 1L) as.character(real$population)[o[2]] else NA_character_,
      runner_up_f1 = if (length(o) > 1L) round(f1[o[2]], 3) else NA_real_,
      pct_catch_all = round(other_pct, 1),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  # THE CALL IS QUALIFIED, NOT ASSERTED. A best match is always available -- the
  # ranking cannot fail to return something -- so the figure has to say how much
  # weight the match carries or it will be read as an annotation when it is
  # sometimes a shrug. The thresholds are deliberately blunt and stated on the
  # figure, because any cut here is a convention rather than a result.
  out$call <- ifelse(
    out$pct_catch_all >= 70, "undescribed by the specification",
    ifelse(is.na(out$f1) | out$f1 < 0.2, "no clear match",
    ifelse(out$f1 >= 0.5, "confident", "partial")))
  # THE REMAINDER CAPS THE CONFIDENCE. F1 is computed over the declared labels
  # only, so a cluster that is 60% catch-all and 40% NK cells scores a high F1
  # against NK -- the NK cells really are nearly all there -- and would be
  # called "confident NK cells" while most of the cluster is cells the
  # specification cannot name. The label is still the best available answer, but
  # "confident" is not, so a cluster more than half remainder is capped at
  # partial.
  out$call[out$pct_catch_all >= 50 & out$call == "confident"] <- "partial"
  out$identity[out$call == "undescribed by the specification"] <- "undescribed"
  if ("phenotype" %in% names(ct))
    out$phenotype <- ct$phenotype[match(out$cluster, as.character(ct$cluster))]
  out$pct_of_cells <- round(100 * out$cells / sum(as.numeric(n_cluster)), 2)
  out[order(-out$cells), , drop = FALSE]
}

#' Figure: which declared population each unsupervised cluster corresponds to
#'
#' The correspondence matrix of (2) above: clusters as rows, declared
#' populations as columns, fill is the share of the cluster. Rows are grouped by
#' the called identity so the answer is readable off the row labels alone, and
#' the catch-all keeps its own column at the right so an undescribed cluster is
#' visibly undescribed.
#'
#' @param ct Long cross-tabulation; see [cluster_identity_table()].
#' @param ident Output of [cluster_identity_table()]; computed if NULL.
#' @param outfile Destination PNG.
#' @param panel_label Panel name for the title; empty for none.
#' @param colors Palette.
#' @keywords internal
fig_cluster_identity <- function(ct, ident, outfile, panel_label = "",
                                 colors = fcs_colors()) {
  if (is.null(ident)) ident <- cluster_identity_table(ct)
  if (is.null(ident) || !nrow(ident)) return(invisible(NULL))
  d <- ct[is.finite(ct$cells) & ct$cells > 0, , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  tot <- tapply(d$cells, d$cluster, sum)
  d$pct <- 100 * d$cells / as.numeric(tot[as.character(d$cluster)])

  # Columns: real populations first, in a stable order, catch-all last. The
  # remainder is not a lineage and must not sit between two that are.
  pops <- unique(as.character(d$population))
  oth  <- pops[is_catch_all_label(pops)]
  real <- sort(setdiff(pops, oth))
  d$population <- factor(as.character(d$population), levels = c(real, sort(oth)))

  # Rows: grouped by the called identity, largest cluster first inside each
  # group, so a reader scanning the row labels reads the answer without
  # decoding the fill.
  ident <- ident[order(ident$identity == "undescribed", ident$identity,
                       -ident$cells), , drop = FALSE]
  # THE SUBSET PHENOTYPE IS IN THE ROW LABEL. Without it the figure answers only
  # "which declared population does this overlap", and the declared populations
  # are lineages -- so every answer reads "B cells" or "T cells" and the subset
  # structure the clustering actually found is invisible. The label is what the
  # cluster EXPRESSES; the facet it sits in is what it OVERLAPS.
  rlab <- sprintf("%s  (n=%s)%s", ident$cluster,
                  format(ident$cells, big.mark = ",", trim = TRUE),
                  if (!is.null(ident$subset_label))
                    ifelse(is.na(ident$subset_label), "",
                           paste0("  -  ", ident$subset_label)) else "")
  names(rlab) <- ident$cluster
  d <- d[as.character(d$cluster) %in% ident$cluster, , drop = FALSE]
  d$row <- factor(unname(rlab[as.character(d$cluster)]), levels = rev(rlab))
  d$grp <- factor(ident$identity[match(as.character(d$cluster), ident$cluster)],
                  levels = unique(ident$identity))

  nk <- nrow(ident); np <- nlevels(d$population)
  # THE WIDTH IS DECIDED HERE, BEFORE THE TEXT IS WRITTEN, because ggplot2 does
  # not wrap a title, subtitle or caption: a long one is drawn as a single line
  # and clipped at the device edge with no warning. Both were losing their ends
  # on the first draft of this figure. cap_wrap() hard-wraps each to the canvas
  # it will actually be drawn on.
  .w <- max(9, 0.46 * np + 6.2)
  .h <- max(5.5, 0.26 * nk + 3.4)
  p <- ggplot2::ggplot(d, ggplot2::aes(population, row, fill = pct)) +
    ggplot2::geom_tile(colour = colors$tile_border, linewidth = 0.25) +
    # The number is printed, not only encoded. This figure is read one row at a
    # time to answer a question about one cluster, and reading a percentage off
    # a colour bar is exactly the operation a reader should not have to do for
    # a value they are going to quote.
    # The ink flips on the dark end of the scale. A fixed dark label is
    # unreadable on the tiles that matter most here -- a cluster that is 98%
    # catch-all is exactly the row a reader needs the number from, and it is
    # also the darkest fill on the figure.
    ggplot2::geom_text(data = d[d$pct >= 5, , drop = FALSE],
                       ggplot2::aes(label = round(pct),
                                    colour = pct > 55),
                       size = 2.5, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey20",
                                            `TRUE` = "white"),
                                 guide = "none") +
    ggplot2::facet_grid(grp ~ ., scales = "free_y", space = "free_y",
                        switch = "y",
                        labeller = ggplot2::label_wrap_gen(width = 18)) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B",
                                 limits = c(0, 100), name = "% of cluster") +
    ggplot2::labs(
      title = paste0("What each unsupervised cluster corresponds to",
                     if (nzchar(panel_label)) paste0(" - ", panel_label) else ""),
      subtitle = cap_wrap(paste(
        "rows are clusters, grouped by the declared population they best match;",
        "fill is the share of the CLUSTER. The match is scored by F1 -- high",
        "only when the cluster is mostly that population AND that population is",
        "mostly in this cluster"), .w),
      caption = cap_wrap(paste(
        "Assignment is descriptive, not a gate: every cluster has a best match,",
        "so read the call in explore_cluster_identity.csv (confident F1>=0.5,",
        "partial, no clear match) before quoting a name. 'undescribed' means",
        "at least 70% of the cluster is the catch-all -- cells inside the parent",
        "gate that no declared population matched -- which is a gap in the",
        "SPECIFICATION, not a property of the cluster. Method after Weber &",
        "Robinson, Cytometry A 2016;89:1084-1096."), .w),
      x = NULL, y = NULL) +
    theme_cyto(9, colors = colors) + theme_panel_borders() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7.5),
      axis.text.y = ggplot2::element_text(size = 7),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 7.5),
      plot.caption = ggplot2::element_text(hjust = 0, size = 7),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid = ggplot2::element_blank())

  # The extra inches are for the row labels: they now carry the subset
  # phenotype as well as the cluster id, which is the longest text on the
  # figure and would otherwise squeeze the tiles.
  safe_ggsave(outfile, plot = p, width = .w + 2.2, height = .h,
              dpi = 200, limitsize = FALSE)
  invisible(outfile)
}
