# SECTION 9c -- EXTERNAL ABSOLUTE CELL COUNTS QC FIGURE
# =============================================================================

#' QC heatmap for --absolute-counts: one tile per sample x population
#'
#' WHAT: log10(cells/uL) as tile colour, one row per matched sample, one
#' column per population, ordered by group where available. Blank tiles are
#' values genuinely absent from the source file for that sample x population,
#' not zero.
#'
#' WHY THIS COMES BEFORE THE GROUP-COMPARISON FIGURE (absolute_counts.png):
#' this is externally supplied data this pipeline did not measure -- a
#' transcription slip, a stray order-of-magnitude error, or a sample that
#' silently failed to match (see the NOTE lines the loader prints) will not
#' look anomalous in a single boxplot the way it jumps out in a full grid.
#' Controls and QC-failed samples are labelled rather than dropped, on the
#' same principle fig_gating_qc() uses: a QC figure that hides the excluded
#' rows cannot be used to catch a problem with them.
#' @param ac The ac.
#' @param outfile Path to write the figure to.
#' @param group_of Named character vector mapping sample_id to group label.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_absolute_counts_qc <- function(ac, outfile, group_of = NULL, dpi = 200, colors = fcs_colors()) {
  d <- ac
  is_ctrl <- if ("is_control" %in% names(d)) d$is_control else FALSE
  qc_stat <- if ("qc_status" %in% names(d)) d$qc_status else NA_character_
  d$tag <- ifelse(is_ctrl, " [control]",
                  ifelse(!is.na(qc_stat) & qc_stat != "pass", " [QC-failed]", ""))
  d$row_label <- paste0(d$sample_id, d$tag)
  if (!is.null(group_of)) {
    d$grp <- unname(group_of[d$sample_id])
    ord <- unique(d[order(d$grp, d$sample_id), c("row_label", "grp")])
  } else {
    ord <- unique(d[order(d$sample_id), "row_label", drop = FALSE])
  }
  d$row_label <- factor(d$row_label, levels = rev(unique(ord$row_label)))
  pop_ord <- stats::aggregate(cells_per_ul ~ population, d, median)
  d$population <- factor(d$population, levels = pop_ord$population[order(-pop_ord$cells_per_ul)])

  n_s <- nlevels(d$row_label); n_p <- nlevels(d$population)
  fig <- ggplot(d, aes(population, row_label, fill = log10(pmax(cells_per_ul, 1)))) +
    geom_tile(colour = colors$tile_border, linewidth = 0.3) +
    scale_fill_viridis_c(name = "log10(cells/\u00b5L)", na.value = colors$na_fill,
                         option = colors$count_viridis) +
    labs(title = "Absolute cell count QC, every matched sample x population",
         subtitle = paste("blank = no value for that sample in the source file;",
                          "[control]/[QC-failed] rows are shown, not excluded, for review"),
         x = NULL, y = NULL) +
    theme_cyto(9, colors = colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 7.5),
          panel.grid = element_blank())
  safe_ggsave(outfile, plot = fig, width = max(8, 0.5 * n_p + 2),
             height = max(4, 0.22 * n_s + 1.5), dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

# =============================================================================
