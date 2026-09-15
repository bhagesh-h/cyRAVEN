# SECTION 7 -- FIGURES
# (validated separately on synthetic 3-population data before integration)
# =============================================================================

# --- shared plot theme -------------------------------------------------------
#' Publication-oriented minimal theme
#' WHY: consistent, legible defaults across every emitted figure; outward ticks
#' and frameless legends read better in print than ggplot2 defaults.
#' @param base_size Base font size in points. Default `11`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
theme_cyto <- function(base_size = 11, colors = fcs_colors()) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = colors$grid_major),
      axis.ticks       = element_line(linewidth = 0.3, colour = colors$axis_ticks),
      axis.ticks.length = unit(-0.12, "lines"),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      # An explicit white ground on the plot object itself, so a figure is
      # readable however it is rendered rather than only when safe_ggsave()
      # writes it. theme_minimal() leaves this transparent.
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      plot.title       = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle    = element_text(colour = colors$subtitle_text, size = rel(0.85))
    )
}

#' Panel borders and gutters for any faceted figure
#'
#' WHY THIS IS NOT IN `theme_cyto()`: putting a border on every panel would box
#' the single-panel figures too, where there is nothing to divide and the frame
#' is just ink. This is for figures whose panels are meant to be COMPARED.
#'
#' WHY A BORDER AT ALL. Two scatters side by side with nothing between them read
#' as one wide cloud: the eye has to find the boundary from a gap in the point
#' density, and where a facet's cells reach the edge of their panel there is no
#' gap to find. A thin border and a wider gutter make the division explicit at a
#' glance. `umap_<marker>_by_<category>.png` has always done this; every other
#' multi-panel figure in the package now does it through this one function, so a
#' report cannot contain framed and unframed versions of the same kind of figure.
#' @param colors Named list of colours; defaults to the package palette. Default `fcs_colors()`.
#' @param spacing Gutter between panels, in lines. Default `0.7`.
#' @export
theme_panel_borders <- function(colors = fcs_colors(), spacing = 0.7) {
  theme(panel.border     = element_rect(colour = colors$axis_ticks,
                                        fill = NA, linewidth = 0.4),
        panel.spacing    = unit(spacing, "lines"),
        strip.background = element_rect(fill = colors$grid_major, colour = NA))
}

#' Order cluster labels k1, k2, ... k30 rather than k1, k10, k11, ... k2
#'
#' Cluster ids are "k" plus an integer, and sorting them as text puts k10 before
#' k2. In a 30-cluster legend that is not a cosmetic problem: the reader looks
#' for k7 between k6 and k8 and finds it near the bottom, and colour-to-cluster
#' matching stops being possible at a glance.
#' @param x Character vector of cluster ids.
#' @return `x` as a factor in natural order.
#' @export
natural_cluster_factor <- function(x) {
  u <- unique(as.character(x))
  n <- suppressWarnings(as.numeric(sub("^[A-Za-z]+", "", u)))
  factor(as.character(x), levels = u[order(is.na(n), n, u)])
}

# =============================================================================
