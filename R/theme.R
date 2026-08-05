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
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      plot.title       = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle    = element_text(colour = colors$subtitle_text, size = rel(0.85))
    )
}

# =============================================================================
