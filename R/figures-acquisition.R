# SECTION 9e -- ACQUISITION-TIME QC FIGURE
# =============================================================================

#' Event rate across the acquisition, with the flagged intervals marked
#'
#' One panel per sample. The line is how many events the instrument recorded in
#' each equal-width slice of the acquisition, and the marked points are the
#' slices whose rate or whose channel medians departed from the rest of that
#' file.
#'
#' WHAT TO LOOK FOR. A clean acquisition is a roughly flat line: the instrument
#' delivered cells at a steady rate and the signal did not move. A sustained
#' trough is a partial clog. A spike is usually a bubble or a sample that was
#' disturbed. A step is a settings change part-way through the tube.
#'
#' WHY THE RATE AND NOT THE SIGNAL IS PLOTTED. Both are tested, and the flag on
#' each point reflects whichever failed, but the rate is the one quantity that is
#' on a single comparable axis for every panel and every instrument. The channel
#' that drove each flag is named in `acquisition_qc.csv` rather than drawn here,
#' because one panel per sample per channel is not a figure anyone reads.
#'
#' @param bins the `bins` element of [run_acquisition_qc()]
#' @param outfile Path to write the figure to.
#' @param summary the `summary` element, used for the subtitle
#' @param dpi Resolution in dots per inch. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. Default `fcs_colors()`.
#' @export
fig_acquisition_qc <- function(bins, outfile, summary = NULL, dpi = 200,
                               colors = fcs_colors()) {
  if (is.null(bins) || !nrow(bins)) {
    log_msg("[fig] no acquisition-time bins, figure skipped")
    return(invisible(NULL))
  }
  d <- bins[is.finite(bins$n_events), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  d$t_mid <- (d$time_from + d$time_to) / 2

  n_bad <- if (!is.null(summary))
    sum(summary$verdict %in% c("unstable", "unstable, minor"), na.rm = TRUE) else
    length(unique(d$sample_id[d$flagged]))
  n_s <- length(unique(d$sample_id))

  fig <- ggplot(d, aes(x = t_mid, y = n_events)) +
    geom_line(linewidth = 0.35, colour = colors$na_fill %||% "grey60") +
    geom_point(data = d[d$flagged, , drop = FALSE], size = 1.1,
               colour = colors$threshold_review %||% "firebrick") +
    facet_wrap(~ sample_id, scales = "free", ncol = 4) +
    labs(title = "Event rate across each acquisition",
         subtitle = paste0(
           "events per equal-width slice of the Time channel; marked slices ",
           "departed from the rest of that file\n",
           n_bad, " of ", n_s, " sample(s) carry at least one flagged slice. ",
           "Nothing is removed: see pct_delta_if_cleaned for what it would cost"),
         x = "acquisition time (file units)", y = "events in slice") +
    theme_cyto(8, colors = colors)
  safe_ggsave(outfile, plot = fig,
              width = min(14, 3.6 * min(4, n_s)),
              height = max(3, 2.4 * ceiling(n_s / 4)), dpi = dpi,
              limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}
