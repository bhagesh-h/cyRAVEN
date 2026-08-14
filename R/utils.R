# =============================================================================
# Small shared helpers
# =============================================================================

#' Default value for `NULL` or empty
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Case- and whitespace-insensitive identifier matching
#'
#' Used everywhere a hand-typed sample map is joined against a clinically
#' exported patient table, where `"HS5_004"` and `"hs5_004"` are the same donor
#' and a plain `merge()` would silently produce all-NA covariates.
#' @noRd
#' @param x A vector of values.
norm_id <- function(x) toupper(trimws(as.character(x)))

#' Save a ggplot at a DPI that cannot exceed the device's raster ceiling
#'
#' @description
#' A drop-in for [ggplot2::ggsave()] that lowers `dpi` -- never raises it -- when
#' the requested canvas would produce a raster larger than the graphics device
#' can allocate (about 32,767 px per side).
#'
#' @details
#' Several figures scale their width or height with the number of samples,
#' populations or markers, so a large cohort or a wide panel can ask for a canvas
#' bigger than the device supports; the run would then die on a figure rather
#' than on the analysis. `ggsave()`'s own `limitsize` guards the opposite mistake
#' -- an accidentally huge size in inches -- and does not help here.
#'
#' Reducing DPI keeps the requested *physical* size, and therefore every relative
#' font, point and line size in the theme, exactly as specified. At extreme scale
#' the figure becomes less sharp; it does not become smaller or more crowded, and
#' the run finishes.
#'
#' The ceiling is `getOption("cyRAVEN.max_raster_px")`, default 30000.
#'
#' @param filename Output path.
#' @param plot A ggplot or patchwork object.
#' @param width,height Canvas size in inches.
#' @param dpi Requested resolution; the effective value may be lower.
#' @param ... Passed to [ggplot2::ggsave()].
#' @return The path, invisibly.
#' @keywords internal
safe_ggsave <- function(filename, plot, width, height, dpi = 300, ...) {
  max_px <- getOption("cyRAVEN.max_raster_px", 30000L)
  dpi_cap <- floor(max_px / max(width, height, 0.1))
  dpi_use <- max(36L, min(dpi, dpi_cap))
  if (dpi_use < dpi)
    log_msg("  NOTE ", basename(filename), ": ", round(width, 1), "x",
            round(height, 1), "in at ", dpi, " dpi would exceed the ",
            max_px, "px raster limit; writing at ", dpi_use,
            " dpi instead (physical size, and everything's relative size,",
            " is unchanged)")
  # WHITE, NOT TRANSPARENT. ggsave() takes its background from the theme, and a
  # theme that leaves plot.background unset writes an RGBA PNG with a
  # transparent ground. That looks fine in a white report and is unreadable
  # anywhere else: click a figure to expand it in a dark viewer, or drop it on a
  # dark slide, and black axis text sits on whatever is behind it. Measured on
  # one run: 11 of 20 figures were written RGBA.
  #
  # Set here rather than only in theme_cyto() because this is the single point
  # every figure passes through, including the few that build their own theme.
  # An explicit bg from the caller still wins.
  dots <- list(...)
  if (!"bg" %in% names(dots)) dots$bg <- "white"
  do.call(ggplot2::ggsave,
          c(list(filename = filename, plot = plot, width = width,
                 height = height, dpi = dpi_use), dots))
  invisible(filename)
}

#' Rows a group comparison may use: real, QC-passing samples only
#'
#' @description
#' Drops unstained controls and staining-QC failures from a table before it
#' reaches a figure or a test.
#'
#' @details
#' Both checks are needed, and `is_control` alone is not enough. Population
#' scoring runs for every gated sample regardless of its staining verdict --
#' the frequency, MFI and functional tables deliberately keep failed and control
#' rows so the CSV exports stay auditable. A sample that FAILED staining (no
#' separable CD45+ mode, or below the minimum CD45 percentage) has
#' `is_control = FALSE` and `qc_status = "failed"`; its percentages carry no more
#' evidence than a control's, and without the second check they would enter every
#' abundance figure and test by default rather than only under
#' `--include-qc-failed`. That flag works by relabelling a forced sample's
#' `qc_status` to `"pass"`, so filtering on `qc_status` here is what gives the
#' flag any meaning for figures and tests rather than for the embedding alone.
#'
#' @param d A table carrying `is_control` and/or `qc_status`, whichever are
#'   present. Tables that carry neither (for example derived population ratios,
#'   which are already built only from passing rows) are returned unchanged.
#' @return `d` with disqualified rows removed.
#' @keywords internal
qc_pass_rows <- function(d) {
  if ("is_control" %in% names(d)) d <- d[!d$is_control, , drop = FALSE]
  if ("qc_status" %in% names(d)) d <- d[d$qc_status == "pass", , drop = FALSE]
  d
}
