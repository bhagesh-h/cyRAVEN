# SECTION 8 -- CONFIG EMISSION
# =============================================================================

#' Write the derived values out as an editable YAML config
#'
#' Every threshold carries its source and a review flag. The config is an
#' OVERRIDE mechanism: the script derives everything at runtime and runs without
#' it. Values under `derived_from_batch` document what THIS batch produced --
#' editing them changes future runs; deleting the file changes nothing.
#' @param path File path.
#' @param derived The derived.
#' @param cofactors The cofactors.
#' @param spec Population specification mapping population name to marker directions. See [default_population_spec()]. Default `default_population_spec()`.
#' @param blocks Functional-marker blocks. See [default_functional_blocks()]. Default `default_functional_blocks()`.
#' @param colmap The colmap. Default `default_column_map()`.
#' @param valmap The valmap. Default `default_value_map()`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `default_colors()`.
#' @keywords internal
write_config <- function(path, derived, cofactors, spec = default_population_spec(),
                         blocks = default_functional_blocks(),
                         colmap = default_column_map(), valmap = default_value_map(),
                         colors = default_colors()) {
  cfg <- list(
    meta = list(
      written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      note = paste("Auto-generated from a derivation run. All values are OPTIONAL",
                   "overrides; the script derives them at runtime when absent.",
                   "Thresholds are on the asinh scale using the cofactor below.")),
    transform = list(method = "asinh",
                     cofactor = list(mode = "auto",
                                     derived_from_batch = as.list(cofactors))),
    gating = list(
      scatter_gate = list(mode = "auto",
                          note = "lower FSC bound = deepest log10 FSC-A density valley"),
      singlet_band = list(mode = "auto", mad_k = 3,
                          note = "median(FSC-H/FSC-A) +/- k*MAD within the scatter gate"),
      staining_qc = list(min_cd45_pct_of_live = 5)),
    thresholds = derived,
    populations = spec,
    functional_blocks = blocks,
    metadata = list(column_map = colmap, value_translations = valmap),
    colors = colors)
  yaml::write_yaml(cfg, path)
  log_msg("wrote config: ", path)
  invisible(cfg)
}

# =============================================================================
# SECTION 9 -- SESSION STATE PERSISTENCE
# =============================================================================

#' Save the complete analysis state to an .RData file
#'
#' WHY: reading and gating a multi-hundred-megabyte batch is the expensive part
#' of this pipeline. Saving the state means re-plotting, re-thresholding or
#' re-embedding never requires re-reading the FCS files.
#'
#' The raw expression matrices are dropped by default (they are large and
#' re-readable from the FCS files); everything derived -- masks, thresholds,
#' coordinates, tables, gate geometry -- is kept. Pass keep_exprs = TRUE to
#' include them for a fully self-contained but much larger file.
#' @param path File path.
#' @param state The state.
#' @param keep_exprs The keep exprs. Default `FALSE`.
#' @keywords internal
save_session <- function(path, state, keep_exprs = FALSE) {
  st <- state
  if (!keep_exprs && !is.null(st$reads))
    st$reads <- lapply(st$reads, function(r) { r$exprs <- NULL; r })
  st$saved_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  st$session_info <- utils::capture.output(utils::sessionInfo())
  save(st, file = path, compress = "xz")
  log_msg("saved session state: ", path, " (",
          round(file.size(path) / 1e6, 1), " MB",
          if (keep_exprs) ", includes raw expression matrices" else
            ", derived state only \u2014 FCS re-read needed for raw events", ")")
  invisible(path)
}

#' Restore a saved session
#' Usage:  st <- load_session("results/session_state.RData")
#'         fig_umap_overview(st$embeddings$panel_1$cells, "new_plot.png")
#' @param path File path.
#' @keywords internal
load_session <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  st <- get("st", envir = e)
  log_msg("restored session from ", path, " (saved ", st$saved_at %||% "?", ")")
  log_msg("  contents: ", paste(names(st), collapse = ", "))
  st
}

# =============================================================================
