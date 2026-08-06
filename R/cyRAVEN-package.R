#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' `cyRAVEN` takes a directory of FCS files and produces a gated, embedded,
#' statistically tested result set. The pipeline runs in nine stages, each of
#' which is also callable on its own:
#'
#' \enumerate{
#'   \item **Read** -- [read_fcs_resolved()] resolves `$PnN`/`$PnS` to marker
#'     symbols and applies the spillover matrix when the file carries one.
#'   \item **Transform** -- [derive_cofactor_pooled()] derives the arcsinh
#'     cofactor from the data rather than assuming one.
#'   \item **Gate** -- [apply_gate_hierarchy()] derives scatter, singlet,
#'     viability and CD45 gates per sample; [density_valley()] places each
#'     marker threshold at that sample's own bimodal split.
#'   \item **Score** -- [score_populations()] evaluates a declarative population
#'     specification (see the `config` vignette).
#'   \item **Embed** -- [run_umap()] builds one shared embedding per marker
#'     panel from a size-balanced subsample.
#'   \item **Test** -- [stats_group_comparison()] for abundance,
#'     [stats_marker_state()] for marker state, both on per-sample values.
#'   \item **Diagnose** -- [batch_mixing_report()], [stats_threshold_drift()],
#'     [stats_confounding()] and [cluster_gate_agreement()].
#'   \item **Explain** -- [explain_cluster()] runs the gating stage backwards:
#'     given cells the specification does not describe, it learns a two-marker
#'     gating strategy that selects them. Descriptive only; it proposes gates
#'     and never alters a scored population.
#'   \item **Report** -- figures and CSVs, plus [write_run_manifest()].
#' }
#'
#' [run_cyraven()] runs all of it. `system.file("scripts", "cyraven.R", package
#' = "cyRAVEN")` is a command-line front end to the same function.
#'
#' @section Statistical stance:
#' Replicates are samples, never cells. Anything computed over pooled cells is
#' labelled descriptive and carries no p-value, because the number of cells is a
#' property of acquisition rather than of the design. See the
#' `vignette("statistics", package = "cyRAVEN")`.
#'
#' @section Options:
#' \describe{
#'   \item{`cyRAVEN.verbose`}{one of `"none"`, `"inform"` (default) or
#'     `"debug"`; controls progress messaging. All messaging goes through
#'     [message()], so [suppressMessages()] also works.}
#'   \item{`cyRAVEN.max_raster_px`}{hard ceiling on any figure's pixel
#'     dimension (default 30000), below the graphics device limit.}
#' }
#'
#' @name cyRAVEN-package
#' @aliases cyRAVEN
#'
#' @import ggplot2
#' @import stats
#' @importFrom grDevices col2rgb convertColor
#' @importFrom grid unit gpar viewport pushViewport popViewport grid.newpage
#'   grid.draw grid.text grid.layout
#' @importFrom utils capture.output head modifyList packageVersion read.csv
#'   sessionInfo str tail write.csv
NULL

# Names that R CMD check cannot resolve statically because they are column names
# evaluated inside ggplot2 aesthetics or data.table expressions. Declaring them
# is the documented remedy; it is not a way of silencing a genuine unbound
# variable, and every name below is a column that exists in the data frame the
# expression is evaluated against.
utils::globalVariables(c(
  ".batch", ".cluster", ".lisi", ".txt_col", "cells_per_ul", "clr", "cond",
  "count", "density", "event_index", "facet", "group", "is_control", "label",
  "lisi", "marker", "md", "mean_val", "median_asinh", "n_cells", "pair",
  "panel", "pct", "pct_of_cd45_pos", "pct_positive", "population",
  "population_label", "qc_status", "sample_id", "sub", "threshold", "umap_1",
  "umap_2", "value", "x", "y",
  # ggplot2 aesthetics evaluated against columns built inside the figure
  # functions themselves, and therefore invisible to static analysis.
  "cd45", "fsc", "lab", "needs_review", "row_label", "ssc", "tag",
  # Uncertainty figures: `lo`/`hi` are the interval bounds and `grp` the study
  # group, all three assigned onto the frame a few lines above the plot call;
  # `u_pct_points` and `term` are columns of the budget table.
  "grp", "hi", "lo", "term", "u_pct_points"
))
