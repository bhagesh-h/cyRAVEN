# =============================================================================
# Package state, options, and messaging
# =============================================================================
#
# WHY THIS FILE EXISTS AT ALL. As a script, this codebase kept two mutable
# globals -- `COLORS` (rewritten by a --config colors: block) and
# `REFERENCE_GROUP` (the cohort every comparison is drawn against) -- and every
# figure function picked them up through a default argument evaluated at call
# time. That is a perfectly reasonable pattern in a script and an illegal one in
# a package: a namespace is locked after loading, so `COLORS <- ...` raises
# "cannot change value of locked binding", and assigning into the global
# environment instead is exactly what a package must never do.
#
# The remedy is the standard one: a private environment owned by the package,
# reached through accessors. Call sites keep the same shape -- a default of
# `colors = fcs_colors()` reads no worse than `colors = COLORS`, and is still
# evaluated at call time, so a re-theme still reaches every figure with no
# argument threading. What changes is that the state is now owned, inspectable,
# settable and restorable rather than ambient.

.cyRAVEN <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .cyRAVEN$colors <- default_colors()
  .cyRAVEN$reference_group <- NULL
  op <- options()
  defaults <- list(
    cyRAVEN.verbose = "inform",
    cyRAVEN.max_raster_px = 30000L
  )
  toset <- !(names(defaults) %in% names(op))
  if (any(toset)) options(defaults[toset])
  invisible()
}

#' The colour palette every figure draws with
#'
#' @description
#' `fcs_colors()` returns the active palette; `set_fcs_colors()` replaces it,
#' merging over the built-in defaults so a partial specification does not blank
#' out the keys it omits.
#'
#' @details
#' Held in a private package environment rather than a global variable, so that
#' re-theming a run cannot leak into another session and the state can be
#' restored (`set_fcs_colors(NULL)`).
#'
#'   `NULL` to restore the defaults.
#' @return `fcs_colors()` returns the active named list of colours.
#'   `set_fcs_colors()` returns the previous value invisibly, so it can be
#'   restored with [on.exit()].
#' @seealso [default_colors()] for the full set of keys and what each controls.
#' @examples
#' old <- set_fcs_colors(list(gate_highlight = "#0072F0"))
#' fcs_colors()$gate_highlight
#' set_fcs_colors(old)
#' @export
fcs_colors <- function() {
  .cyRAVEN$colors %||% default_colors()
}

#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()].
#' @rdname fcs_colors
#' @export
set_fcs_colors <- function(colors) {
  old <- .cyRAVEN$colors
  .cyRAVEN$colors <- if (is.null(colors)) default_colors()
                     else utils::modifyList(default_colors(), colors)
  invisible(old)
}

#' The cohort every comparison is drawn against
#'
#' @description
#' The reference group is figure-wide state: it decides which cohort is plotted
#' first, which one is drawn unfilled, which one keeps the reference colour, and
#' which one every statistical comparison is made against. Setting it once is
#' what keeps those four answers consistent across a result set.
#'
#' @return `reference_group()` returns the active label or `NULL`.
#'   `set_reference_group()` returns the previous value invisibly.
#' @examples
#' old <- set_reference_group("Healthy controls")
#' reference_group()
#' set_reference_group(old)
#' @export
reference_group <- function() {
  .cyRAVEN$reference_group
}

#' @param group The group.
#' @rdname reference_group
#' @export
set_reference_group <- function(group) {
  old <- .cyRAVEN$reference_group
  .cyRAVEN$reference_group <- group
  invisible(old)
}

# --- messaging ---------------------------------------------------------------
#
# WHY message() AND NOT cat(). A package must not write to stdout: cat() cannot
# be suppressed by suppressMessages(), cannot be captured by the calling code,
# and pollutes the output of anything that runs the pipeline programmatically.
# message() goes to stderr, is suppressible, and is what R CMD check expects.
#
# WHY A LEVEL AND NOT A verbose = TRUE ARGUMENT. Threading a flag through forty
# functions makes every signature longer and still cannot be set once for a whole
# run. An option can, and rOpenSci's guidance is explicit that verbosity belongs
# at package level rather than in each call.

fcs_verbosity <- function() {
  v <- getOption("cyRAVEN.verbose", "inform")
  if (!is.character(v) || length(v) != 1L || !v %in% c("none", "inform", "debug"))
    "inform" else v
}

#' @noRd
log_msg <- function(...) {
  if (fcs_verbosity() == "none") return(invisible(NULL))
  message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
}

#' @noRd
log_step <- function(...) {
  if (fcs_verbosity() == "none") return(invisible(NULL))
  bar <- strrep("=", 78)
  message("\n", bar, "\n", ..., "\n", bar)
}

#' @noRd
log_debug <- function(...) {
  if (fcs_verbosity() != "debug") return(invisible(NULL))
  message("[debug] ", ...)
}
