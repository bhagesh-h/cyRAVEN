# SECTION 5 -- DECLARATIVE POPULATION DEFINITIONS
# =============================================================================

#' Built-in population spec, derived from the supplied gating strategy
#'
#' This is DATA, not code: each population is a list of (marker, direction)
#' requirements evaluated against the resolved thresholds, all within the CD45+
#' parent. Users add or edit populations in the YAML config without touching R.
#' A population whose markers are absent from the panel is reported UNAVAILABLE.
#' This is a LITERAL transcription of the 15 populations in the gating strategy
#' document, in document order. Each entry is (marker, direction) with direction
#' one of "above", "below" or "intermediate". Nothing is added, renamed or
#' inferred: populations the document does not list (e.g. dendritic cells,
#' double-negative T cells, NKT cells) are deliberately absent, and gates the
#' document specifies as double-positive are kept double-positive.
#'
#' `lineage` marks the three monocyte subsets, because the document scopes the
#' functional marker blocks by monocyte vs non-monocyte membership.
#' @export
default_population_spec <- function() {
  list(
    "CD14 pos"                        = list(CD14 = "above"),
    "CD16 pos"                        = list(CD16 = "above"),
    "CD14 int CD16 int"               = list(CD14 = "intermediate",
                                            CD16 = "intermediate"),
    "CD3 pos"                         = list(CD3 = "above"),
    "CD3 pos CD4 pos"                 = list(CD3 = "above", CD4 = "above"),
    "CD3 pos CD4 pos CD127 lo CD25 hi"= list(CD3 = "above", CD4 = "above",
                                            CD127 = "below", CD25 = "above"),
    "CD3 pos CD4 pos Vd1 pos"         = list(CD3 = "above", CD4 = "above",
                                            `TCR-Vd1` = "above"),
    "CD3 pos CD4 pos Vd2 pos"         = list(CD3 = "above", CD4 = "above",
                                            `TCR-Vd2` = "above"),
    "CD3 pos CD8 pos"                 = list(CD3 = "above", CD8 = "above"),
    "CD3 pos CD8 pos Vd1 pos"         = list(CD3 = "above", CD8 = "above",
                                            `TCR-Vd1` = "above"),
    "CD3 pos CD8 pos Vd2 pos"         = list(CD3 = "above", CD8 = "above",
                                            `TCR-Vd2` = "above"),
    "CD3 pos Vd1 pos"                 = list(CD3 = "above", `TCR-Vd1` = "above"),
    "CD3 pos Vd2 pos"                 = list(CD3 = "above", `TCR-Vd2` = "above"),
    "CD56 pos NKG2D pos"              = list(CD56 = "above", NKG2D = "above"),
    "CD19 pos"                        = list(CD19 = "above")
  )
}

#' The three monocyte subsets, named exactly as in the population spec.
#' WHY a named helper: the document scopes functional blocks as "every subset
#' EXCEPT the monocytes" and "the monocyte subsets only", so both blocks need the
#' same list and it must stay consistent if the spec is edited.
#' @keywords internal
monocyte_populations <- function() c("CD14 pos", "CD16 pos", "CD14 int CD16 int")

#' Derive the upper bound of an "intermediate" gate
#'
#' WHAT: among the marker-positive cells, find the valley separating the
#' intermediate mode from the bright mode; that valley is the upper bound.
#' WHY: "CD14 int CD16 int" is a genuine three-level gate (negative /
#' intermediate / bright), and monocyte subsetting depends on it. Returns NA when
#' the positive fraction has no second valley -- the population is then reported
#' UNAVAILABLE rather than being silently merged with the bright subset.
#' @param tmat The tmat.
#' @param thr The thr.
#' @param parent The parent.
#' @param spec Population specification mapping population name to marker directions. See [default_population_spec()].
#' @keywords internal
derive_intermediate_bounds <- function(tmat, thr, parent, spec) {
  want <- unique(unlist(lapply(spec, function(d)
    names(d)[vapply(d, function(x) identical(x, "intermediate"), logical(1))])))
  want <- intersect(want, colnames(tmat))
  out <- setNames(rep(NA_real_, length(want)), want)
  for (mk in want) {
    lo <- thr[[mk]]
    if (!is.finite(lo)) next
    x <- tmat[parent & tmat[, mk] > lo, mk]
    v <- density_valley(x, bins = 180L, smooth = 4, min_gap_frac = 0.08)
    if (is.finite(v) && v > lo) out[[mk]] <- v
  }
  as.list(out)
}

#' Functional marker blocks, transcribed from the document's four scoping rules
#'
#' Scope is expressed as a RULE, not a hardcoded population list, so that adding
#' a population to the config automatically scores it under the right blocks:
#'   exclude  = populations to omit (everything else is included)
#'   require  = only populations whose definition requires this marker "above"
#'
#' Document rules, verbatim in intent:
#'   exhaustion (LAG-3/TIM-3/PD-1/CTLA-4) -- "in each cell subset except for
#'     monocytes CD14, CD16 and CD14 int CD16 int"
#'   homing (CCR5/CCR9/CXCR3)            -- "in CD3 pos subsets"
#'   activation (CD69/HLA-DR/CD25/CD57/CD38) -- "in CD3 pos subsets"
#'   monocyte (HLA-DR/BTN2A2/BTN3A1/2/3) -- "in CD14 pos, CD16 pos, and
#'     CD14 int CD16 int subsets"
#' @export
default_functional_blocks <- function() {
  list(
    exhaustion = list(markers = c("LAG-3", "TIM-3", "PD-1", "CTLA-4"),
                      exclude = monocyte_populations()),
    homing     = list(markers = c("CCR5", "CCR9", "CXCR3"),
                      require = "CD3"),
    activation = list(markers = c("CD69", "HLA-DR", "CD25", "CD57", "CD38"),
                      require = "CD3"),
    monocyte   = list(markers = c("HLA-DR", "BTN2A2", "BTN3A1/2/3"),
                      populations = monocyte_populations()),
    # EXTENSION, not from the baseline gating-strategy document: intracellular
    # cytokine / degranulation markers, present in the second panel of the test
    # batch. Scoped to lymphocyte subsets (all non-monocyte populations) because
    # these are effector-function readouts. Markers absent from a panel are
    # skipped, so this block is simply inert on the baseline panel.
    cytokine   = list(markers = c("IFNg", "TNF", "IL-17a", "Granzyme B",
                                  "Perforin", "CD107a"),
                      exclude = monocyte_populations())
  )
}

#' Resolve a functional block's scope to concrete population names
#' Precedence: explicit `populations`, else `require`, else all minus `exclude`.
#' @param block The block.
#' @param spec Population specification mapping population name to marker directions. See [default_population_spec()].
#' @param available The available.
#' @keywords internal
resolve_block_populations <- function(block, spec, available) {
  if (!is.null(block$populations)) return(intersect(block$populations, available))
  if (!is.null(block$require)) {
    hit <- names(spec)[vapply(spec, function(d)
      isTRUE(d[[block$require]] == "above"), logical(1))]
    return(intersect(hit, available))
  }
  setdiff(available, block$exclude %||% character(0))
}

#' Evaluate the population spec against one file's thresholded markers
#'
#' @param tmat matrix of transformed marker values (cells x markers, named cols)
#' @param thr named numeric vector of thresholds
#' @param parent logical mask (the CD45+ parent)
#' @param spec Population specification mapping population name to marker directions. See [default_population_spec()]. Default `default_population_spec()`.
#' @return list(masks, labels, unavailable)
#' @param hi_thr named numeric: upper bound for "intermediate" requirements.
#'   Derived by derive_intermediate_bounds(); a marker requested as
#'   "intermediate" without an upper bound makes its population UNAVAILABLE
#'   rather than silently collapsing to "above".
#' @export
score_populations <- function(tmat, thr, parent, spec = default_population_spec(),
                              hi_thr = NULL) {
  avail <- colnames(tmat)
  masks <- list(); unavailable <- list()
  for (nm in names(spec)) {
    req <- spec[[nm]]
    # Check only the DIRECT marker requirements here. "any_of" is a reserved key
    # holding a nested disjunction, not a marker name, and its members are checked
    # separately below with different semantics: a direct requirement that is
    # missing makes the population unavailable, whereas a disjunction survives as
    # long as one member remains. Including the literal key here marks every
    # population carrying a disjunction as "marker not in panel".
    direct <- setdiff(names(req), "any_of")
    miss <- setdiff(direct, avail)
    miss <- c(miss, direct[!is.na(match(direct, avail)) & !is.finite(thr[direct])])
    miss <- unique(miss[!is.na(miss)])
    if (length(miss)) {
      unavailable[[nm]] <- paste("marker(s) not in panel:", paste(miss, collapse = ", "))
      next
    }
    # "any_of" requirements: a list under the reserved key any_of, itself a
    # (marker -> direction) map, satisfied when ANY of its members is satisfied.
    #
    # WHY the DSL needs this: real gating strategies contain genuine disjunctions
    # -- "NK cells = CD3- CD19- and (CD16+ OR CD56+)" is standard, because the two
    # receptors mark overlapping but non-identical NK subsets and the population
    # is defined as their union. Expressing that as two separate populations
    # would double-count the CD16+CD56+ cells and, since each cell takes a single
    # most-specific label, would split one population across two legend entries.
    any_req <- req[["any_of"]]
    req_simple <- req[setdiff(names(req), "any_of")]
    if (!is.null(any_req)) {
      amiss <- setdiff(names(any_req), avail)
      amiss <- c(amiss, names(any_req)[!is.na(match(names(any_req), avail)) &
                                       !is.finite(thr[names(any_req)])])
      # A disjunction survives partial absence: CD16+ OR CD56+ is still a usable
      # gate when only CD56 is in the panel. Report the loss, do not drop the
      # population -- but if EVERY member is missing there is nothing left.
      any_req <- any_req[!names(any_req) %in% amiss]
      if (!length(any_req)) {
        unavailable[[nm]] <- paste("no any_of marker available:",
                                   paste(unique(amiss), collapse = ", "))
        next
      }
      if (length(amiss))
        log_msg("  NOTE ", nm, ": any_of reduced to ",
                paste(names(any_req), collapse = "/"), " (missing ",
                paste(unique(amiss), collapse = ", "), ")")
    }
    req <- req_simple
    int_mk <- names(req)[vapply(req, function(d) identical(d, "intermediate"), logical(1))]
    no_hi <- int_mk[!vapply(int_mk, function(mk)
      is.finite((hi_thr %||% list())[[mk]] %||% NA_real_), logical(1))]
    if (length(no_hi)) {
      unavailable[[nm]] <- paste("no upper bound derivable for intermediate marker(s):",
                                 paste(no_hi, collapse = ", "))
      next
    }
    m <- parent
    for (mk in names(req)) {
      m <- m & switch(req[[mk]],
        above        = tmat[, mk] > thr[[mk]],
        below        = tmat[, mk] < thr[[mk]],
        intermediate = tmat[, mk] > thr[[mk]] & tmat[, mk] < hi_thr[[mk]],
        stop("unknown direction '", req[[mk]], "' for marker ", mk,
             " (use above/below/intermediate)"))
    }
    if (!is.null(any_req)) {
      any_m <- rep(FALSE, nrow(tmat))
      for (mk in names(any_req))
        any_m <- any_m | switch(any_req[[mk]],
          above = tmat[, mk] > thr[[mk]],
          below = tmat[, mk] < thr[[mk]],
          stop("any_of supports only above/below (marker ", mk, ")"))
      m <- m & any_m
    }
    masks[[nm]] <- m
  }
  # most-specific label per cell: more requirements = deeper definition.
  # An any_of group counts as ONE requirement, not as its member count: a
  # disjunction constrains a cell less than the same number of AND terms, so
  # counting members would rank a loose gate as more specific than a tight one.
  depth <- vapply(spec[names(masks)], function(d)
    length(setdiff(names(d), "any_of")) + (!is.null(d[["any_of"]])), integer(1))
  ord <- names(masks)[order(depth, decreasing = TRUE)]
  labels <- rep(NA_character_, nrow(tmat))
  for (nm in ord) labels[is.na(labels) & masks[[nm]]] <- nm
  labels[is.na(labels) & parent] <- "Other CD45+"
  if (length(unavailable))
    log_msg("  UNAVAILABLE populations: ", paste(names(unavailable), collapse = "; "))
  list(masks = masks, labels = labels, unavailable = unavailable)
}

# =============================================================================
