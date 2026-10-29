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
#'
#' ONE DEPARTURE IN HOW, NOT WHAT. The regulatory T cell population keeps the
#' document's name and markers, but its CD25 and CD127 cuts are derived within
#' the CD3+ CD4+ T cells (`parent:`) rather than across every CD45+ cell. A
#' CD25 cut across all leukocytes falls wherever the CD25+ cells of every
#' lineage put it, and above it sit activated conventional T cells as well as
#' regulatory ones; the consensus gate is CD25+ CD127lo within CD4+ T cells
#' (Santegoets et al. 2015; Liu et al. 2006; Seddiki et al. 2006).
#' @export
default_population_spec <- function() {
  resolve_population_parents(list(
    "CD14 pos"                        = list(CD14 = "above"),
    "CD16 pos"                        = list(CD16 = "above"),
    "CD14 int CD16 int"               = list(CD14 = "intermediate",
                                            CD16 = "intermediate"),
    "CD3 pos"                         = list(CD3 = "above"),
    "CD3 pos CD4 pos"                 = list(CD3 = "above", CD4 = "above"),
    "CD3 pos CD4 pos CD127 lo CD25 hi"= list(parent = "CD3 pos CD4 pos",
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
  ))
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
    names(d)[vapply(d, function(x)
      identical(x, "intermediate") || identical(x, "bright"), logical(1))])))
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
#' @param scoped Cuts derived within parent populations, from
#'   [derive_scoped_cuts()], or NULL. Used for every marker the specification's
#'   `cut_scope` attribute (see [resolve_population_parents()]) places within a
#'   parent; without it those markers fall back to the CD45+ cut.
#' @param quiet Suppress the per-population NOTEs about `any_of` groups losing a
#'   member. The gate-uncertainty bootstrap re-scores every population once per
#'   replicate, and the note is a property of the SPECIFICATION against the
#'   PANEL -- identical on every replicate, because neither changes. Emitted
#'   from there it produced 2,120 of 2,399 log lines on a 20-sample run, 88% of
#'   the log saying one of two things, which buries the NOTEs that do differ
#'   between samples. Callers that score once leave it FALSE and see it once per
#'   sample; callers that score in a loop pass TRUE.
score_populations <- function(tmat, thr, parent, spec = default_population_spec(),
                              hi_thr = NULL, quiet = FALSE, scoped = NULL) {
  avail <- colnames(tmat)
  masks <- list(); unavailable <- list()
  cut_scope <- attr(spec, "cut_scope")
  for (nm in names(spec)) {
    req <- spec[[nm]]
    # Cuts for this population: the CD45+ ones, with any marker its
    # specification places within a parent population replaced by the cut
    # derived on that parent's events. See resolve_population_parents().
    sc <- cut_scope[[nm]]
    thr_pop <- thr; hi_pop <- hi_thr %||% list()
    if (length(sc) && !is.null(scoped)) {
      for (mk in names(sc)[!is.na(sc)]) {
        ws <- scoped[[sc[[mk]]]]
        thr_pop[[mk]] <- if (mk %in% names(ws$thr)) ws$thr[[mk]] else NA_real_
        hi_pop[[mk]] <- ws$hi[[mk]] %||% NA_real_
      }
    }
    # Check only the DIRECT marker requirements here. "any_of" is a reserved key
    # holding a nested disjunction, not a marker name, and its members are checked
    # separately below with different semantics: a direct requirement that is
    # missing makes the population unavailable, whereas a disjunction survives as
    # long as one member remains. Including the literal key here marks every
    # population carrying a disjunction as "marker not in panel".
    direct <- setdiff(names(req), "any_of")
    miss <- setdiff(direct, avail)
    miss <- c(miss, direct[!is.na(match(direct, avail)) & !is.finite(thr_pop[direct])])
    miss <- unique(miss[!is.na(miss)])
    if (length(miss)) {
      within <- miss[miss %in% names(sc)[!is.na(sc)] & miss %in% avail]
      unavailable[[nm]] <- if (length(within) && !is.null(scoped))
        paste0("no cut derivable within ", sc[[within[1]]], " for: ",
               paste(within, collapse = ", "))
      else paste("marker(s) not in panel:", paste(miss, collapse = ", "))
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
                                       !is.finite(thr_pop[names(any_req)])])
      # A disjunction survives partial absence: CD16+ OR CD56+ is still a usable
      # gate when only CD56 is in the panel. Report the loss, do not drop the
      # population -- but if EVERY member is missing there is nothing left.
      any_req <- any_req[!names(any_req) %in% amiss]
      if (!length(any_req)) {
        unavailable[[nm]] <- paste("no any_of marker available:",
                                   paste(unique(amiss), collapse = ", "))
        next
      }
      if (length(amiss) && !isTRUE(quiet))
        log_msg("  NOTE ", nm, ": any_of reduced to ",
                paste(names(any_req), collapse = "/"), " (missing ",
                paste(unique(amiss), collapse = ", "), ")")
    }
    req <- req_simple
    # "bright" needs the same upper bound as "intermediate" -- it is the band
    # ABOVE it. Both are checked together, or a bright requirement with no
    # derivable bound would fall through and be scored as if it were "above",
    # silently returning the whole positive population under a name that claims
    # only its upper mode.
    int_mk <- names(req)[vapply(req, function(d)
      identical(d, "intermediate") || identical(d, "bright"), logical(1))]
    no_hi <- int_mk[!vapply(int_mk, function(mk)
      is.finite(hi_pop[[mk]] %||% NA_real_), logical(1))]
    if (length(no_hi)) {
      unavailable[[nm]] <- paste("no upper bound derivable for intermediate marker(s):",
                                 paste(no_hi, collapse = ", "))
      next
    }
    m <- parent
    for (mk in names(req)) {
      m <- m & switch(req[[mk]],
        above        = tmat[, mk] > thr_pop[[mk]],
        below        = tmat[, mk] < thr_pop[[mk]],
        intermediate = tmat[, mk] > thr_pop[[mk]] & tmat[, mk] < hi_pop[[mk]],
        # The band above the intermediate one: CD56-bright NK cells and the
        # like, where the canonical subset is the upper mode of a bimodal
        # positive population rather than the whole of it.
        bright       = tmat[, mk] > hi_pop[[mk]],
        stop("unknown direction '", req[[mk]], "' for marker ", mk,
             " (use above/below/intermediate/bright)"))
    }
    if (!is.null(any_req)) {
      any_m <- rep(FALSE, nrow(tmat))
      for (mk in names(any_req))
        any_m <- any_m | switch(any_req[[mk]],
          above = tmat[, mk] > thr_pop[[mk]],
          below = tmat[, mk] < thr_pop[[mk]],
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
  labels[is.na(labels) & parent] <- catch_all_label()
  if (length(unavailable))
    log_msg("  UNAVAILABLE populations: ", paste(names(unavailable), collapse = "; "))
  list(masks = masks, labels = labels, unavailable = unavailable)
}

# =============================================================================

#' Resolve `parent:` declarations in a population specification
#'
#' WHY. A cut on the whole CD45+ parent sits wherever the lineages that
#' dominate that mixture put a valley, which is not where a subset is defined.
#' Regulatory T cells are CD25+ CD127lo WITHIN CD4+ T cells (Santegoets et al.
#' 2015), and the CD25 valley across all leukocytes is not the one that
#' separates them from activated conventional CD4 T cells (Seddiki et al.
#' 2006); monocyte subsets are split by CD14 and CD16 within monocytes (Abeles
#' et al. 2012). Gating templates express this as a hierarchy, each gate
#' derived on its parent's events (openCyto; Finak et al. 2014).
#'
#' WHAT. A population may name another in its `parent:` key:
#'
#' ```yaml
#' Tregs:
#'   parent: CD4 T cells
#'   CD25: above
#'   CD127: below
#' ```
#'
#' Its own markers' cuts are then derived on the parent's events (see
#' [derive_scoped_cuts()]), and it is scored within the parent. The parent's
#' conditions are merged into the child here, so the definition every other
#' part of the package reads is the complete one; attribute `cut_scope` records,
#' per population and marker, which parent's events that marker's cut comes
#' from (NA for the CD45+ cut), and `parents` the declared parent of each
#' population.
#'
#' A child may not re-declare a marker its parent requires in a different
#' direction, and parent and child cannot both carry an `any_of` group: both
#' would leave the child's relation to its parent undefined.
#'
#' @param spec population specification as read from YAML.
#' @return the specification with every `parent:` resolved, carrying
#'   attributes `cut_scope` and `parents`; unchanged when no population names a
#'   parent.
#' @export
resolve_population_parents <- function(spec) {
  if (!length(spec)) return(spec)
  has_par <- vapply(spec, function(d) is.list(d) && !is.null(d[["parent"]]), logical(1))
  if (!any(has_par)) return(spec)
  resolved <- list(); scope <- list(); visiting <- character(0)
  marker_keys <- function(d) unique(c(setdiff(names(d), "any_of"), names(d[["any_of"]])))
  res <- function(nm) {
    if (!is.null(resolved[[nm]])) return(resolved[[nm]])
    if (nm %in% visiting)
      stop("population specification: the parent chain through '", nm,
           "' loops back on itself", call. = FALSE)
    d <- spec[[nm]]
    p <- d[["parent"]]; d[["parent"]] <- NULL
    own <- marker_keys(d)
    if (is.null(p)) {
      resolved[[nm]] <<- d
      scope[[nm]] <<- stats::setNames(rep(NA_character_, length(own)), own)
      return(d)
    }
    p <- as.character(p)[1]
    if (!p %in% names(spec))
      stop("population '", nm, "' names parent '", p, "', which the ",
           "specification does not declare", call. = FALSE)
    visiting <<- c(visiting, nm)
    pd <- res(p)
    visiting <<- setdiff(visiting, nm)
    merged <- pd[setdiff(names(pd), "any_of")]
    sc <- scope[[p]]
    for (k in setdiff(names(d), "any_of")) {
      if (k %in% names(merged)) {
        if (!identical(as.character(merged[[k]]), as.character(d[[k]])))
          stop("population '", nm, "' requires ", k, " ", d[[k]], " but its ",
               "parent '", p, "' requires ", k, " ", merged[[k]], call. = FALSE)
        next                     # restated: still the parent's condition and cut
      }
      merged[[k]] <- d[[k]]
      sc[[k]] <- p
    }
    if (!is.null(d[["any_of"]]) && !is.null(pd[["any_of"]]))
      stop("population '", nm, "' and its parent '", p, "' both carry an any_of ",
           "group; a population holds one disjunction", call. = FALSE)
    if (!is.null(d[["any_of"]])) {
      merged[["any_of"]] <- d[["any_of"]]
      for (k in names(d[["any_of"]])) sc[[k]] <- p
    } else if (!is.null(pd[["any_of"]])) merged[["any_of"]] <- pd[["any_of"]]
    resolved[[nm]] <<- merged; scope[[nm]] <<- sc
    merged
  }
  for (nm in names(spec)) res(nm)
  out <- resolved[names(spec)]
  attr(out, "cut_scope") <- scope[names(spec)]
  attr(out, "parents") <- vapply(spec, function(d)
    as.character(d[["parent"]] %||% NA_character_)[1], "")
  out
}

#' Derive the cuts a specification places within parent populations
#'
#' For each parent named in the specification's `cut_scope` (see
#' [resolve_population_parents()]), the parent is scored in this sample and
#' every marker a child places within it is cut on the parent's events alone,
#' with the same rule the CD45+ cuts use (`resolve`). Bounds for `bright` and
#' `intermediate` requirements are derived within the parent too. Parents are
#' processed outermost first, so a parent that is itself scoped is scored with
#' its own scoped cuts.
#'
#' A parent with fewer than `min_events` events gives no cut: its children are
#' reported unavailable in that sample rather than cut on a handful of cells.
#'
#' @param tmat transformed marker matrix.
#' @param thr named CD45+ cuts.
#' @param parent the CD45+ mask.
#' @param spec resolved specification.
#' @param hi_thr CD45+ upper bounds, from [derive_intermediate_bounds()].
#' @param resolve `function(marker, x)` returning a [resolve_threshold()]
#'   result for values `x`. Default: `resolve_threshold()` with its defaults.
#' @param min_events fewest parent events a cut is derived from. Default 200,
#'   the least the tail rule accepts.
#' @return list keyed by parent population, each `list(thr, hi, details,
#'   n_parent)`, or NULL when the specification scopes nothing.
#' @export
derive_scoped_cuts <- function(tmat, thr, parent, spec, hi_thr = NULL,
                               resolve = NULL, min_events = 200L) {
  cs <- attr(spec, "cut_scope")
  scopes <- unique(stats::na.omit(unlist(cs, use.names = FALSE)))
  if (!length(scopes)) return(NULL)
  resolve <- resolve %||% function(m, x) resolve_threshold(m, x)
  pars <- attr(spec, "parents")
  depth <- function(p) {
    d <- 0L
    while (!is.null(p) && p %in% names(pars) && !is.na(pars[[p]])) {
      p <- pars[[p]]; d <- d + 1L
    }
    d
  }
  scopes <- scopes[order(vapply(scopes, depth, integer(1)))]
  out <- list()
  for (P in scopes) {
    sp <- score_populations(tmat, thr, parent, spec, hi_thr = hi_thr, quiet = TRUE,
                            scoped = out)
    m <- sp$masks[[P]]
    n <- if (is.null(m)) 0L else sum(m)
    kids <- names(cs)[vapply(cs, function(x) any(x %in% P), logical(1))]
    mk <- intersect(unique(unlist(lapply(kids, function(k) names(cs[[k]])[cs[[k]] %in% P]))),
                    colnames(tmat))
    thrP <- stats::setNames(rep(NA_real_, length(mk)), mk); det <- list()
    for (k in mk) {
      if (n < min_events) {
        det[[k]] <- list(threshold = NA_real_, needs_review = TRUE,
                         source = sprintf("none: %s holds %d event(s), fewer than %d",
                                          P, n, as.integer(min_events)))
        next
      }
      rr <- resolve(k, tmat[m, k])
      thrP[[k]] <- rr$threshold; det[[k]] <- rr
    }
    # bright and intermediate bounds, within the parent, for the markers placed there
    sub <- lapply(kids, function(k) spec[[k]][intersect(names(spec[[k]]),
                                                       names(cs[[k]])[cs[[k]] %in% P])])
    hiP <- if (n >= min_events) derive_intermediate_bounds(tmat, thrP, m, sub) else list()
    out[[P]] <- list(thr = thrP, hi = hiP, details = det, n_parent = n)
  }
  out
}

# =============================================================================
# THE CATCH-ALL LABEL
# =============================================================================
#
# WHY THIS IS ONE DEFINITION AND NOT A STRING REPEATED AROUND THE PACKAGE. The
# label given to cells inside the parent gate that match no declared population
# is not a population -- it is the remainder, and every consumer that reasons
# about coverage has to be able to tell it apart from a real lineage. That test
# was written out by hand in several places and, in one of them, the list of
# names to treat as "not a real population" simply did not contain the name this
# package actually assigns. The consequence was silent and total: explore mode's
# findings table, whose entire purpose is to report the clusters the
# specification fails to describe, scored every cluster as 0% unlabelled and
# declared all of them covered -- including clusters that were 100% remainder --
# on a run whose specification described under a tenth of the cells.

#' The label assigned to cells inside the parent gate matching no population
#' @return A single string.
#' @keywords internal
catch_all_label <- function() "Other CD45+"

#' Is this population label the catch-all remainder rather than a real lineage?
#'
#' Matches the label this package assigns plus the conventional spellings an
#' imported or hand-edited table may carry, so a caller reasoning about coverage
#' cannot miss the remainder on a technicality of spelling.
#' @param x Character vector of population labels.
#' @return Logical vector, `TRUE` where the label is a catch-all.
#' @keywords internal
is_catch_all_label <- function(x) {
  x <- as.character(x)
  is.na(x) | !nzchar(trimws(x)) |
    grepl("^[[:space:]]*(other|unclassified|unlabelled|unlabeled|unassigned|none|ungated)([[:space:]]|$|[+-])",
          x, ignore.case = TRUE)
}
