# =============================================================================
# IMMUNE SUBSETS: the layer below the lineages
# =============================================================================
#
# THE GAP THIS CLOSES. A specification that names T cells, B cells, NK cells and
# monocytes describes the lineages and stops. Every question a clinical
# immunophenotyping study actually asks lives one level below that: not "how
# many T cells" but "how many are activated, exhausted or senescent"; not "how
# many monocytes" but "how many have lost HLA-DR". A run whose finest label is
# "CD4 T cells" cannot answer any of them, and the unsupervised clusters have
# nothing finer to be matched against, so cluster identity comes back coarse.
#
# WHAT THIS IS NOT. It is not a fixed taxonomy. The subsets that can be defined
# are a property of the PANEL, not of immunology, so this proposes only the ones
# the supplied markers actually resolve and says nothing about the rest. A
# definition that needs a marker the panel does not carry is omitted rather than
# approximated, because an approximated subset carries a real subset's name.
#
# WHAT A TYPICAL PANEL CANNOT DO, and why it is worth stating rather than
# quietly skipping:
#
#   - Naive / central memory / effector memory / TEMRA need CD45RA together with
#     CCR7 (or CD62L). Without both, any "memory" label is guesswork. cyCONDOR's
#     own reference annotation splits its CD4 and CD8 compartments precisely on
#     CD45RA and CD127, which is why its taxonomy is finer than what a panel
#     lacking CD45RA can reach.
#   - Classical / intermediate / non-classical monocytes are defined by CD14
#     against CD16. Without CD16 the monocyte compartment can be split on
#     activation state but not into its three canonical subsets.
#   - Senescence is conventionally CD28- CD57+. With CD57 alone the label is
#     named "CD57+" rather than "senescent", because CD57 without CD28 marks a
#     late-differentiated population that is not the same set.
#
# SOURCES FOR THE DEFINITIONS USED BELOW. Monocyte HLA-DR as the immunoparalysis
# readout, and the CD14/HLA-DR combination as the way it is measured, follow the
# mHLA-DR literature (Front Immunol 2023;14:1130214). CD38+HLA-DR+ as the
# activated-and-exhausted T cell phenotype, and its association with severity,
# follows the CD38+HLA-DR+ T cell work (PMC12074963). CD38-high monocytes as an
# early sepsis population follow PMC12199442. CD56 bright against dim as the two
# NK compartments is standard. The rest -- PD-1, TIM-3, LAG-3, CTLA-4 for
# exhaustion; CD69 for early activation; CXCR3, CCR5, CCR9 for homing -- are the
# conventional single-marker readouts and are named for the marker, not for an
# inferred function.

#' Immune subset definitions resolvable from a given marker panel
#'
#' Proposes the subsets the supplied markers can actually define, in the same
#' format as [default_population_spec()]: a named list of populations, each a
#' list of marker -> direction requirements.
#'
#' @param markers Character vector of marker names present in the panel.
#' @param parent Marker gating the parent population, added to every definition
#'   so subsets sit inside the same parent as the lineages. `NULL` to omit.
#' @return Named list of population definitions; empty when the panel resolves
#'   none.
#' @examples
#' immune_subset_spec(c("CD3", "CD4", "CD8", "CD38", "HLA-DR"))
#' @export
immune_subset_spec <- function(markers, parent = "CD45") {
  markers <- as.character(markers)
  has <- function(...) all(c(...) %in% markers)
  out <- list()
  add <- function(name, def) {
    # Every requirement must be satisfiable by the panel, or the subset is not
    # proposed at all. Silently dropping one requirement would leave a broader
    # population wearing the narrower one's name.
    if (!all(names(def) %in% markers)) return(invisible(NULL))
    if (!is.null(parent) && parent %in% markers)
      def <- c(stats::setNames(list("above"), parent), def)
    out[[name]] <<- def
    invisible(NULL)
  }

  # ---- the base lineages ----------------------------------------------------
  # WHY THE LINEAGES ARE CANDIDATES TOO. A cluster of plain CD4 T cells -- CD3
  # and CD4 positive, and nothing else over the cut -- matches no QUALIFIED
  # subset, and would come back unnamed while being one of the easiest clusters
  # in the run to name. Including the lineages means every cluster gets the
  # finest name its profile supports: the qualified subset where one fits,
  # the lineage where none does, and NA only when the profile matches neither.
  # Most-specific-wins does the ranking, so a lineage never outranks a subset
  # that also holds.
  if (has("CD3", "CD4", "CD8")) {
    add("CD4 T cells", list(CD3 = "above", CD4 = "above", CD8 = "below"))
    add("CD8 T cells", list(CD3 = "above", CD8 = "above", CD4 = "below"))
  }
  if (has("CD3")) add("T cells", list(CD3 = "above"))
  if (has("CD3", "CD19")) add("B cells", list(CD3 = "below", CD19 = "above"))
  if (has("CD3", "CD56")) add("NK cells", list(CD3 = "below", CD56 = "above"))
  if (has("CD3", "CD56")) add("NKT cells", list(CD3 = "above", CD56 = "above"))
  if (has("CD3", "CD14")) add("Monocytes", list(CD3 = "below", CD14 = "above"))
  if (has("CD3", "CD4", "CD25", "CD127"))
    add("Regulatory T cells", list(CD3 = "above", CD4 = "above",
                                  CD25 = "above", CD127 = "below"))
  for (.vv in c("TCR-Vd1", "TCR-Vd2")) if (has("CD3", .vv))
    add(paste0(sub("^TCR-", "", .vv), " T cells"),
        c(list(CD3 = "above"), stats::setNames(list("above"), .vv)))

  # ---- CD4 T cell subsets ---------------------------------------------------
  if (has("CD3", "CD4", "CD8")) {
    base4 <- list(CD3 = "above", CD4 = "above", CD8 = "below")
    add("CD4 T cells CD38+ HLA-DR+ (activated)",
        c(base4, list(CD38 = "above", `HLA-DR` = "above")))
    add("CD4 T cells CD69+ (early activated)", c(base4, list(CD69 = "above")))
    add("CD4 T cells PD-1+ (exhausted)",       c(base4, list(`PD-1` = "above")))
    add("CD4 T cells TIM-3+",                  c(base4, list(`TIM-3` = "above")))
    add("CD4 T cells LAG-3+",                  c(base4, list(`LAG-3` = "above")))
    add("CD4 T cells CTLA-4+",                 c(base4, list(`CTLA-4` = "above")))
    add("CD4 T cells CXCR3+ (Th1-like)",       c(base4, list(CXCR3 = "above")))
    add("CD4 T cells CCR5+",                   c(base4, list(CCR5 = "above")))
    add("CD4 T cells CCR9+ (gut-homing)",      c(base4, list(CCR9 = "above")))
    add("CD4 T cells CD57+",                   c(base4, list(CD57 = "above")))
  }

  # ---- CD8 T cell subsets ---------------------------------------------------
  if (has("CD3", "CD8", "CD4")) {
    base8 <- list(CD3 = "above", CD8 = "above", CD4 = "below")
    add("CD8 T cells CD38+ HLA-DR+ (activated)",
        c(base8, list(CD38 = "above", `HLA-DR` = "above")))
    add("CD8 T cells CD69+ (early activated)", c(base8, list(CD69 = "above")))
    add("CD8 T cells PD-1+ (exhausted)",       c(base8, list(`PD-1` = "above")))
    add("CD8 T cells TIM-3+",                  c(base8, list(`TIM-3` = "above")))
    add("CD8 T cells LAG-3+",                  c(base8, list(`LAG-3` = "above")))
    # CD57 WITHOUT CD28 IS NOT "SENESCENT". Named for the marker.
    add("CD8 T cells CD57+ (late-differentiated)",
        c(base8, list(CD57 = "above")))
    add("CD8 T cells CXCR3+",                  c(base8, list(CXCR3 = "above")))
    add("CD8 T cells CCR5+",                   c(base8, list(CCR5 = "above")))
    add("CD8 T cells CCR9+ (gut-homing)",      c(base8, list(CCR9 = "above")))
  }

  # ---- regulatory T cell subsets --------------------------------------------
  if (has("CD3", "CD4", "CD25", "CD127")) {
    treg <- list(CD3 = "above", CD4 = "above", CD25 = "above",
                 CD127 = "below")
    add("Regulatory T cells HLA-DR+ (activated)",
        c(treg, list(`HLA-DR` = "above")))
    add("Regulatory T cells CTLA-4+", c(treg, list(`CTLA-4` = "above")))
  }

  # ---- gamma-delta subsets --------------------------------------------------
  for (.v in c("TCR-Vd1", "TCR-Vd2")) {
    if (!has("CD3", .v)) next
    nm <- sub("^TCR-", "", .v)
    g <- c(list(CD3 = "above"), stats::setNames(list("above"), .v))
    add(paste0(nm, " T cells CD69+ (activated)"), c(g, list(CD69 = "above")))
    add(paste0(nm, " T cells PD-1+ (exhausted)"), c(g, list(`PD-1` = "above")))
    add(paste0(nm, " T cells CD57+"),             c(g, list(CD57 = "above")))
    add(paste0(nm, " T cells NKG2D+"),            c(g, list(NKG2D = "above")))
  }

  # ---- NK subsets -----------------------------------------------------------
  if (has("CD3", "CD56")) {
    nk <- list(CD3 = "below", CD56 = "above")
    # BRIGHT AGAINST DIM. The dim subset is the intermediate band between the
    # positivity threshold and the valley separating the two modes; the bright
    # subset is above that valley. Both become UNAVAILABLE rather than wrong
    # when no second valley is found -- see derive_intermediate_bounds(). Note
    # that "above" would NOT give bright: it is every CD56-positive cell, dim
    # included, which is simply the NK gate under a narrower name.
    add("NK cells CD56 bright", list(CD3 = "below", CD56 = "bright"))
    add("NK cells CD56 dim",    list(CD3 = "below", CD56 = "intermediate"))
    add("NK cells CD57+ (mature)", c(nk, list(CD57 = "above")))
    add("NK cells NKG2D+",         c(nk, list(NKG2D = "above")))
    add("NK cells HLA-DR+ (activated)", c(nk, list(`HLA-DR` = "above")))
  }

  # ---- B cell subsets -------------------------------------------------------
  if (has("CD3", "CD19")) {
    b <- list(CD3 = "below", CD19 = "above")
    # CD38+ without CD27 is plasmablast-LIKE, not plasmablast. Named for what
    # was measured.
    add("B cells CD38+", c(b, list(CD38 = "above")))
    add("B cells CD25+", c(b, list(CD25 = "above")))
    add("B cells HLA-DR+", c(b, list(`HLA-DR` = "above")))
  }

  # ---- monocyte subsets -----------------------------------------------------
  if (has("CD3", "CD14")) {
    m <- list(CD3 = "below", CD14 = "above")
    # The immunoparalysis readout. Already in most hand-written specs; included
    # here so a panel that has the markers gets it without being told.
    add("Monocytes HLA-DR low (immunoparalysis)",
        c(m, list(`HLA-DR` = "below")))
    add("Monocytes CD38+", c(m, list(CD38 = "above")))
    add("Monocytes CD69+", c(m, list(CD69 = "above")))
    add("Monocytes BTN3A1/2/3+", c(m, list(`BTN3A1/2/3` = "above")))
    add("Monocytes BTN2A2+",     c(m, list(BTN2A2 = "above")))
    # The canonical three monocyte subsets, available only with CD16.
    if (has("CD16")) {
      add("Classical monocytes (CD14+ CD16-)",
          list(CD3 = "below", CD14 = "above", CD16 = "below"))
      add("Intermediate monocytes (CD14+ CD16+)",
          list(CD3 = "below", CD14 = "intermediate", CD16 = "intermediate"))
      add("Non-classical monocytes (CD14low CD16+)",
          list(CD3 = "below", CD14 = "below", CD16 = "above"))
    }
  }

  out
}

#' Report which canonical subsets the panel cannot reach, and why
#'
#' Written to the log so a reader knows the absence of a subset is a property of
#' the panel rather than of the cohort. A missing label that is never mentioned
#' reads as a population that was looked for and not found.
#'
#' @param markers Character vector of marker names present in the panel.
#' @return Character vector of one sentence per unreachable group.
#' @keywords internal
unreachable_subsets <- function(markers) {
  markers <- as.character(markers)
  msg <- character(0)
  need <- function(what, req) {
    miss <- setdiff(req, markers)
    if (length(miss))
      msg <<- c(msg, sprintf("%s: needs %s, missing %s", what,
                             paste(req, collapse = " + "),
                             paste(miss, collapse = ", ")))
  }
  need("naive / central memory / effector memory / TEMRA T cells",
       c("CD45RA", "CCR7"))
  need("classical / intermediate / non-classical monocytes", c("CD14", "CD16"))
  need("senescent T cells (CD28- CD57+)", c("CD28", "CD57"))
  need("plasmablasts (CD19+ CD27++ CD38++)", c("CD19", "CD27", "CD38"))
  need("memory / naive B cells", c("CD19", "CD27", "IgD"))
  msg
}

# =============================================================================
# ANNOTATING EXPLORE CLUSTERS WITH SUBSETS
# =============================================================================
#
# WHY THE CLUSTER PROFILE AND NOT A GATE. The subsets must not be added to the
# declared specification: doing that re-scores every cell, changes every
# frequency, and asks a circular question of a subset named after the marker
# that defines it. A cluster, though, already carries its own marker profile --
# the fraction of its cells positive for each channel, in
# explore_cluster_profile.csv -- and that is exactly what an annotator needs.
# This is how interactive cluster-annotation tools work: CyCadas
# (Bioinformatics 2024;40:btae595) builds cell types by selecting positive and
# negative markers and reading off which clusters satisfy them, and cyCONDOR's
# metaclustering() is a hand-written version of the same mapping.
#
# So the declared analysis is untouched, and the label attaches to the cluster
# rather than to the cell.
#
# THE MATCH IS ALL-OR-NOTHING PER DEFINITION, and the most specific satisfied
# definition wins, exactly as score_populations() ranks its own labels: a
# definition with more requirements describes a narrower population, so where
# two both hold, the narrower one is the better name. Reporting the margin
# alongside is what stops that being read as certainty -- a cluster that is 51%
# CD69-positive satisfies "CD69 above" on the same footing as one that is 99%,
# and only the margin distinguishes them.

#' Annotate clusters with the immune subset their marker profile matches
#'
#' @param prof Cluster profile: one row per cluster, with a `cluster` column and
#'   one `frac_pos.<marker>` column per marker, as written to
#'   explore_cluster_profile.csv.
#' @param subsets Named list of subset definitions; defaults to the ones the
#'   profile's own markers can resolve.
#' @param pos_cut Fraction of a cluster's cells that must be positive for a
#'   marker before the cluster counts as positive for it.
#' @return data.frame with one row per cluster: `cluster`, `subset_label`,
#'   `subset_markers` (the requirements that were checked), `n_requirements`,
#'   `margin` (the smallest distance from `pos_cut` across those requirements)
#'   and `subset_alternatives`.
#' @keywords internal
annotate_clusters_with_subsets <- function(prof, subsets = NULL, pos_cut = 0.5) {
  if (is.null(prof) || !nrow(prof) || !"cluster" %in% names(prof))
    return(NULL)
  fp <- grep("^frac_pos[.]", names(prof), value = TRUE)
  if (!length(fp)) return(NULL)
  mk <- sub("^frac_pos[.]", "", fp)
  if (is.null(subsets)) subsets <- immune_subset_spec(mk, parent = NULL)
  # An "intermediate" or "bright" requirement needs the upper bound derived from
  # the intensity distribution, which a positivity fraction does not carry.
  # Dropped rather than approximated: treating bright as merely positive would
  # label a whole compartment after its upper mode.
  subsets <- Filter(function(d) all(vapply(d, function(x)
    x %in% c("above", "below"), logical(1))), subsets)
  if (!length(subsets)) return(NULL)

  rows <- lapply(seq_len(nrow(prof)), function(i) {
    v <- suppressWarnings(as.numeric(unlist(prof[i, fp])))
    names(v) <- mk
    hit <- character(0); nreq <- integer(0); marg <- numeric(0)
    for (nm in names(subsets)) {
      d <- subsets[[nm]]
      if (!all(names(d) %in% mk)) next
      val <- v[names(d)]
      if (anyNA(val)) next
      ok <- vapply(seq_along(d), function(j)
        if (identical(d[[j]], "above")) val[j] >= pos_cut else val[j] < pos_cut,
        logical(1))
      if (!all(ok)) next
      hit  <- c(hit, nm)
      nreq <- c(nreq, length(d))
      marg <- c(marg, min(abs(val - pos_cut)))
    }
    if (!length(hit))
      return(data.frame(cluster = prof$cluster[i],
                        subset_label = NA_character_,
                        subset_markers = NA_character_,
                        n_requirements = 0L, margin = NA_real_,
                        subset_alternatives = NA_character_,
                        stringsAsFactors = FALSE))
    # Most specific first, then the best-separated of equals.
    o <- order(-nreq, -marg)
    best <- hit[o[1]]
    data.frame(
      cluster = prof$cluster[i],
      subset_label = best,
      subset_markers = paste(sprintf("%s %s", names(subsets[[best]]),
                                     ifelse(unlist(subsets[[best]]) == "above",
                                            "+", "-")), collapse = " "),
      n_requirements = nreq[o[1]],
      margin = round(marg[o[1]], 3),
      subset_alternatives = if (length(o) > 1L)
        paste(utils::head(hit[o[-1]], 3), collapse = "; ") else NA_character_,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
