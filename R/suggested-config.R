# THE NEXT SPECIFICATION, ASSEMBLED FROM WHAT EXPLORE FOUND
# =============================================================================
#
# WHAT THIS IS FOR. --maybe-learn already tells you two things: which declared
# populations span several clusters, and which clusters nothing in the
# specification covers. Both land in spec_gaps.csv as prose. Acting on them
# meant opening that file, finding each named cluster in
# explore_suggested_spec_<panel>.yaml, and pasting its marker list into the
# config by hand -- and until the vocabulary was fixed the suggested spec was
# written in "pos"/"neg", which no parser in this package reads, so the paste
# produced a population that silently scored nothing.
#
# This assembles the obvious next config instead: the declared specification
# exactly as given, plus one draft entry per uncovered cluster, in the format
# --config actually takes.
#
# WHAT IT IS NOT. It is not a validated specification and nothing in the run
# that writes it uses it. A cluster is not a population: the entries are
# unnamed, they may be two populations or half of one, and several of them are
# usually debris or doublets. The package's position is that a specification is
# a hypothesis declared before the data are examined -- a draft derived FROM the
# data is a weaker thing, and calling it a specification without curating it
# would be the exact circularity explore exists to avoid. Hence the header the
# file carries, and hence "suggested" in its name.

#' Render one population's marker requirements as config YAML lines
#' @param name Population name.
#' @param req Named list of marker -> direction.
#' @param comment Optional comment line placed above the entry.
#' @keywords internal
.pop_yaml_lines <- function(name, req, comment = NULL) {
  if (!length(req)) return(character(0))
  q <- function(x) if (grepl("^[A-Za-z][A-Za-z0-9_]*$", x)) x else paste0('"', x, '"')
  out <- character(0)
  if (!is.null(comment)) out <- c(out, paste0("  # ", comment))
  out <- c(out, paste0("  ", q(name), ":"))
  for (m in names(req)) {
    v <- req[[m]]
    if (is.list(v)) {
      # an any_of block, carried through unchanged
      out <- c(out, "    any_of:")
      for (mm in names(v)) out <- c(out, sprintf("      %s: %s", q(mm), v[[mm]]))
    } else {
      out <- c(out, sprintf("    %s: %s", q(m), as.character(v)))
    }
  }
  out
}

#' Write the declared specification plus a draft entry per uncovered cluster
#'
#' @param cfg The parsed run config, whose `populations` block is copied
#'   through unchanged.
#' @param ex What `run_explore()` returned: `dir` and `gaps`.
#' @param outfile Destination YAML.
#' @return `list(n_declared, n_added, path)`, or NULL when there is nothing to add.
#' @export
write_suggested_config <- function(cfg, ex, outfile) {
  declared <- cfg$populations %||% list()
  gaps <- ex$gaps
  if (is.null(gaps) || !nrow(gaps)) return(NULL)
  miss <- gaps[gaps$issue == "cluster no declared population covers", , drop = FALSE]
  if (!nrow(miss)) return(NULL)

  # The cluster definitions come from the files explore has already written,
  # rather than being recomputed here: one source for what a cluster is, so the
  # config cannot disagree with explore_suggested_spec_<panel>.yaml about it.
  sugg <- list()
  for (f in list.files(ex$dir, pattern = "^explore_suggested_spec.*[.]yaml$",
                       full.names = TRUE)) {
    tag <- sub("^explore_suggested_spec_?", "", sub("[.]yaml$", "", basename(f)))
    y <- try(yaml::yaml.load_file(f), silent = TRUE)
    if (inherits(y, "try-error") || is.null(y$populations)) next
    for (k in names(y$populations))
      sugg[[paste0(tag, "::", k)]] <- y$populations[[k]]
  }

  added <- list(); notes <- character(0)
  for (i in seq_len(nrow(miss))) {
    pan <- as.character(miss$panel[i]); k <- as.character(miss$subject[i])
    key <- paste0(pan, "::", k)
    req <- sugg[[key]] %||% sugg[[paste0("::", k)]]
    if (is.null(req)) next
    # Qualified by panel: k9 in panel_1 and k9 in panel_2 are different cells,
    # and an unqualified name would silently merge them.
    nm <- if (nzchar(pan) && !is.na(pan)) paste0("explore_", pan, "_", k)
          else paste0("explore_", k)
    added[[nm]] <- req
    notes[nm] <- as.character(miss$detail[i])
  }
  if (!length(added)) return(NULL)

  lines <- c(
    "# DRAFT specification, assembled by cyRAVEN from --explore --maybe-learn.",
    "#",
    "# It is the declared specification this run used, unchanged, followed by",
    "# one entry per cluster that spec_gaps.csv reported nothing covered. The",
    "# run that wrote this file did NOT use it.",
    "#",
    "# Before running it:",
    "#   - a cluster is not a population. Merge entries that are one thing,",
    "#     drop the ones that are debris or doublets, and give the rest names.",
    "#   - the explore_* entries are unnamed by design; naming them is the",
    "#     judgement this file cannot make for you.",
    "#   - spec_gaps.csv also lists declared populations that span several",
    "#     clusters. Those are NOT split here, because splitting one declared",
    "#     population into several is a decision about what the populations",
    "#     are, not a gap to be filled mechanically.",
    "#",
    "# Then: --config suggested_config_next_run.yaml, after editing.",
    "",
    "populations:")

  for (nm in names(declared))
    lines <- c(lines, .pop_yaml_lines(nm, declared[[nm]]), "")
  lines <- c(lines,
             "  # ---- drafts, from clusters no declared population covered ----")
  for (nm in names(added))
    lines <- c(lines, .pop_yaml_lines(nm, added[[nm]], comment = notes[[nm]]), "")

  # Everything else in the config travels with it: dropping functional_blocks
  # or ratios would silently shrink the next run's output relative to this one.
  for (blk in setdiff(names(cfg), "populations")) {
    y <- try(yaml::as.yaml(stats::setNames(list(cfg[[blk]]), blk)), silent = TRUE)
    if (!inherits(y, "try-error"))
      lines <- c(lines, "", strsplit(y, "\n", fixed = TRUE)[[1]])
  }

  writeLines(lines, outfile)
  list(n_declared = length(declared), n_added = length(added), path = outfile)
}
