# SECTION 13 -- RUN REPORT
# =============================================================================
#
# WHY THIS FILE EXISTS. A run writes several dozen tables and figures into a
# directory. The reading ORDER is the product: the diagnostics article sets out
# which checks can invalidate which results, and the statistics article sets out
# that a difference is read as adjusted p, then effect size, then against the
# gate's own uncertainty, then against the counting uncertainty. A directory
# listing enforces none of that. It presents a failed staining QC and a headline
# p-value as two files of equal standing, sorted alphabetically.
#
# The field says the same thing about adoption. The blockers reported for
# automated analysis are accessibility and user confidence rather than accuracy
# (Popp et al. 2025, Cytometry A 107:189).
#
# WHY IT IS WRITTEN BY HAND RATHER THAN THROUGH rmarkdown. Rendering Markdown to
# HTML needs pandoc, which is a system binary rather than an R package, and the
# runtime image deliberately does not carry one. A report that only works outside
# the container would be a report the documented execution path cannot produce.
# The HTML here is assembled directly, so it needs nothing that is not already
# present.
#
# WHY FIGURES ARE REFERENCED RATHER THAN EMBEDDED. The report is written into the
# results directory, beside the files it describes, and is meant to be read there.
# Embedding would need a base64 encoder this package does not have and would
# multiply the size of every figure by four to hide a dependency that does not
# exist. Moving the report elsewhere without its directory breaks the images, and
# that is the intended behaviour: the report is an index of a result, not a
# substitute for it.

#' Escape text for inclusion in HTML
#' @param x character vector
#' @keywords internal
html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

#' Render a data.frame as an HTML table
#' @param d a data.frame
#' @param max_rows rows beyond which the table is truncated, with a note
#' @keywords internal
html_table <- function(d, max_rows = 20L) {
  if (is.null(d) || !nrow(d)) return("<p class='none'>No rows.</p>")
  note <- ""
  if (nrow(d) > max_rows) {
    note <- sprintf("<p class='none'>Showing %d of %d rows; the file has the rest.</p>",
                    max_rows, nrow(d))
    d <- utils::head(d, max_rows)
  }
  hdr <- paste0("<th>", html_escape(names(d)), "</th>", collapse = "")
  body <- apply(d, 1, function(r)
    paste0("<tr>", paste0("<td>", html_escape(r), "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead><tr>", hdr, "</tr></thead><tbody>",
         paste(body, collapse = ""), "</tbody></table>", note)
}

#' One report section: a question, its evidence, and where the numbers live
#' @keywords internal
report_section <- function(outdir, title, question, figures = character(0),
                           tables = character(0), body = NULL, max_rows = 20L) {
  fig_present <- figures[file.exists(file.path(outdir, figures))]
  tab_present <- tables[file.exists(file.path(outdir, tables))]
  if (!length(fig_present) && !length(tab_present) && is.null(body)) return("")

  h <- c(sprintf("<section><h2>%s</h2>", html_escape(title)),
         sprintf("<p class='q'>%s</p>", question))
  if (!is.null(body)) h <- c(h, body)
  for (f in fig_present)
    h <- c(h, sprintf("<figure><img src='%s' alt='%s'/><figcaption>%s</figcaption></figure>",
                      f, html_escape(f), html_escape(f)))
  for (t in tab_present) {
    d <- tryCatch(utils::read.csv(file.path(outdir, t), stringsAsFactors = FALSE),
                  error = function(e) NULL)
    h <- c(h, sprintf("<h3><a href='%s'>%s</a></h3>", t, html_escape(t)),
           html_table(d, max_rows = max_rows))
  }
  paste(c(h, "</section>"), collapse = "\n")
}

#' Write a single HTML report indexing everything a run produced
#'
#' Presents the outputs in the order the documentation says they must be read,
#' because each stage can invalidate the ones after it. Sections whose files are
#' absent are omitted rather than shown empty.
#'
#' @param outdir the results directory, which is also where the report is written
#' @param opt parsed options, used for the header
#' @param verdicts per-sample staining verdicts, used for the summary banner
#' @return the path, invisibly, or NULL when there is nothing to report
#' @export
write_run_report <- function(outdir, opt = NULL, verdicts = NULL) {
  if (!dir.exists(outdir)) return(invisible(NULL))
  path <- file.path(outdir, "report.html")

  css <- paste(
    "body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;",
    "max-width:1100px;margin:2rem auto;padding:0 1.5rem;line-height:1.55;color:#1a1a1a}",
    "h1{font-size:1.7rem;margin-bottom:.2rem}h2{font-size:1.25rem;margin-top:2.5rem;",
    "border-bottom:1px solid #ddd;padding-bottom:.3rem}h3{font-size:.95rem;",
    "font-family:ui-monospace,SFMono-Regular,Menlo,monospace;margin-top:1.5rem}",
    "p.q{color:#555;font-style:italic;margin-top:-.2rem}",
    "p.none{color:#777;font-size:.85rem}",
    "table{border-collapse:collapse;width:100%;font-size:.8rem;margin:.5rem 0}",
    "th,td{border:1px solid #ddd;padding:.25rem .5rem;text-align:left}",
    "th{background:#f4f4f4}tr:nth-child(even){background:#fafafa}",
    "figure{margin:1rem 0}img{max-width:100%;border:1px solid #eee}",
    "figcaption{font-size:.75rem;color:#777;font-family:ui-monospace,monospace}",
    "a{color:#0b5}.banner{padding:.8rem 1rem;border-radius:4px;margin:1rem 0}",
    ".ok{background:#eefaf0;border:1px solid #b6e3c2}",
    ".warn{background:#fff6e5;border:1px solid #f0d9a8}",
    ".stop{background:#fdecea;border:1px solid #f5b5ae}",
    sep = "")

  # The banner. A run whose staining QC excluded samples says so before anything
  # else, because every frequency below is computed without them.
  banner <- ""
  if (!is.null(verdicts) && length(verdicts)) {
    st <- vapply(verdicts, function(v) v$qc_status %||% "pass", character(1))
    ctl <- vapply(verdicts, function(v) isTRUE(v$is_control), logical(1))
    nfail <- sum(st == "failed" & !ctl)
    banner <- if (nfail)
      sprintf("<div class='banner stop'><b>%d of %d sample(s) failed staining QC and are excluded from every test below.</b> See staining_qc.csv for the reason each was excluded.</div>",
              nfail, length(st))
    else
      sprintf("<div class='banner ok'>All %d declared sample(s) passed staining QC.</div>",
              sum(!ctl))
  }

  head <- c(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>cyRAVEN run report</title>",
    sprintf("<style>%s</style>", css), "</head><body>",
    "<h1>cyRAVEN run report</h1>",
    sprintf("<p class='q'>%s &middot; cyRAVEN %s</p>",
            html_escape(format(Sys.time(), tz = "UTC", usetz = TRUE)),
            html_escape(tryCatch(as.character(utils::packageVersion("cyRAVEN")),
                                 error = function(e) "unknown"))),
    banner,
    "<div class='banner warn'>Read the sections in the order given. Each one can",
    "invalidate the sections after it, so a result taken from the bottom without",
    "the top is not supported by this run.</div>")

  s <- c(
    report_section(outdir, "1. Did the gates land in the right place?",
      "A threshold placed on a distribution shoulder rather than a density minimum is visible here and in no table.",
      figures = c("recon_diagnostics.png", "gating_qc.png")),

    report_section(outdir, "2. Was the acquisition stable?",
      "A file whose event rate or signal moved part-way through was two instruments over its run, and one threshold suits neither half.",
      figures = "acquisition_qc.png",
      tables = c("acquisition_qc.csv", "acquisition_qc_impact.csv")),

    report_section(outdir, "3. Is the staining usable?",
      "A sample with no resolvable CD45+ mode has no usable parent gate, and its percentages are fractions of an arbitrary slice.",
      tables = "staining_qc.csv"),

    report_section(outdir, "4. Do the populations express what they were declared to express?",
      "Identity is declared here rather than inferred, so this figure can falsify the gate that produced it.",
      figures = c("population_marker_heatmap.png", "cohort_composition_heatmap.png")),

    report_section(outdir, "5. Where did the thresholds come from?",
      "A quantile fallback is not invalid, but it carries no evidence of separation. A marker that falls back across most samples is not resolving in this panel.",
      tables = c("thresholds_used.csv", "threshold_scale_qc.csv",
                 "fmo_agreement.csv", "spreading_receivers.csv"),
      max_rows = 12L),

    report_section(outdir, "6. How well determined are they, and were enough cells counted?",
      "Placement precision and counting sufficiency are separate guarantees. A cut through a wide empty gap is well determined however few events lie beyond it.",
      figures = c("frequency_uncertainty.png", "uncertainty_budget.png",
                  "detection_limits.png"),
      tables = "uncertainty_budget.csv", max_rows = 10L),

    report_section(outdir, "7. What are the populations?",
      "Report pct_of_cd45_pos. count is an event count, set by acquisition duration, and is not a cell number.",
      figures = c("population_frequencies.png", "umap_overview.png"),
      tables = "population_frequencies.csv", max_rows = 12L),

    report_section(outdir, "8. Do the groups differ?",
      "Read adjusted p, then Cliff's delta, then difference_over_gate_u, then difference_over_total_u. A result surviving all four is one the run supports.",
      figures = c("group_comparison.png", "marker_state.png"),
      tables = c("group_comparison_stats.csv", "compositional_concordance.csv"),
      max_rows = 15L),

    report_section(outdir, "9. Is the contrast confounded?",
      "A variable confounds only when it both differs between groups and associates with the outcome. Batch correction is refused where batch and group are not separable.",
      figures = c("threshold_drift.png", "batch_diagnostic.png"),
      tables = c("threshold_drift_stats.csv", "confounding_diagnostics.csv",
                 "batch_group_confounding.csv", "marker_batch_drift.csv"),
      max_rows = 10L),

    report_section(outdir, "10. Does this run agree with the accepted one?",
      "The within-run peer check finds one deviant tube. It cannot find a cohort that moved as a whole, because the peer median moves with it.",
      tables = c("specification_conformance.csv", "specification_changes.csv"),
      max_rows = 12L),

    report_section(outdir, "11. Provenance",
      "What produced this folder.",
      body = paste0(
        "<p>Full detail in <a href='run_manifest.txt'>run_manifest.txt</a>",
        if (file.exists(file.path(outdir, "miflowcyt.md")))
          " and <a href='miflowcyt.md'>miflowcyt.md</a>" else "",
        ".</p>")))

  s <- s[nzchar(s)]
  writeLines(c(head, s, "</body></html>"), path)
  log_msg("wrote report.html (", length(s), " section(s), ",
          "figures referenced from this directory)")
  invisible(path)
}
