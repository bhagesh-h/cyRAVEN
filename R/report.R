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
# WHY EVERYTHING IS EMBEDDED. The report is one file that carries every figure
# and every table inside it, as data URIs and as JSON. It references nothing.
# That is what makes it the thing you can attach to an email, put in a
# supplement, or archive on its own and still have the whole result: a report
# whose images live beside it becomes a page of broken icons the moment it is
# moved, and a result that cannot survive being moved is not a record.
#
# The cost is size, and it is a real one. Figures are embedded once, at full
# resolution, and displayed scaled down by CSS; the same bytes serve the
# on-screen figure and the full-resolution download, so nothing is stored twice.
# The final size is logged, because a 40 MB HTML file should not be a surprise.
#
# WHY THE TABLES ARE JSON RATHER THAN MARKUP. Every table is searchable, can be
# paged at 10/50/100/all rows, and can be exported to CSV. Doing that over
# pre-rendered <tr> elements means the export re-parses the DOM and the search
# hides rather than filters. Carrying the data and rendering from it makes all
# three operations read the same array, so what you export is what you filtered.

#' Size above which a table is named rather than embedded
#'
#' Set with `options(cyRAVEN.report_table_max_mb = )`. The default of 8 MB sits
#' above every summary table a run writes and below the per-cell exports, whose
#' row count is the number of cells rather than the number of samples.
#' @keywords internal
report_table_max_bytes <- function() {
  mb <- getOption("cyRAVEN.report_table_max_mb", 8)
  if (!is.numeric(mb) || length(mb) != 1L || is.na(mb) || mb <= 0) mb <- 8
  mb * 1024^2
}

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

#' Escape a string for embedding inside a <script> block as JSON
#'
#' `</script>` anywhere inside the data would end the block early, so the slash
#' of any closing tag is escaped. JSON treats `\/` as `/`, so the data is
#' unchanged by it.
#' @param x character vector.
#' @keywords internal
json_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\r", " ", x, fixed = TRUE)
  x <- gsub("\n", " ", x, fixed = TRUE)
  x <- gsub("\t", " ", x, fixed = TRUE)
  # Control characters JSON forbids unescaped.
  x <- gsub("[\001-\037]", " ", x)
  gsub("</", "<\\/", x, fixed = TRUE)
}

#' Render a data.frame as a JSON array-of-arrays with a header
#' @param d a data.frame
#' @keywords internal
json_table <- function(d) {
  cols <- paste0("[", paste0("\"", json_escape(names(d)), "\"", collapse = ","), "]")
  if (!nrow(d)) return(paste0("{\"cols\":", cols, ",\"rows\":[]}"))
  cells <- lapply(d, function(col) paste0("\"", json_escape(col), "\""))
  rows <- do.call(paste, c(cells, sep = ","))
  paste0("{\"cols\":", cols, ",\"rows\":[",
         paste0("[", rows, "]", collapse = ","), "]}")
}

#' One line saying what a table contains
#'
#' WHY IT IS A LOOKUP AND NOT DERIVED. A column list is not a description: a
#' reader who meets `compositional_concordance.csv` needs to be told that it
#' classifies each result against both parameterisations, and no amount of
#' inspecting its headers says that. The lines below are the same one-liners the
#' `outputs` vignette carries, so the report and the documentation cannot drift
#' into saying different things about the same file.
#'
#' A run with several marker panels writes `group_comparison_stats_<panel>.csv`,
#' so a name that is not found is tried once more with its last underscore segment
#' removed. An unknown table gets no line rather than a vague one -- a description
#' that could be about any table is worse than none, and its absence is a visible
#' prompt to add it here.
#' @param file table file name, with or without a directory.
#' @return A single string, or "" when the table is not known.
#' @keywords internal
report_table_note <- function(file) {
  d <- c(
    # ---- quality control ----
    "staining_qc.csv" = "Per-sample staining verdict and the reason for any exclusion.",
    "thresholds_used.csv" = "Every threshold this run derived: value, how it was derived, the cofactor, and whether the sample is an outlier against its peers.",
    "threshold_scale_qc.csv" = "Per marker, the panel median threshold and each sample's robust z against it.",
    "acquisition_qc.csv" = "Per sample, whether the instrument behaved consistently through the acquisition, and which intervals were flagged.",
    "acquisition_qc_bins.csv" = "One row per time interval per sample: event rate and how far each channel departed from the run.",
    "acquisition_qc_impact.csv" = "How far each population's frequency would move if the flagged acquisition intervals were dropped.",
    "fmo_agreement.csv" = "Each derived cut against its FMO-anchored equivalent, expressed in units of the cut's own uncertainty.",
    "spreading_pairs.csv" = "Per channel pair, how much wider the receiving channel's negative population becomes.",
    "spreading_receivers.csv" = "Per marker, the total spreading it receives and whether that explains a threshold fallback.",
    "calibration.csv" = "Bead calibration fit per channel, written by --calibration-beads.",
    "gate_counts.csv" = "Event counts surviving each level of the gate hierarchy, per sample.",
    "populations_unavailable.csv" = "One row per sample and population the panel cannot score, with the reason it cannot.",
    # ---- uncertainty ----
    "threshold_uncertainty.csv" = "Per sample and marker: the sampling and method components of the threshold's uncertainty, their quadrature sum, and how well the valley resolved.",
    "uncertainty_budget.csv" = "How much each threshold contributes to each population's frequency uncertainty, per sample.",
    # ---- conformance ----
    "specification_conformance.csv" = "This run's thresholds against an accepted baseline's, scaled by the baseline's own spread, with a verdict per marker.",
    "specification_conformance_populations.csv" = "The same baseline comparison for population frequencies rather than thresholds.",
    "specification_changes.csv" = "Populations added, removed or redefined since the baseline run.",
    # ---- abundance ----
    "population_frequencies.csv" = "The run's primary result: one row per sample and population, with the percentage of parent and the event count behind it.",
    "population_marker_mfi.csv" = "Median transformed intensity and percent positive, per sample, population and marker.",
    "functional_markers.csv" = "Percent positive for the functional marker blocks declared in the config, per sample and population.",
    "population_ratios.csv" = "Derived population ratios, where the config declares any.",
    "absolute_counts.csv" = "Measured cells per microlitre per sample and population, from --absolute-counts.",
    "absolute_counts_raw.csv" = "The supplied counts workbook flattened exactly as read, before any header or name interpretation. Read this first when a count looks wrong.",
    "absolute_counts_stats.csv" = "Between-group tests on measured cells per microlitre rather than on percentages.",
    # ---- inference ----
    "group_comparison_stats.csv" = "Between-group abundance test per population: medians, Cliff's delta, raw and BH-adjusted p, and the uncertainty the difference is worth in.",
    "compositional_clr_stats.csv" = "The same abundance tests run on centred log-ratios, because percentages of one parent cannot all move independently.",
    "compositional_concordance.csv" = "Each population classified by whether the percentage and the log-ratio parameterisations agree about it.",
    "marker_state_stats.csv" = "Differential marker state per population and marker, on per-sample medians.",
    "functional_markers_stats.csv" = "Between-group tests on functional marker positivity.",
    "population_ratios_stats.csv" = "Between-group tests on the declared population ratios.",
    "paired_comparison_stats.csv" = "Paired test per population for a repeated-measures design, using only the units present in both conditions.",
    "covariate_adjusted_stats.csv" = "Rank ANCOVA: the group difference after adjusting for the named covariates.",
    "subcluster_marker_shifts.csv" = "Marker shifts between compartments on pooled events. Descriptive, with effect sizes and no p-values.",
    "design_feasibility.csv" = "Which group comparisons this design can support and why the rest cannot. Read before any p-value.",
    "parametric_tests.csv" = "The t-test or ANOVA equivalent of each rank test, with the normality and equal-variance assumptions it needs recorded beside it.",
    "posthoc_tests.csv" = "Pairwise comparisons by Games-Howell, Tukey HSD and Dunn, for three or more groups.",
    "normality_tests.csv" = "Shapiro-Wilk per population and group, and Brown-Forsythe across groups: the evidence for using rank tests.",
    "statistical_methods.csv" = "Every method commonly reported in this literature, whether this run computed it, and why.",
    # ---- clinical ----
    "clinical_association.csv" = "One row per population and clinical variable: the test, n, the effect with its bootstrap interval, raw and BH-adjusted p, and an underpowered flag.",
    "clinical_association_markers.csv" = "The same clinical association against per-sample median marker intensity, collapsed across populations.",
    # ---- diagnostics ----
    "threshold_drift_stats.csv" = "Whether the thresholds themselves separate by group, which would make an abundance difference definitional rather than biological.",
    "confounding_diagnostics.csv" = "Per covariate: whether it differs between groups, whether it associates with the outcome, and the verdict from both.",
    "batch_mixing_stats.csv" = "iLISI batch mixing in the embedding against a permutation null.",
    "batch_group_confounding.csv" = "Cramer's V between batch and study group, and whether batch correction was allowed or refused on it.",
    "marker_batch_drift.csv" = "Earth Mover's distance between batches per marker, in analysis units and scaled by the marker's own MAD.",
    "threshold_batch_drift.csv" = "The threshold drift test grouped by acquisition batch instead of by study group.",
    "batch_correction.csv" = "Whether a batch correction ran, which method fitted it, and the reason. Written on a refusal too.",
    # ---- clustering ----
    "unsupervised_clusters.csv" = "Cluster assignment per cell, from the unsupervised cross-check on the gate specification.",
    "cluster_gate_agreement_clusters.csv" = "Per cluster, the declared population that dominates it and how pure that match is.",
    "cluster_gate_agreement_populations.csv" = "Per declared population, how much of it the unsupervised clustering recovers.",
    "subcluster_k_selection.csv" = "The silhouette-selected number of subclusters, from --auto-subcluster-k.",
    "cluster_gate_proposals.csv" = "Learned two-marker gate geometry for clusters no declared population covers, with held-out metrics.",
    "cluster_gate_polygons.csv" = "Polygon vertices of those proposed gates, in transformed units.",
    # ---- external labels ----
    "external_label_gates.csv" = "Per supplied external label, the learned gating strategy and its held-out metrics at each depth.",
    "external_label_polygons.csv" = "Polygon vertices of the learned strategies, on the analysis scale.",
    "gate_transferability.csv" = "Precision, recall and F1 on each donor, from a strategy refitted with that donor withheld.",
    "gate_transferability_summary.csv" = "Minimum, median, maximum and IQR of F1 across donors.",
    # ---- explore ----
    "explore_cluster_profile.csv" = "Per cluster: size, phenotype string, fraction positive and median per channel.",
    "explore_qc_clusters.csv" = "The cluster-level gate: the call for each cluster and the basis for it.",
    "explore_cluster_abundance.csv" = "Per donor and cluster, with counting uncertainty and the limits of detection and quantification.",
    "explore_cluster_stats.csv" = "Donor-level group tests on cluster abundance, carrying the batch and group confounding verdict.",
    "explore_cells.csv" = "One row per cell: sample, event index, cluster and UMAP coordinates.",
    "explore_vs_populations.csv" = "Cross-tabulation of unsupervised clusters against the declared population labels.",
    "explore_findings.csv" = "Clusters that no declared population covers -- what the specification is missing.",
    "explore_population_split.csv" = "Declared labels that span several clusters -- what the specification lumps together.",
    "explore_provenance.csv" = "Every choice explore mode made, and the basis for each.",
    "spec_gaps.csv" = "Gaps between the population specification and what the data contains.",
    # ---- auxiliary ----
    "patient_metadata_english.csv" = "The patient table after column mapping and value translation, as the pipeline read it.",
    "population_codes.csv" = "Numeric codes for populations, samples and cohorts, for writing back into FCS keywords.",
    "sample_map_template.csv" = "A sample map skeleton for this directory, written by --write-sample-map.",
    "cells_umap.csv" = "One row per cell: embedding coordinates with sample, population and panel.")
  f <- basename(file)
  if (f %in% names(d)) return(unname(d[f]))
  # Panel suffix: group_comparison_stats_panel2.csv -> group_comparison_stats.csv.
  g <- sub("_[^_]+\\.csv$", ".csv", f)
  if (g != f && g %in% names(d))
    return(paste0(unname(d[g]), " This copy covers one marker panel only."))
  ""
}

#' Legacy HTML table renderer
#'
#' Retained because the interactive tables are rendered from JSON in the
#' browser; this is the static fallback used when a table is too small for
#' searching to be worth the controls.
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

#' The package logo as a data URI, or NA when it cannot be found
#'
#' An installed package keeps `man/figures/` at `help/figures/`, while a source
#' tree loaded with pkgload keeps the original path, so both are tried rather
#' than assumed. Returns NA rather than failing: a report is worth writing
#' without its logo, and this runs at the end of an analysis that may have taken
#' an hour.
#' @keywords internal
report_logo_uri <- function() {
  for (p in c(system.file("help", "figures", "logo.png", package = "cyRAVEN"),
              system.file("man", "figures", "logo.png", package = "cyRAVEN"))) {
    if (nzchar(p) && file.exists(p)) return(file_data_uri(p, "image/png"))
  }
  NA_character_
}

#' The per-marker embeddings, pooled first and then split by category
#'
#' WHY A HELPER. A run writes one figure per marker and then one per marker per
#' category, so the count is a property of the cohort rather than something a
#' section can name in advance. Returned as paths relative to `outdir`, which is
#' what [report_section()] and the catch-all both work in.
#'
#' Ordering is deliberate: the pooled view of a marker comes before its splits,
#' so the reader sees the whole embedding before any comparison drawn on it.
#' @param outdir The run directory.
#' @keywords internal
marker_umap_files <- function(outdir) {
  d <- file.path(outdir, "marker_umaps_by_group")
  if (!dir.exists(d)) return(character(0))
  f <- list.files(d, "[.]png$")
  if (!length(f)) return(character(0))
  split_by <- grepl("_by_", f, fixed = TRUE)
  file.path("marker_umaps_by_group", c(sort(f[!split_by]), sort(f[split_by])))
}

#' Build one report section, embedding its figures and tables
#'
#' @param outdir the results directory.
#' @param id anchor id, unique per section.
#' @param title section heading, stated as what the section reports.
#' @param description what the figures and tables below show and how to read them.
#' @param figures figure filenames, embedded when present.
#' @param tables table filenames, embedded when present.
#' @param body optional raw HTML inserted before the figures.
#' @param open whether the section starts expanded.
#' @return list(html, nav, n_fig, n_tab, bytes, used)
#' @keywords internal
report_section <- function(outdir, id, title, description, figures = character(0),
                           tables = character(0), body = NULL, open = FALSE) {
  fig_present <- figures[file.exists(file.path(outdir, figures))]
  tab_present <- tables[file.exists(file.path(outdir, tables))]

  # A per-cell export carries one row per CELL, so on a real cohort it is
  # hundreds of thousands of rows. It is read by software rather than by a
  # person scrolling a report, and embedding it would multiply the file size for
  # something nobody reads there. Such a table is NAMED with its size instead of
  # being dropped, so its absence is a stated fact rather than a silent gap.
  oversize <- character(0)
  if (length(tab_present)) {
    big <- vapply(tab_present, function(t)
      file.size(file.path(outdir, t)) > report_table_max_bytes(), logical(1))
    oversize <- tab_present[big]
    tab_present <- tab_present[!big]
  }

  # A section with nothing to show is omitted rather than rendered empty: an
  # empty heading reads as "this check found nothing", which is a different
  # claim from "this check did not run".
  if (!length(fig_present) && !length(tab_present) && !length(oversize) &&
      is.null(body))
    return(list(html = "", nav = "", n_fig = 0L, n_tab = 0L, bytes = 0,
                used = character(0)))

  h <- c(sprintf("<section class='sec' id='%s'>", id),
         sprintf("<details%s><summary><span class='chev' aria-hidden='true'></span>",
                 if (open) " open" else ""),
         sprintf("<span class='sec-title'>%s</span>", html_escape(title)),
         sprintf("<span class='sec-count'>%s</span></summary>",
                 html_escape(paste0(length(fig_present), " fig / ",
                                    length(tab_present), " tab"))),
         "<div class='sec-body'>",
         sprintf("<p class='q'>%s</p>", description))
  if (!is.null(body)) h <- c(h, body)

  nav <- character(0)
  bytes <- 0
  for (f in fig_present) {
    p <- file.path(outdir, f)
    uri <- file_data_uri(p, "image/png")
    if (is.na(uri)) next
    bytes <- bytes + file.size(p)
    # The id keeps the full relative path so two figures with the same file name
    # in different directories cannot collide. Everything the READER sees is the
    # base name: a caption reading marker_umaps_by_group/umap_CD3.png tells them
    # about this report's directory layout, which is not what they are looking
    # at, and the same string was being used as the download file name.
    fid <- paste0("fig-", gsub("[^A-Za-z0-9]+", "-", sub("[.]png$", "", f)))
    fb  <- basename(f)
    h <- c(h, sprintf(paste0(
      "<figure class='fig' id='%s'>",
      "<div class='figbox'><img src='%s' alt='%s' loading='lazy' ",
      "onclick='cyZoom(this)'/>",
      "<button class='zoom' title='Zoom' onclick='cyZoom(this.parentNode.querySelector(\"img\"))'>&#9974;</button>",
      "</div>",
      "<figcaption><span class='fn'>%s</span>",
      "<a class='dl' download='%s' href='%s'>Full resolution PNG</a>",
      "</figcaption></figure>"),
      fid, uri, html_escape(fb), html_escape(fb), html_escape(fb), uri))
    nav <- c(nav, sprintf("<a class='nav-fig' href='#%s'>%s</a>", fid,
                          html_escape(sub("[.]png$", "", fb))))
  }

  for (t in tab_present) {
    p <- file.path(outdir, t)
    d <- tryCatch(utils::read.csv(p, stringsAsFactors = FALSE,
                                  check.names = FALSE),
                  error = function(e) NULL)
    if (is.null(d)) next
    bytes <- bytes + file.size(p)
    tid <- paste0("tab-", gsub("[^A-Za-z0-9]+", "-", sub("[.]csv$", "", t)))
    note <- report_table_note(t)
    # EACH TABLE IS ITS OWN TOGGLE, CLOSED. A section can carry a dozen tables of
    # several hundred rows each, and opening the section to look at one figure
    # used to unroll all of them: the figures ended up separated by yards of
    # scrollable grid. The name and the row count stay visible when it is closed,
    # so the section still says what it holds, and the description is inside --
    # a reader deciding whether to open it needs the name, and a reader who has
    # opened it needs the sentence.
    #
    # It is also why the body renders lazily: cyInit() fills only the tables whose
    # toggle is open, and the rest are built the first time they are opened. On a
    # run with forty tables that is most of the report's load time.
    h <- c(h, sprintf(paste0(
      "<div class='tab' id='%s'>",
      "<details class='tabdet'><summary>",
      "<span class='chev' aria-hidden='true'></span>",
      "<h3>%s</h3><span class='dim'>%d rows &times; %d cols</span></summary>",
      "%s",
      "<div class='tabhead'>",
      "<input class='search' type='search' placeholder='Search this table' ",
      "oninput='cyFilter(\"%s\")' aria-label='Search %s'/>",
      "<select class='pagesel' onchange='cyFilter(\"%s\")' aria-label='Rows to show'>",
      "<option value='10'>10 rows</option><option value='50' selected>50 rows</option>",
      "<option value='100'>100 rows</option><option value='0'>All rows</option></select>",
      "<button class='dl' onclick='cyCsv(\"%s\")'>Export CSV</button>",
      "</div><div class='tabwrap'><table></table></div>",
      "<p class='none shown'></p></details></div>"),
      tid, html_escape(basename(t)), nrow(d), ncol(d),
      if (nzchar(note)) sprintf("<p class='tabnote'>%s</p>", html_escape(note)) else "",
      tid, html_escape(basename(t)), tid, tid))
    h <- c(h, sprintf("<script type='application/json' id='%s-data'>%s</script>",
                      tid, json_table(d)))
    nav <- c(nav, sprintf("<a class='nav-tab' href='#%s'>%s</a>", tid,
                          html_escape(basename(t))))
  }

  for (t in oversize) {
    p <- file.path(outdir, t)
    nr <- tryCatch(length(readLines(p, warn = FALSE)) - 1L,
                   error = function(e) NA_integer_)
    h <- c(h, sprintf(paste0(
      "<p class='none'><b>%s</b> is in the results directory but not embedded ",
      "here: %s row(s), %s MB. It is a per-cell export, one row per cell, read ",
      "by software rather than read in a report; carrying it would multiply the ",
      "size of this file. Every summary derived from it is embedded above.</p>"),
      html_escape(basename(t)), format(nr, big.mark = ","),
      format(round(file.size(p) / 1024^2, 1), nsmall = 1)))
  }

  h <- c(h, "</div></details></section>")
  list(html = paste(h, collapse = "\n"),
       # The per-figure and per-table links live in their own container so the
       # sidebar can show section titles alone until a section is opened. A
       # twelve-section run lists over a hundred entries otherwise, and the
       # structure of the document -- which is what a contents list is for --
       # disappears into them.
       nav = paste0(sprintf(paste0("<div class='nav-sec' data-sec='%s'>",
                                   "<a class='nav-h' href='#%s'>%s</a>"),
                            id, id, html_escape(title)),
                    "<div class='nav-kids'>", paste(nav, collapse = ""),
                    "</div></div>"),
       n_fig = length(fig_present), n_tab = length(tab_present), bytes = bytes,
       # `oversize` belongs in `used`. It was handled by this section -- named
       # with its size rather than embedded -- which is exactly what "covered"
       # means for a table too large to read here. Leaving it out sent every
       # per-cell export to the catch-all as well, so cells_umap.csv was named
       # twice in every report and section 11 was never empty, which destroyed
       # the one invariant that section exists to provide: if it has anything in
       # it, a named section is missing an entry.
       used = c(fig_present, tab_present, oversize))
}

#' Stylesheet for the run report
#' @keywords internal
report_css <- function() {
  paste0(
  ":root{--fg:#1a1a1a;--mut:#5b6470;--line:#e2e5ea;--bg:#fff;--panel:#f7f8fa;",
  "--accent:#0a7d4a;--warn:#8a6100;--stop:#a4231c;",
  "--ui:ui-sans-serif,system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;",
  "--mono:ui-monospace,SFMono-Regular,'SF Mono',Menlo,Consolas,monospace}",
  # One font stack for every element, set once on the root and inherited. Tables,
  # buttons, inputs and selects do NOT inherit font by default in any browser,
  # so they are named explicitly rather than left to the user-agent stylesheet.
  "*{box-sizing:border-box}",
  "html{font-family:var(--ui);font-size:15px;color:var(--fg);background:var(--bg)}",
  "body{margin:0;line-height:1.55;font-family:var(--ui)}",
  "button,input,select,table,th,td{font-family:var(--ui);font-size:inherit;color:inherit}",
  "code,kbd,.mono,.fn,.dim{font-family:var(--mono)}",
  # Layout: fixed sidebar, scrolling main.
  ".wrap{display:flex;align-items:flex-start;min-height:100vh}",
  "nav.side{position:sticky;top:0;flex:0 0 268px;height:100vh;overflow-y:auto;",
  "border-right:1px solid var(--line);background:var(--panel);padding:1rem .75rem}",
  "nav.side h2{font-size:.78rem;text-transform:uppercase;letter-spacing:.08em;",
  "color:var(--mut);margin:.2rem 0 .6rem .4rem}",
  ".nav-sec{margin-bottom:.35rem}",
  "nav.side a{display:block;text-decoration:none;color:var(--fg);border-radius:5px;",
  "padding:.22rem .45rem;font-size:.84rem;overflow:hidden;text-overflow:ellipsis;",
  "white-space:nowrap}",
  "nav.side a:hover{background:#e8ebf0}",
  "a.nav-h{font-weight:600;font-size:.87rem}",
  # The entry for the section being read. Without it the sidebar looks identical
  # at the top of section 1 and the bottom of section 9, so a twelve-section
  # report gives no cue where you are or how much is left. Tinted fill plus an
  # accent rule down the left edge, which survives being scrolled past at a
  # glance better than colour alone.
  ".nav-sec.cur>a.nav-h{background:#e6f1eb;color:var(--accent);",
  "box-shadow:inset 3px 0 0 var(--accent)}",
  # The contents list shows section titles only until a section is opened, and
  # then shows that section's figures and tables. It mirrors the document rather
  # than duplicating it: over a hundred entries at once buries the twelve-line
  # structure a reader actually navigates by.
  ".nav-kids{display:none}",
  ".nav-sec.open>.nav-kids{display:block}",
  # A caret on the title, turned down when the section is open, so the sidebar
  # says which sections are expanded even where the entry is scrolled past.
  "a.nav-h:before{content:'\\25B8';color:var(--mut);display:inline-block;",
  "width:.85em;transition:transform .12s ease}",
  ".nav-sec.open>a.nav-h:before{transform:rotate(90deg)}",
  "a.nav-fig,a.nav-tab{padding-left:1.35rem;color:var(--mut);font-family:var(--mono);",
  "font-size:.76rem}",
  "a.nav-fig:before{content:'\\25A6\\00a0';color:var(--accent)}",
  "a.nav-tab:before{content:'\\2263\\00a0';color:var(--mut)}",
  "main{flex:1 1 auto;min-width:0;max-width:1180px;padding:1.5rem 2rem 5rem}",
  "h1{font-size:1.6rem;margin:0 0 .15rem}",
  # Masthead: mark, identity, document-level actions, on one rule.
  "header.masthead{display:flex;align-items:flex-start;gap:1.3rem;",
  "flex-wrap:wrap;padding:0 0 1.1rem;margin:0 0 .4rem;",
  "border-bottom:2px solid var(--line)}",
  "img.brand{width:104px;height:104px;flex:0 0 auto;object-fit:contain}",
  ".mast-id{flex:1 1 340px;min-width:0}",
  ".mast-id h1{margin:0 0 .1rem;line-height:1.1}",
  "p.meta{color:var(--mut);margin:0;font-size:.86rem}",
  # The note sits under the run line inside the same column, so it reads as part
  # of the masthead rather than as the first content of the document.
  ".mast-note{margin:.6rem 0 0}",
  ".mast-actions{flex:0 0 auto;display:flex;gap:.45rem;align-items:center}",
  # Collapsible sections.
  "section.sec{border:1px solid var(--line);border-radius:8px;margin:.7rem 0;",
  "background:var(--bg);overflow:hidden}",
  "summary{cursor:pointer;padding:.7rem .9rem;display:flex;align-items:center;",
  "gap:.6rem;background:var(--panel);user-select:none;list-style:none}",
  "summary::-webkit-details-marker{display:none}",
  "summary:hover{background:#eef1f5}",
  ".chev{width:0;height:0;border-left:6px solid var(--mut);",
  "border-top:4.5px solid transparent;border-bottom:4.5px solid transparent;",
  "transition:transform .15s;flex:0 0 auto}",
  # DIRECT CHILD, NOT ANY DESCENDANT. Tables are now toggles of their own, so a
  # bare `details[open] .chev` rotated every table's chevron as soon as the
  # SECTION containing them opened -- a dozen closed tables all showing an open
  # arrow.
  "details[open] > summary .chev{transform:rotate(90deg)}",
  ".sec-title{font-weight:600;flex:1 1 auto}",
  ".sec-count{color:var(--mut);font-size:.76rem;font-family:var(--mono)}",
  ".sec-body{padding:.4rem 1rem 1.1rem}",
  "p.q{color:var(--mut);margin:.4rem 0 1rem}",
  "p.none{color:var(--mut);font-size:.82rem}",
  # Figures: every figure occupies the same box whatever its native aspect
  # ratio, so the page does not jump between a wide strip and a tall grid.
  "figure.fig{margin:1rem 0 1.4rem}",
  ".figbox{position:relative;height:420px;border:1px solid var(--line);",
  "border-radius:6px;background:var(--panel);display:flex;align-items:center;",
  "justify-content:center;overflow:hidden}",
  ".figbox img{max-width:100%;max-height:100%;width:auto;height:auto;",
  "object-fit:contain;cursor:zoom-in;display:block}",
  "button.zoom{position:absolute;top:.4rem;right:.4rem;border:1px solid var(--line);",
  "background:rgba(255,255,255,.92);border-radius:5px;cursor:pointer;",
  "padding:.15rem .4rem;font-size:.95rem;line-height:1}",
  "button.zoom:hover{background:#fff}",
  "figcaption{display:flex;justify-content:space-between;align-items:center;",
  "gap:1rem;margin-top:.4rem;font-size:.78rem;color:var(--mut)}",
  "a.dl,button.dl{font-size:.76rem;color:var(--accent);text-decoration:none;",
  "border:1px solid var(--accent);border-radius:5px;padding:.16rem .5rem;",
  "background:none;cursor:pointer;white-space:nowrap}",
  "a.dl:hover,button.dl:hover{background:var(--accent);color:#fff}",
  # Lightbox.
  # --lbbar-h is shared by the toolbar's height and the image's top margin, so
  # the two cannot drift apart.
  "#lb{position:fixed;inset:0;background:rgba(12,14,18,.93);display:none;",
  "z-index:99;overflow:auto;cursor:zoom-out;--lbbar-h:2.6rem}",
  "#lb.on{display:block}",
  # THE IMAGE STARTS BELOW THE TOOLBAR, NOT UNDER IT. The toolbar is fixed so it
  # stays reachable while a zoomed image is scrolled, which also means it is
  # painted over the image rather than above it. Every figure carries its title
  # and subtitle as the first thing inside the PNG, so with the image at y=0 the
  # bar covered exactly the line saying what the figure is and what to read from
  # it -- the reader had to scroll a zoomed image up to identify it. A top
  # margin of the bar's own height clears it at the default view, and the bar
  # still overlays the middle of a scrolled image, which is what a fixed toolbar
  # is for.
  "#lb img{display:block;margin:var(--lbbar-h) auto 0;max-width:none}",
  "#lbbar{position:fixed;top:0;left:0;right:0;height:var(--lbbar-h);",
  "box-sizing:border-box;padding:.5rem .9rem;display:flex;",
  "gap:.5rem;align-items:center;background:rgba(12,14,18,.86);color:#eef1f5;",
  "font-size:.8rem;z-index:100}",
  "#lbbar button{background:none;color:#eef1f5;border:1px solid #556;",
  "border-radius:5px;padding:.15rem .55rem;cursor:pointer}",
  "#lbbar button:hover{background:#2a2f38}",
  "#lbname{font-family:var(--mono);flex:1 1 auto;overflow:hidden;",
  "text-overflow:ellipsis;white-space:nowrap}",
  # Tables. Each is a toggle inside its section, so its summary is styled a step
  # quieter than a section's: a thinner bar, no background of its own until
  # hovered, and the file name in the mono face it is referred to by everywhere
  # else.
  ".tab{margin:.55rem 0}",
  ".tabdet{border:1px solid var(--line);border-radius:6px;overflow:hidden}",
  ".tabdet > summary{padding:.4rem .6rem;background:var(--bg);gap:.5rem}",
  ".tabdet > summary:hover{background:var(--panel)}",
  ".tabdet[open] > summary{background:var(--panel);",
  "border-bottom:1px solid var(--line)}",
  ".tabdet summary h3{margin:0;font-size:.82rem;font-family:var(--mono);",
  "font-weight:600;text-transform:none;letter-spacing:0;color:var(--fg);",
  "flex:0 1 auto}",
  # The one-line description of what the table holds, between its name and the
  # grid itself.
  ".tabnote{color:var(--mut);font-size:.79rem;margin:.5rem .6rem .1rem;",
  "max-width:70ch}",
  ".tabhead{display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;",
  "margin:.45rem .6rem .35rem}",
  ".tabdet > .tabwrap,.tabdet > p.shown{margin-left:.6rem;margin-right:.6rem}",
  ".tabdet > p.shown{margin-bottom:.5rem}",
  ".dim{color:var(--mut);font-size:.74rem;flex:1 1 auto}",
  "input.search,select.pagesel{border:1px solid var(--line);border-radius:5px;",
  "padding:.2rem .45rem;font-size:.78rem;background:#fff}",
  "input.search{width:12rem}",
  ".tabwrap{max-height:26rem;overflow:auto;border:1px solid var(--line);",
  "border-radius:6px}",
  ".tabwrap table{border-collapse:collapse;width:100%;font-size:.78rem}",
  ".tabwrap th,.tabwrap td{border-bottom:1px solid var(--line);padding:.26rem .5rem;",
  "text-align:left;white-space:nowrap}",
  ".tabwrap th{background:var(--panel);position:sticky;top:0;font-weight:600;",
  "cursor:pointer}",
  ".tabwrap th:hover{background:#e8ebf0}",
  ".tabwrap tr:nth-child(even) td{background:#fbfcfd}",
  # Banners.
  # Diagnosis block on a failed run.
  "h3{font-size:.9rem;margin:1.1rem 0 .3rem;text-transform:uppercase;",
  "letter-spacing:.05em;color:var(--mut)}",
  ".sec-body h3:first-of-type{margin-top:.2rem}",
  "pre{font-family:var(--mono);font-size:.79rem;line-height:1.45;overflow-x:auto;",
  "border:1px solid var(--line);border-radius:6px;padding:.6rem .75rem;",
  "background:var(--panel);margin:.2rem 0;white-space:pre-wrap;word-break:break-word}",
  "pre.err{background:#fdecea;border-color:#f5b5ae;color:#7d1a15;white-space:pre-wrap}",
  "pre.log{max-height:22rem;overflow-y:auto;white-space:pre}",
  "pre.cmd{background:#14181d;border-color:#14181d;color:#e6edf3}",
  ".banner{padding:.7rem .9rem;border-radius:6px;margin:.8rem 0;font-size:.9rem}",
  ".ok{background:#eefaf0;border:1px solid #b6e3c2}",
  ".warn{background:#fff6e5;border:1px solid #f0d9a8}",
  ".stop{background:#fdecea;border:1px solid #f5b5ae}",
  "@media print{nav.side{display:none}details{open:true}}",
  # The report is read on laptops and projected in meetings; below 900px the
  # sidebar becomes a normal block above the content rather than disappearing.
  "@media(max-width:900px){.wrap{display:block}nav.side{position:static;height:auto;",
  "width:auto;border-right:none;border-bottom:1px solid var(--line)}",
  "main{padding:1rem}.figbox{height:300px}}")
}

#' Behaviour for the run report
#' @keywords internal
report_js <- function() {
  paste0(
  "var CY={};\n",
  # Tables are parsed once on load from their JSON blocks and kept in memory,
  # so search, paging and export all read the same array.
  "function cyInit(){document.querySelectorAll(\"script[type='application/json']\")",
  ".forEach(function(s){CY[s.id.replace(/-data$/,'')]=JSON.parse(s.textContent);});",
  # RENDER ONLY WHAT IS OPEN. Every table is a closed toggle, so building all of
  # their grids at load would spend the time and then hide the result. cyShow()
  # builds one the first time it is opened and marks it, so reopening is free.
  "Object.keys(CY).forEach(function(id){CY[id].sort=-1;",
  "var b=document.getElementById(id);var t=b&&b.querySelector('details');",
  "if(!t||t.open)cyShow(id);});",
  "cySync();cyCur();}\n",
  "function cyShow(id){if(!CY[id]||CY[id].drawn)return;CY[id].drawn=1;cyFilter(id);}\n",
  # A table link in the contents opens the table it points at. Without this the
  # browser jumps to a collapsed toggle and the reader sees the name they clicked
  # and no table, which reads as a broken link.
  "function cyReveal(el){for(var n=el;n;n=n.parentNode){",
  "if(n.tagName==='DETAILS'&&!n.open)n.open=true;}",
  "var b=el.closest?el.closest('.tab'):null;if(b)cyShow(b.id);}\n",
  "document.addEventListener('click',function(e){",
  "var a=e.target&&e.target.closest?e.target.closest('a.nav-tab,a.nav-fig'):null;",
  "if(!a)return;var h=a.getAttribute('href')||'';if(h.charAt(0)!=='#')return;",
  "var el=document.getElementById(h.slice(1));if(!el)return;cyReveal(el);",
  "setTimeout(function(){el.scrollIntoView({block:'start'});},0);});\n",
  "function cyRows(id){var d=CY[id];var box=document.getElementById(id);",
  "var q=box.querySelector('input.search').value.trim().toLowerCase();",
  "var rows=d.rows;",
  "if(q){var terms=q.split(/\\s+/);rows=rows.filter(function(r){",
  "var hay=r.join(' ').toLowerCase();",
  "return terms.every(function(t){return hay.indexOf(t)>-1;});});}",
  # Sorting is numeric when every value in the column parses as a number, and
  # lexical otherwise. Deciding per column rather than per cell stops a column
  # of p-values sorting as text because one cell reads "NA".
  "if(d.sort>=0){var c=d.sort,dir=d.dir||1;",
  "var num=rows.every(function(r){return r[c]===''||!isNaN(parseFloat(r[c]));});",
  "rows=rows.slice().sort(function(a,b){var x=a[c],y=b[c];",
  "if(num){x=parseFloat(x);y=parseFloat(y);",
  "if(isNaN(x))return 1;if(isNaN(y))return -1;return (x-y)*dir;}",
  "return x.localeCompare(y)*dir;});}",
  "return rows;}\n",
  "function cyFilter(id){var d=CY[id];var box=document.getElementById(id);",
  "var n=parseInt(box.querySelector('select.pagesel').value,10);",
  "var rows=cyRows(id);var total=rows.length;",
  "var show=(n===0)?rows:rows.slice(0,n);",
  "var h='<thead><tr>'+d.cols.map(function(c,i){",
  "var mark=(d.sort===i)?(d.dir===1?' \\u25B2':' \\u25BC'):'';",
  "return \"<th onclick='cySort(\\\"\"+id+\"\\\",\"+i+\")'>\"+cyEsc(c)+mark+'</th>';",
  "}).join('')+'</tr></thead><tbody>';",
  "h+=show.map(function(r){return '<tr>'+r.map(function(v){",
  "return '<td>'+cyEsc(v)+'</td>';}).join('')+'</tr>';}).join('');",
  "box.querySelector('table').innerHTML=h+'</tbody>';",
  "box.querySelector('p.shown').textContent='Showing '+show.length+' of '+total+",
  "(total===d.rows.length?'':' matching')+' row'+(total===1?'':'s')+",
  "'. Export CSV writes what is shown after search and sorting.';}\n",
  "function cySort(id,i){var d=CY[id];",
  "if(d.sort===i){d.dir=(d.dir===1)?-1:1;}else{d.sort=i;d.dir=1;}cyFilter(id);}\n",
  "function cyEsc(s){return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;')",
  ".replace(/>/g,'&gt;');}\n",
  # Export writes exactly what is on screen after search and sorting, which is
  # why it reads cyRows() rather than the untouched array.
  "function cyCsv(id){var d=CY[id];var box=document.getElementById(id);",
  "var n=parseInt(box.querySelector('select.pagesel').value,10);",
  "var rows=cyRows(id);if(n>0)rows=rows.slice(0,n);",
  "var q=function(v){v=String(v);return /[\",\\n]/.test(v)?'\"'+v.replace(/\"/g,'\"\"')+'\"':v;};",
  "var csv=[d.cols.map(q).join(',')].concat(rows.map(function(r){",
  "return r.map(q).join(',');})).join('\\n');",
  "var b=new Blob([csv],{type:'text/csv;charset=utf-8'});",
  "var a=document.createElement('a');a.href=URL.createObjectURL(b);",
  "a.download=box.querySelector('h3').textContent;document.body.appendChild(a);",
  "a.click();document.body.removeChild(a);URL.revokeObjectURL(a.href);}\n",
  # Zoom opens the embedded full-resolution image; there is no second copy and
  # no request, so it works from a file:// path and with no network.
  "var LBZ=1;\n",
  "function cyZoom(img){var lb=document.getElementById('lb');",
  "var i=document.getElementById('lbimg');i.src=img.src;LBZ=1;cyZset(1);",
  "document.getElementById('lbname').textContent=img.getAttribute('alt')||'';",
  "document.getElementById('lbdl').href=img.src;",
  "document.getElementById('lbdl').download=img.getAttribute('alt')||'figure.png';",
  "lb.classList.add('on');document.body.style.overflow='hidden';}\n",
  "function cyZset(z){LBZ=Math.min(8,Math.max(0.1,z));var i=document.getElementById('lbimg');",
  "i.style.width=(LBZ*100)+'%';document.getElementById('lbpct').textContent=",
  "Math.round(LBZ*100)+'%';}\n",
  "function cyClose(){document.getElementById('lb').classList.remove('on');",
  "document.body.style.overflow='';}\n",
  "document.addEventListener('keydown',function(e){",
  "if(e.key==='Escape')cyClose();",
  "if(document.getElementById('lb').classList.contains('on')){",
  "if(e.key==='+'||e.key==='=')cyZset(LBZ*1.25);",
  "if(e.key==='-')cyZset(LBZ/1.25);}});\n",
  # Expand all reaches the table toggles too -- they are part of the document now,
  # and a button that left half the content folded would not be expanding all.
  "function cyAll(open){document.querySelectorAll('details').forEach(function(d){",
  "d.open=open;});",
  "if(open)Object.keys(CY).forEach(cyShow);",
  "cySync();cyCur();}\n",
  # THE SIDEBAR MIRRORS THE DOCUMENT. A section's figure and table links are
  # shown exactly when that section is open, so Expand all and Collapse all move
  # both at once and the contents list never claims a state the page is not in.
  "function cySync(){document.querySelectorAll('section.sec').forEach(function(s){",
  "var d=s.querySelector('details');if(!d)return;",
  "var n=document.querySelector(\".nav-sec[data-sec='\"+s.id+\"']\");",
  "if(n)n.classList.toggle('open',d.open);});}\n",
  # A title in the contents TOGGLES its section: open it if closed, close it if
  # open. The default anchor jump is suppressed, because scrolling to a section
  # in the same gesture that collapses it puts the reader somewhere they did not
  # ask to be; opening still scrolls, which is what a contents entry is for.
  "document.addEventListener('click',function(e){",
  "var a=e.target&&e.target.closest?e.target.closest('a.nav-h'):null;if(!a)return;",
  "var h=a.getAttribute('href')||'';if(h.charAt(0)!=='#')return;",
  "var s=document.getElementById(h.slice(1));if(!s)return;",
  "var d=s.querySelector('details');if(!d)return;",
  "e.preventDefault();d.open=!d.open;",
  "if(d.open)s.scrollIntoView({block:'start'});",
  "cySync();cyCur();});\n",
  # WHICH SECTION AM I IN. The section whose top has most recently passed a line
  # near the top of the viewport, which is the one filling the screen. Chosen by
  # scan rather than IntersectionObserver because a collapsed section is only a
  # heading tall, so several intersect the viewport at once and "most visible"
  # picks the wrong one; "last one I have scrolled past" does not have that
  # failure mode and behaves the same whether sections are open or closed.
  "function cyCur(){var s=document.querySelectorAll('section.sec'),best=null;",
  "for(var i=0;i<s.length;i++){if(s[i].getBoundingClientRect().top<=140)best=s[i];}",
  "if(!best&&s.length)best=s[0];",
  "document.querySelectorAll('.nav-sec.cur').forEach(function(n){",
  "n.classList.remove('cur');});",
  "if(!best)return;",
  "var a=document.querySelector(\".nav-sec a.nav-h[href='#\"+best.id+\"']\");",
  "if(!a)return;a.parentNode.classList.add('cur');",
  # Keep the marked entry visible in a sidebar that lists every figure and table
  # and so runs well past one screen. The nav's own scrollTop is set directly
  # rather than calling scrollIntoView, which would also scroll the page and
  # fight the reader.
  "var nv=a.closest('nav.side');if(!nv)return;var at=a.offsetTop;",
  "if(at<nv.scrollTop||at+a.offsetHeight>nv.scrollTop+nv.clientHeight)",
  "nv.scrollTop=Math.max(0,at-nv.clientHeight/3);}\n",
  "window.addEventListener('scroll',cyCur,{passive:true});\n",
  "window.addEventListener('resize',cyCur,{passive:true});\n",
  # Opening or closing a section moves every section below it, so the marked
  # entry has to be recomputed then too.
  "document.addEventListener('toggle',function(e){",
  "if(e.target.tagName!=='DETAILS')return;",
  "if(e.target.open){var b=e.target.parentNode;",
  "if(b&&b.classList&&b.classList.contains('tab'))cyShow(b.id);}",
  "cySync();cyCur();},true);\n",
  "document.addEventListener('DOMContentLoaded',cyInit);\n")
}

#' Write a single self-contained HTML report of everything a run produced
#'
#' Presents the outputs in the order the documentation says they must be read,
#' because each stage can invalidate the ones after it. Sections whose files are
#' absent are omitted rather than shown empty.
#'
#' The file references nothing: every figure is embedded at full resolution as a
#' data URI and every table as JSON, so the report can be moved, emailed or
#' archived on its own. Figures can be zoomed and downloaded at full resolution;
#' tables can be searched, sorted, paged at 10/50/100/all rows and exported to
#' CSV.
#'
#' When `failure` is supplied the same report is written for a run that stopped
#' early, with a diagnosis section first and every output produced up to the
#' failure embedded below it. The partial output is usually where the evidence
#' is, so it is kept rather than discarded.
#'
#' @param outdir the results directory, which is also where the report is written
#' @param opt parsed options, used for the header
#' @param verdicts per-sample staining verdicts, used for the summary banner
#' @param failure a condition or message; when given, the report is written as
#'   the record of a failed run rather than a completed one.
#' @return the path, invisibly, or NULL when there is nothing to report
#' @export
write_run_report <- function(outdir, opt = NULL, verdicts = NULL,
                             failure = NULL) {
  if (!dir.exists(outdir)) return(invisible(NULL))
  path <- file.path(outdir, "report.html")
  failed <- !is.null(failure)

  # The banner. A run whose staining QC excluded samples says so before anything
  # else, because every frequency below is computed without them.
  banner <- ""
  if (failed) banner <- failure_block(failure, outdir, opt)
  else if (!is.null(verdicts) && length(verdicts)) {
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

  secs <- list(
    report_section(outdir, "s1", "Gate placement",
      paste("Per sample, the scatter gate boundary and every marker threshold,",
            "each drawn on the density it was derived from. A cut sitting in the",
            "trough between two modes is determined by the data. A cut on the",
            "flank of a single mode is a quantile fallback: not invalid, but",
            "carrying no evidence that the two populations separate. That",
            "distinction is visible here and in no downstream table, which is why",
            "this section comes first."),
      figures = c("recon_diagnostics.png", "gating_qc.png"), open = TRUE),

    report_section(outdir, "s2", "Acquisition stability",
      paste("The Time channel binned into equal-width intervals, tracking the",
            "event rate and each channel's median. A sustained trough is a partial",
            "clog, a spike is usually a bubble, a step is a settings change; any",
            "of them makes one file two instruments over its run, so a single",
            "threshold suits neither half. acquisition_qc.csv gives one verdict",
            "per file. acquisition_qc_impact.csv states how far each population",
            "would move if the flagged intervals were dropped, and that is the",
            "number the decision rests on: compare it against the same",
            "population's gate uncertainty in section 6. Nothing is removed",
            "unless --drop-unstable-events was given."),
      figures = "acquisition_qc.png",
      tables = c("acquisition_qc.csv", "acquisition_qc_impact.csv",
                 "acquisition_qc_bins.csv")),

    report_section(outdir, "s3", "Staining quality control",
      paste("One verdict per sample. A sample with no resolvable CD45+ mode has",
            "no usable parent gate, so its percentages are fractions of an",
            "arbitrary scatter region rather than of leukocytes. Such samples are",
            "excluded from every test below and the exclusion is recorded here,",
            "unless --include-qc-failed was given, in which case the verdict",
            "column still records it while qc_status reads pass."),
      tables = "staining_qc.csv"),

    report_section(outdir, "s4", "Phenotype concordance",
      paste("Measured marker intensity across the declared populations, z-scored",
            "over the run. Identity here is declared before the data are",
            "examined, so this figure serves the inverse function of its",
            "equivalent in clustering-first analysis: a population that does not",
            "express the markers its own definition requires falsifies the gate",
            "that produced it. The cohort heatmap shows each group's share of",
            "every population after normalising the groups to a common cell",
            "count, so an uneven split is not an artefact of unequal group",
            "sizes."),
      figures = c("population_marker_heatmap.png", "cohort_composition_heatmap.png")),

    report_section(outdir, "s5", "Threshold provenance and spillover spreading",
      paste("thresholds_used.csv records how every cut was obtained: a density",
            "minimum, a quantile fallback, an unstained or fluorescence-minus-one",
            "control, or a manual override. spreading_receivers.csv pairs each",
            "marker's fallback rate with how much wider its negative population",
            "becomes when another channel is bright. Reading the two columns",
            "together separates a cut that failed for an optical reason, which no",
            "gating strategy recovers and which a different fluorochrome",
            "assignment would fix, from one that failed for want of positive",
            "events or a titration problem."),
      tables = c("thresholds_used.csv", "threshold_scale_qc.csv",
                 "fmo_agreement.csv", "spreading_receivers.csv",
                 "spreading_pairs.csv")),

    report_section(outdir, "s6", "Gate uncertainty and detection limits",
      paste("Each frequency carries two separate uncertainties. The gate",
            "uncertainty is propagated from the thresholds behind the population",
            "and says how far the number moves when the cut moves;",
            "uncertainty_budget.csv attributes it to the individual markers, so",
            "the gate worth fixing is named. The counting uncertainty comes from",
            "the number of events observed. They are different guarantees: a cut",
            "through a wide empty gap is well determined however few events lie",
            "beyond it. detection_limits.png classifies every population against",
            "the limits of detection and quantification set by its own parent",
            "gate, and a population below them cannot be recovered by re-gating."),
      figures = c("frequency_uncertainty.png", "uncertainty_budget.png",
                  "detection_limits.png"),
      tables = c("uncertainty_budget.csv", "threshold_uncertainty.csv")),

    report_section(outdir, "s7", "Population abundance and the shared embedding",
      paste("pct_of_cd45_pos is the reportable quantity. count is an event count,",
            "set by how long the tube was run, and is not a cell number. The UMAP",
            "is computed once across all samples, so positions are comparable",
            "between them; per-sample embeddings would not be. The unsupervised",
            "clustering is computed without reference to the specification, which",
            "is what lets it contradict it: a cluster dominated by no declared",
            "label is a population the specification does not describe, and a",
            "declared label spread thinly across many clusters covers several",
            "distinct phenotypes."),
      figures = c("population_frequencies.png", "umap_overview.png",
                  "umap_overview_by_group.png", "umap_markers.png",
                  "umap_density.png", "umap_density_by_group.png",
                  "umap_multigraph_overlay.png", "unsupervised_clusters.png"),
      tables = c("population_frequencies.csv", "population_marker_mfi.csv",
                 "gate_counts.csv", "cluster_gate_agreement_populations.csv",
                 "cluster_gate_agreement_clusters.csv",
                 "subcluster_marker_shifts.csv", "unsupervised_clusters.csv",
                 "cells_umap.csv")),

    # The per-marker embeddings sit in their own directory because a run writes
    # one per marker and again per category, which is dozens of files. They are
    # the same shared embedding as section 7 with one marker's intensity on it,
    # so they belong beside it. Numbered 7b rather than renumbering everything
    # below, which would break every anchor an existing report already carries.
    report_section(outdir, "s7b", "The embedding, one marker at a time",
      paste("The same shared embedding as above, coloured by a single marker's",
            "intensity and drawn at full size. Each marker appears once pooled",
            "across every sample, then again split by each category the sample",
            "sheet carries, so a difference between groups can be read as a",
            "shift in where the signal sits rather than only as a number in a",
            "table. Read the pooled panel first: a split panel shows a subset of",
            "the same cells in the same positions, so anything that looks like a",
            "new structure there is a difference in density, not in phenotype."),
      figures = marker_umap_files(outdir)),

    report_section(outdir, "s8", "Between-group differences",
      paste("Tests are on per-sample values with donors as the replicates, not on",
            "pooled cells. Read each row in four steps: the adjusted p-value,",
            "then Cliff's delta as the effect size, then difference_over_gate_u,",
            "which expresses the difference in units of the gate's own",
            "uncertainty, then difference_over_total_u, which adds the counting",
            "uncertainty. A difference many times the gate uncertainty is not an",
            "artefact of threshold placement, but it is still not a finding until",
            "it survives correction across populations. Compositional results",
            "test the same abundances as centred log-ratios, because percentages",
            "of one parent cannot all move independently."),
      # group_differences.png comes FIRST: it is the whole comparison on one pair
      # of axes and it is how a reader picks which of the panels below to read
      # carefully. population_trajectories.png is last because it only exists for
      # a repeated-measures design and answers a different question -- direction
      # of travel rather than difference at a timepoint.
      figures = c("group_differences.png", "group_comparison.png",
                  "marker_state.png",
                  "functional_markers.png", "population_ratios.png",
                  "absolute_counts.png", "absolute_counts_qc.png",
                  "population_trajectories.png"),
      # design_feasibility comes FIRST because it decides which of the tests
      # below exist at all: a comparison it rules out is absent from the results
      # rather than negative, and a reader who meets that table after the
      # p-values has already drawn the wrong conclusion from their absence.
      # The four tables after it were reaching section 11, the catch-all, whose
      # own comment says a non-empty section 11 means a named section is missing
      # an entry. These are those entries.
      tables = c("design_feasibility.csv",
                 "group_comparison_stats.csv", "compositional_concordance.csv",
                 "compositional_clr_stats.csv", "marker_state_stats.csv",
                 "functional_markers.csv", "functional_markers_stats.csv",
                 "population_ratios.csv", "population_ratios_stats.csv",
                 "absolute_counts.csv", "absolute_counts_stats.csv",
                 "absolute_counts_raw.csv",
                 # These two were reaching the catch-all, whose whole purpose is
                 # to be empty. A paired test and a covariate-adjusted test are
                 # between-group results and belong beside the others.
                 "paired_comparison_stats.csv", "covariate_adjusted_stats.csv",
                 "subcluster_marker_shifts.csv",
                 "normality_tests.csv", "parametric_tests.csv",
                 "posthoc_tests.csv", "statistical_methods.csv")),

    report_section(outdir, "s8b", "Clinical variables",
      paste("A severity score, a laboratory value or an outcome flag against",
            "every population and every marker. Numeric variables are tested",
            "with Spearman's rho, two-level variables with Wilcoxon and Cliff's",
            "delta, more levels with Kruskal-Wallis, and p-values are adjusted",
            "within each variable rather than across all of them: each variable",
            "is its own question asked of every population.",
            "Read the effect before the asterisk. On a small cohort these tests",
            "detect only very large effects, so a null result says little and",
            "the `underpowered` column marks every test run on fewer than ten",
            "samples. This is association and not survival analysis: a 28-day",
            "flag is a two-group comparison, which is what it is, and no",
            "time-to-event model is fitted because the sheet carries no",
            "follow-up time.",
            "The figures are in reading order. The correlation between the",
            "variables comes first, because adjusting within each variable",
            "assumes they are separate questions and that figure is what shows",
            "whether they are: two variables strongly correlated with each other",
            "are one question asked twice. Then the association heatmap across",
            "every variable, then the cohort landscape, then per variable the",
            "effect sizes with their bootstrap intervals and the points behind",
            "them."),
      # Ordered explicitly rather than left to the glob: the glob would sort
      # clinical_association.png ahead of the correlogram the section tells the
      # reader to read first, and the per-variable files after both, which is the
      # order the description promises. The glob still runs last so a variable
      # named in the sheet cannot produce a figure no section names.
      # [A-Za-z0-9_], not [a-z0-9_]. The file name keeps the sheet column's case --
      # gsub("[^A-Za-z0-9]+", "_", cv) -- so a column written SOFA produces
      # clinical_SOFA.png, which a lower-case-only pattern missed. The figure then
      # reached the catch-all section instead of this one, which is the section
      # that explains how to read it.
      figures = unique(c("clinical_variables_correlation.png",
                  "clinical_association.png", "clinical_association_markers.png",
                  "clinical_landscape.png",
                  list.files(outdir, "^clinical_effects_[A-Za-z0-9_]+[.]png$"),
                  list.files(outdir, "^clinical_[A-Za-z0-9_]+[.]png$"))),
      tables = c("clinical_association.csv",
                 "clinical_association_markers.csv")),

    report_section(outdir, "s9", "Confounding and batch structure",
      paste("A variable confounds only when it both differs between the groups",
            "and associates with the outcome, and the two conditions are reported",
            "separately rather than combined into one verdict. threshold_drift",
            "asks whether the cuts themselves track the study group, which would",
            "make a difference definitional. batch_group_confounding.csv reports",
            "Cramer's V between batch and group: correction is refused above the",
            "configured threshold, because at that level of association removing",
            "the batch effect and removing the finding are the same operation."),
      # populations_by_batch asks a different question from batch_diagnostic:
      # whether the REPORTED NUMBERS step with batch, rather than whether the
      # embedding separates by it.
      figures = c("threshold_drift.png", "batch_diagnostic.png",
                  "populations_by_batch.png"),
      tables = c("threshold_drift_stats.csv", "confounding_diagnostics.csv",
                 "batch_group_confounding.csv", "batch_mixing_stats.csv",
                 "marker_batch_drift.csv")),

    report_section(outdir, "s10", "Agreement with an accepted baseline",
      paste("This run's per-sample values against a previously accepted run,",
            "written when --baseline was given. The within-run peer check finds",
            "one deviant tube against its peers. It cannot find a cohort that",
            "moved as a whole, because the peer median moves with it; only a",
            "baseline from a different run can."),
      tables = c("specification_conformance.csv", "specification_changes.csv")))

  # ---- completeness sweep ---------------------------------------------------
  # Every section above names its files, which is what puts them in reading
  # order. A file this run wrote but no section names would be invisible here
  # while appearing to be covered, so whatever is left over is collected rather
  # than dropped. If this section is ever non-empty for a standard run, a named
  # section is missing an entry.
  # RECURSIVE, deliberately. This section is the safety net for exactly the
  # failure described above, and a non-recursive listing defeated it for
  # anything written into a subdirectory. Two whole directories were therefore
  # absent from every report while appearing to be covered:
  # marker_umaps_by_group/, which is dozens of figures on every run, and
  # explore/. The first now has its own section; the second lands here.
  # Provenance, section 12, names patient_metadata_english.csv -- but it is
  # built AFTER this sweep, so its table was never in `used` and the catch-all
  # embedded the same file a second time, at full size, in every report. Name it
  # here rather than reordering the sections: the order below is the reading
  # order, and it should not bend around a bookkeeping detail.
  provenance_tables <- "patient_metadata_english.csv"
  used <- c(unlist(lapply(secs, `[[`, "used")), provenance_tables)
  rest_fig <- setdiff(list.files(outdir, "[.]png$", recursive = TRUE), used)
  rest_tab <- setdiff(list.files(outdir, "[.]csv$", recursive = TRUE), used)
  # An HTML output cannot be embedded as a figure or a table, so it would be
  # invisible here however the listing is done. Naming it is the most this
  # report can honestly do without ceasing to be self-contained.
  other_html <- setdiff(list.files(outdir, "[.]html$", recursive = TRUE),
                        "report.html")
  secs <- c(secs, list(report_section(outdir, "s11", "Further outputs",
    paste("Everything else this run wrote, in no particular order. These files",
          "are produced by optional flags or by stages the sections above do not",
          "cover; each is documented in the Output article."),
    figures = rest_fig, tables = rest_tab,
    body = if (length(other_html)) paste0(
      "<p>This run also wrote ",
      paste(sprintf("<span class='mono'>%s</span>", html_escape(other_html)),
            collapse = ", "),
      ". A report cannot embed another report without ceasing to be one file, ",
      "so it is named here and left in the results directory beside this one.</p>")
      else NULL)))

  # Provenance describes the other sections, so it is only worth writing when
  # there are some. A directory holding nothing produces no report rather than
  # a page whose sole content is a note about how to read the content.
  if (length(Filter(function(s) nzchar(s$html), secs)) || failed)
  secs <- c(secs, list(report_section(outdir, "s12", "Provenance",
    paste("What produced this folder, and what this file is."),
    tables = provenance_tables,
    body = paste0(
      "<p>This report is self-contained. Every figure above is embedded at full",
      " resolution and every table is embedded in full, so the file can be",
      " moved, attached to an email or archived on its own with nothing lost;",
      " it references no other file and needs no network. Click any figure to",
      " zoom, or use its download link for the original PNG. Any table can be",
      " searched, sorted by clicking a column heading, shown 10, 50, 100 or all",
      " rows at a time, and exported to CSV exactly as filtered and sorted.</p>",
      "<p>The complete run record stays in the results directory as ",
      "<span class='mono'>run_manifest.txt</span>",
      if (file.exists(file.path(outdir, "miflowcyt.md")))
        " and <span class='mono'>miflowcyt.md</span>" else "",
      ", which record the package versions, the invocation and every option in",
      " force.</p>"))))

  secs <- Filter(function(s) nzchar(s$html), secs)
  # A failed run that produced nothing still gets a report: the diagnosis is
  # the point of it, and an empty directory is exactly the case where the user
  # has least else to go on.
  if (!length(secs) && !failed) return(invisible(NULL))

  nfig <- sum(vapply(secs, `[[`, integer(1), "n_fig"))
  ntab <- sum(vapply(secs, `[[`, integer(1), "n_tab"))

  head <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>cyRAVEN run report</title>",
    sprintf("<style>%s</style>", report_css()), "</head><body>",
    "<div id='lb' onclick='if(event.target.id===\"lb\")cyClose()'>",
    "<div id='lbbar'><span id='lbname'></span>",
    "<button onclick='cyZset(LBZ/1.25)'>&minus;</button>",
    "<span id='lbpct'>100%</span>",
    "<button onclick='cyZset(LBZ*1.25)'>+</button>",
    "<button onclick='cyZset(1)'>Fit</button>",
    "<a class='dl' id='lbdl' download href='#'>Download</a>",
    "<button onclick='cyClose()'>Close</button></div>",
    "<img id='lbimg' alt='' onclick='event.stopPropagation()'/></div>",
    "<div class='wrap'>",
    "<nav class='side'><h2>Contents</h2>",
    paste(vapply(secs, `[[`, character(1), "nav"), collapse = "\n"),
    "</nav><main>",
    # The masthead is a flex row: mark, identity, actions. The two buttons used
    # to be a float:right span INSIDE the note below, which is why they hung off
    # its bottom-right corner and overlapped its border -- a float is taken out
    # of flow, so a note long enough to wrap pushed them past its own box. They
    # are controls for the whole document, not part of the note, and now sit
    # where that is legible.
    # One masthead block: the mark, and to its right the name, then what this
    # file is and when, then the note on how to read it. Everything identifying
    # the document sits in one column against the logo instead of trailing down
    # the page as three unrelated blocks.
    local({
      logo <- report_logo_uri()
      paste0(
        "<header class='masthead'>",
        if (!is.na(logo))
          sprintf("<img class='brand' src='%s' alt='cyRAVEN'/>", logo) else "",
        "<div class='mast-id'>",
        sprintf("<h1>cyRAVEN%s</h1>", if (failed) " &mdash; FAILED" else ""),
        sprintf(paste0("<p class='meta'>Run report &middot; %s &middot; ",
                       "cyRAVEN %s &middot; %d figures, %d tables, ",
                       "all embedded</p>"),
                html_escape(format(Sys.time(), tz = "UTC", usetz = TRUE)),
                html_escape(tryCatch(as.character(utils::packageVersion("cyRAVEN")),
                                     error = function(e) "unknown")),
                nfig, ntab),
        sprintf("<div class='banner warn mast-note'>%s</div>",
          if (failed)
            paste("What follows is everything the run wrote before it stopped,",
                  "in the order a completed run is read in. Stages after the",
                  "failure are absent, so a section missing here did not run",
                  "rather than finding nothing.")
          else
            paste("Read the sections in the order given. Each one can invalidate",
                  "the sections after it, so a result taken from the bottom",
                  "without the top is not supported by this run.")),
        "</div>",
        "<div class='mast-actions'>",
        "<button class='dl' onclick='cyAll(true)'>Expand all</button>",
        "<button class='dl' onclick='cyAll(false)'>Collapse all</button>",
        "</div></header>")
    }),
    banner)

  writeLines(c(head, vapply(secs, `[[`, character(1), "html"),
               "</main></div>", sprintf("<script>%s</script>", report_js()),
               "</body></html>"), path)
  sz <- file.size(path)
  log_msg("wrote report.html (", if (failed) "FAILED run, " else "",
          length(secs), " section(s), ", nfig,
          " figure(s) and ", ntab, " table(s) embedded, ",
          format(round(sz / 1024^2, 1), nsmall = 1), " MB, ",
          "self-contained: it references no other file)")
  invisible(path)
}
