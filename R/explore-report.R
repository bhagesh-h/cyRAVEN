# A separate self-contained report for explore mode.
#
# Separate rather than a section in report.html, for one reason: the isolation
# guarantee. --explore must not alter a single byte of the declared
# deliverables, and report.html is one of them. A section inside it would change
# that file on every run with the flag set, which is exactly the kind of quiet
# coupling this feature is supposed to avoid.
#
# So explore/explore_report.html stands on its own, built from the same
# machinery: every figure embedded at full resolution, every table embedded and
# searchable, nothing loaded from a network.

#' Write the explore-mode report
#'
#' @param ex_dir The explore directory, `<outdir>/explore`.
#' @param opt Option list, for the invocation line.
#' @return Path to the report, invisibly, or NULL when there is nothing to show.
#' @keywords internal
write_explore_report <- function(ex_dir, opt = NULL) {
  if (!dir.exists(ex_dir)) return(invisible(NULL))
  path <- file.path(ex_dir, "explore_report.html")
  # --report-images reaches this report too, so the two files a run writes
  # are either both self-contained or both light.
  .op <- options(cyRAVEN.report_images = opt$report_images %||% report_image_mode())
  on.exit(options(.op), add = TRUE)

  has <- function(...) {
    f <- c(...)
    f[file.exists(file.path(ex_dir, f))]
  }
  # Same lookup, but tolerant of the panel suffix. A run that resolves to more
  # than one panel writes explore_cluster_stats_panel_1.csv rather than
  # explore_cluster_stats.csv, so an exact-name lookup finds nothing and a
  # section gated on one disappears from the report entirely -- on exactly the
  # runs that have the most in them.
  has_stem <- function(...) {
    all_f <- list.files(ex_dir)
    out <- character(0)
    for (s in c(...)) {
      stem <- sub("[.][^.]+$", "", s)
      ext  <- sub("^.*[.]", "", s)
      out <- c(out, all_f[grepl(paste0("^", stem, "(_.*)?[.]", ext, "$"), all_f)])
    }
    unique(out)
  }
  # Exactly stem.ext or stem_panel_N.ext. has_stem() is looser -- any suffix --
  # which would let explore_cluster_mfi_heatmap.png also collect the by-group
  # variant that its own section names, embedding it twice.
  has_panels <- function(stem, ext = "png")
    sort(list.files(ex_dir, paste0("^", stem, "(_panel_[0-9]+)?[.]", ext, "$")))
  # Figures in the by-group subdirectory are referenced by relative path, which
  # report_section() resolves against ex_dir exactly like a top-level figure.
  bg <- list.files(file.path(ex_dir, "explore_marker_umaps_by_group"),
                   pattern = "\\.png$")
  bg <- if (length(bg)) file.path("explore_marker_umaps_by_group", bg) else character(0)

  prov <- NULL
  pp <- file.path(ex_dir, "explore_provenance.csv")
  if (file.exists(pp)) prov <- utils::read.csv(pp, stringsAsFactors = FALSE)
  basis <- if (!is.null(prov)) {
    b <- prov$value[prov$item == "positivity_basis"]
    if (length(b)) b[1] else NA_character_
  } else NA_character_
  thr_based <- identical(basis, "per-sample thresholds")

  secs <- list()

  secs <- c(secs, list(report_section(ex_dir, "e1",
    "What this run found",
    paste("Unsupervised clusters over every eligible channel. No population",
          "specification was used to produce them and the parent gate was not",
          "applied, so anything here that the declared analysis missed is",
          "genuinely a finding rather than a restatement.",
          if (thr_based)
            paste("Each cluster is named from the fraction of its cells above",
                  "THAT SAMPLE'S OWN threshold for each marker, so the",
                  "phenotype is a measurement rather than a reading off a",
                  "colour scale.")
          else
            paste("No declared pipeline ran, so cluster naming falls back to",
                  "pooled medians -- what a standalone clusterer has. Run with",
                  "a --config to name them against per-sample thresholds.")),
    figures = has("explore_umap_clusters.png", "explore_cluster_heatmap.png"),
    tables = has("explore_cluster_profile.csv", "explore_findings.csv"),
    open = TRUE)))

  secs <- c(secs, list(report_section(ex_dir, "e1b",
    "Composition of the embedding",
    paste("Every cluster in its own colour on the shared embedding, with a",
          "legend that says what each colour is: the cluster, the immune subset",
          "its marker profile matches where one does, its measured phenotype,",
          "and its share of the cells kept by the QC gate -- overall and within",
          "each study group. The phenotype is a MEASUREMENT, the markers this",
          "cluster's cells sit above or below their thresholds on, not a",
          "definition anyone wrote. The table adds each cluster's spread of",
          "shares across donors."),
    figures = has_panels("explore_umap_composition"),
    tables = has_panels("explore_cluster_composition", "csv"))))

  secs <- c(secs, list(report_section(ex_dir, "e2",
    "Quality control, decided from the data",
    paste("Explore mode draws single cells from inside the scatter and",
          "singlet gates and applies no CD45 or viability gate of its own, so",
          "dead cells and non-leukocytes are still among them. The gate works",
          "on whole clusters rather than events: cluster coarsely, then judge",
          "each cluster by its marker profile. Read the `call` and `basis`",
          "columns before anything else -- a run whose gate went wrong produces",
          "internally consistent numbers that are unusable."),
    tables = has("explore_qc_clusters.csv"))))

  if (length(has_stem("explore_cluster_median_heatmap.png")))
    secs <- c(secs, list(report_section(ex_dir, "e2b",
      "The same clusters as scaled median expression",
      paste("The heatmap above asks what share of each cluster is above its",
            "own sample's threshold. This one asks how bright each cluster is",
            "for each marker relative to the other clusters, which needs no",
            "threshold at all. That matters when the thresholds are weak: a",
            "cluster can read uniformly negative on the first figure and still",
            "be clearly separated on this one, and the difference is a",
            "statement about the gate rather than about the cells. It is also",
            "the view an unsupervised tool with no gating step shows, so it is",
            "the figure to compare against FlowSOM star plots or cyCONDOR's",
            "cluster_marker_heatmap.png. Rows and columns are ordered by",
            "hierarchical clustering, so adjacency means resemblance."),
      figures = has_stem("explore_cluster_median_heatmap.png"))))

  secs <- c(secs, list(report_section(ex_dir, "e2c",
    "Marker intensity per cluster (MFI)",
    paste("Median fluorescence intensity per cluster and marker in instrument",
          "units: each sample's median within the cluster first, then the median",
          "over samples, so the donor with the most events does not set it. The",
          "heatmap prints the value in every tile and colours it on the analysis",
          "scale rescaled per marker, rows ordered by similarity, with each",
          "cluster's share of the cells beside it. In the dot plot, size is the",
          "percent positive and colour the scaled intensity. Percent positive is",
          "against each sample's own threshold where the declared run supplied",
          "one, and otherwise the cluster's pooled call."),
    figures = c(has_panels("explore_cluster_mfi_heatmap"),
                has_panels("explore_cluster_mfi_dotplot"),
                has_panels("explore_cluster_mfi_heatmap_by_group")),
    tables = has_panels("explore_cluster_mfi", "csv"))))

  secs <- c(secs, list(report_section(ex_dir, "e3",
    "Abundance per donor",
    paste("One value per sample, never one per cell. Cluster frequencies carry",
          "the counting uncertainty of the events behind them and are",
          "classified against the limits of detection and quantification, on",
          "the same basis as the declared populations."),
    tables = has("explore_cluster_abundance.csv"))))

  secs <- c(secs, list(report_section(ex_dir, "e4",
    "Between-group differences",
    paste("Kruskal-Wallis across groups and Wilcoxon rank-sum against the",
          "reference, on per-donor frequencies, with Cliff's delta and",
          "Benjamini-Hochberg across every cluster tested. Where the run",
          "measured batch and group to be confounded, that verdict is carried",
          "in this table as its own columns: a q-value here means nothing if",
          "the groups were acquired on separate days."),
    tables = has("explore_cluster_stats.csv"),
    figures = has("explore_umap_by_group.png"))))

  # Written only when --total-counts supplied one, so a run without it produces
  # the same report it always did.
  if (length(has_stem("explore_cluster_stats_absolute.csv",
                      "explore_absolute_vs_share.png",
                      "explore_total_counts_qc.png")))
    secs <- c(secs, list(report_section(ex_dir, "e4b",
      "Absolute cell numbers",
      paste("Everything above this section is compositional: a share of the",
            "events acquired, and shares are constrained to sum to 100. That",
            "constraint means no frequency table can separate one cluster",
            "expanding from every other cluster contracting -- the composition",
            "is identical either way, and the centred log-ratio does not fix",
            "it either. Multiplying each share by that acquisition's own",
            "externally measured total lifts the constraint, because cell",
            "numbers are free to all move the same way.",
            "READ THE QC FIGURE FIRST: everything here inherits the errors of",
            "those external totals, and a yield entered in the wrong unit is",
            "obvious there and invisible in the derived table.",
            "The concordance table is the point of the pair -- a cluster",
            "significant on cell number but not on share is a compartment that",
            "changed size; the reverse is a redistribution at constant size.",
            "These are DUAL-PLATFORM numbers, one measurement from this run",
            "multiplied by one from an instrument it never saw, and the",
            "published interlaboratory CVs for that route are roughly 20-33%",
            "against 10-16% for single-platform bead counting. Treat a",
            "difference smaller than that as unresolved."),
      tables = c(has_stem("explore_cluster_stats_absolute.csv"),
                 has_stem("explore_cluster_count_concordance.csv")),
      figures = c(has_stem("explore_total_counts_qc.png"),
                  has_stem("explore_absolute_vs_share.png")))))

  if (length(bg))
    secs <- c(secs, list(report_section(ex_dir, "e5",
      "Each marker, split by group",
      paste("One panel per group per marker, on the shared embedding. This",
            "answers whether a marker is expressed in a DIFFERENT PLACE",
            "between groups, as opposed to more or less of it. Cells were",
            "equalised per sample before embedding, so panel density is",
            "comparable and not an artefact of group size."),
      figures = bg)))

  secs <- c(secs, list(report_section(ex_dir, "e6",
    "Against the declared specification",
    paste("Where the two views agree and where they do not. The identity",
          "figure is the one to read first: rows are clusters, grouped by the",
          "declared population each best matches, so it answers 'this cluster",
          "is what?' directly. The match is scored by F1 -- high only when the",
          "cluster is mostly that population AND that population is mostly in",
          "this cluster -- because the largest-overlap label alone names a",
          "cluster after a population that barely overlaps it. A cluster",
          "mostly made of the catch-all is what the specification MISSED, and",
          "is shown as undescribed rather than being given a name. A declared",
          "population spanning several clusters is what it lumped together --",
          "the abundance of the whole may be flat while a subset inside it",
          "moves."),
    figures = has_stem("explore_cluster_identity.png"),
    # explore_findings.csv is NOT named here: section 2 already names it, and a
    # file named by two sections is embedded twice, at full size, in a report
    # that is already the largest artefact the run writes.
    tables = c(has_stem("explore_cluster_identity.csv"),
               has("explore_vs_populations.csv",
                   "explore_population_split.csv")))))

  secs <- c(secs, list(report_section(ex_dir, "e7",
    "Marker expression over the embedding",
    paste("The sanity check that the map is organised by biology. Islands",
          "should light up for the markers that define them; if no marker",
          "structures the map, the clustering below it means little."),
    figures = has("explore_umap_markers.png"))))

  secs <- c(secs, list(report_section(ex_dir, "e8",
    "How this run was produced",
    paste("Every choice explore mode made, and on what basis. Where the",
          "declared pipeline supplied something -- the estimated transform,",
          "per-sample thresholds, the confounding verdict -- this table says",
          "so; where it fell back to what a standalone clusterer would do, it",
          "says that instead."),
    tables = has("explore_provenance.csv"),
    body = paste0(
      "<p>A draft population specification is written alongside this report as ",
      "<span class='mono'>explore_suggested_spec.yaml</span>. It is a starting ",
      "point, not a result: merge the clusters that are one population, drop ",
      "the debris and the doublets, give them real names, then run the ",
      "supervised path with <span class='mono'>--config</span>, which is where ",
      "they acquire per-sample thresholds, propagated uncertainty and the six ",
      "specification checks.</p>",
      "<p>This report is self-contained. It references no other file and needs ",
      "no network.</p>"))))

  secs <- report_number(secs)
  if (!length(secs)) return(invisible(NULL))

  nfig <- sum(vapply(secs, `[[`, integer(1), "n_fig"))
  ntab <- sum(vapply(secs, `[[`, integer(1), "n_tab"))

  masthead <- paste0(
    "<h1>cyRAVEN explore report</h1>",
    sprintf("<p class='q'>%s &middot; cyRAVEN %s &middot; %d sections, %d figures, %d tables, %s</p>",
            html_escape(format(Sys.time(), tz = "UTC", usetz = TRUE)),
            html_escape(tryCatch(as.character(utils::packageVersion("cyRAVEN")),
                                 error = function(e) "unknown")),
            length(secs), nfig, ntab,
            if (identical(report_image_mode(), "link"))
              "figures linked from this directory" else "all embedded"))
  banner <- paste0(
    sprintf("<div class='banner %s'>%s</div>",
            if (thr_based) "ok" else "warn",
            if (thr_based)
              paste("Unsupervised discovery, with clusters named against each",
                    "sample's own marker thresholds. Nothing here was used to",
                    "produce the declared results in the parent directory.")
            else
              paste("Standalone unsupervised discovery. No declared pipeline",
                    "ran, so clusters are named from pooled medians rather than",
                    "per-sample thresholds.")),
    "<div class='banner warn'>Explore mode is a HYPOTHESIS GENERATOR. A cluster",
    " is not a population until it has been declared, thresholded per sample",
    " and checked.</div>")
  html <- report_page(secs, "cyRAVEN explore report", masthead, banner)

  writeLines(html, path, useBytes = TRUE)
  log_msg("wrote ", path, " (", length(secs), " numbered sections, ", nfig,
          " figures, ", ntab, " tables, ",
          if (identical(report_image_mode(), "link"))
            "figures linked from this directory" else "self-contained", ")")
  invisible(path)
}
