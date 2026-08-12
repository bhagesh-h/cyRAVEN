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

  has <- function(...) {
    f <- c(...)
    f[file.exists(file.path(ex_dir, f))]
  }
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
    "1. What this run found",
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

  secs <- c(secs, list(report_section(ex_dir, "e2",
    "2. Quality control, decided from the data",
    paste("Explore mode embeds every event it is given, and on ungated whole",
          "blood most events are not leukocytes. The gate works on whole",
          "clusters rather than events: cluster coarsely, then judge each",
          "cluster by its marker profile. Read the `call` and `basis` columns",
          "before anything else -- a run whose gate went wrong produces",
          "internally consistent numbers that are unusable."),
    tables = has("explore_qc_clusters.csv"))))

  secs <- c(secs, list(report_section(ex_dir, "e3",
    "3. Abundance per donor",
    paste("One value per sample, never one per cell. Cluster frequencies carry",
          "the counting uncertainty of the events behind them and are",
          "classified against the limits of detection and quantification, on",
          "the same basis as the declared populations."),
    tables = has("explore_cluster_abundance.csv"))))

  secs <- c(secs, list(report_section(ex_dir, "e4",
    "4. Between-group differences",
    paste("Kruskal-Wallis across groups and Wilcoxon rank-sum against the",
          "reference, on per-donor frequencies, with Cliff's delta and",
          "Benjamini-Hochberg across every cluster tested. Where the run",
          "measured batch and group to be confounded, that verdict is carried",
          "in this table as its own columns: a q-value here means nothing if",
          "the groups were acquired on separate days."),
    tables = has("explore_cluster_stats.csv"),
    figures = has("explore_umap_by_group.png"))))

  if (length(bg))
    secs <- c(secs, list(report_section(ex_dir, "e5",
      "5. Each marker, split by group",
      paste("One panel per group per marker, on the shared embedding. This",
            "answers whether a marker is expressed in a DIFFERENT PLACE",
            "between groups, as opposed to more or less of it. Cells were",
            "equalised per sample before embedding, so panel density is",
            "comparable and not an artefact of group size."),
      figures = bg)))

  secs <- c(secs, list(report_section(ex_dir, "e6",
    "6. Against the declared specification",
    paste("Where the two views agree and where they do not. A cluster whose",
          "cells are mostly unlabelled is what the specification missed. A",
          "declared population spanning several clusters is what it lumped",
          "together -- the abundance of the whole may be flat while a subset",
          "inside it moves."),
    tables = has("explore_vs_populations.csv", "explore_population_split.csv"))))

  secs <- c(secs, list(report_section(ex_dir, "e7",
    "7. Marker expression over the embedding",
    paste("The sanity check that the map is organised by biology. Islands",
          "should light up for the markers that define them; if no marker",
          "structures the map, the clustering below it means little."),
    figures = has("explore_umap_markers.png"))))

  secs <- c(secs, list(report_section(ex_dir, "e8",
    "8. How this run was produced",
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

  secs <- Filter(function(s) nzchar(s$html), secs)
  if (!length(secs)) return(invisible(NULL))

  nfig <- sum(vapply(secs, `[[`, integer(1), "n_fig"))
  ntab <- sum(vapply(secs, `[[`, integer(1), "n_tab"))

  html <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>cyRAVEN explore report</title>",
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
    "<h1>cyRAVEN explore report</h1>",
    sprintf("<p class='q'>%s &middot; cyRAVEN %s &middot; %d figures, %d tables, all embedded</p>",
            html_escape(format(Sys.time(), tz = "UTC", usetz = TRUE)),
            html_escape(tryCatch(as.character(utils::packageVersion("cyRAVEN")),
                                 error = function(e) "unknown")),
            nfig, ntab),
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
    sprintf(paste0("<div class='banner warn'>Explore mode is a HYPOTHESIS",
                   " GENERATOR. A cluster is not a population until it has been",
                   " declared, thresholded per sample and checked.",
                   "<span style='float:right'>",
                   "<button class='dl' onclick='cyAll(true)'>Expand all</button> ",
                   "<button class='dl' onclick='cyAll(false)'>Collapse all</button>",
                   "</span></div>")),
    paste(vapply(secs, `[[`, character(1), "html"), collapse = "\n"),
    "</main></div>",
    sprintf("<script>%s</script>", report_js()),
    "</body></html>")

  writeLines(html, path, useBytes = TRUE)
  log_msg("wrote ", path, " (", nfig, " figures, ", ntab, " tables, self-contained)")
  invisible(path)
}
