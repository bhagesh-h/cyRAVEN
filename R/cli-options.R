build_option_list <- function() list(
  optparse::make_option("--files", type = "character", default = NULL,
              help = "FCS files: comma-separated paths or a glob in quotes"),
  optparse::make_option("--dir", type = "character", default = NULL,
              help = "directory to search for .fcs files (alternative to --files)"),
  optparse::make_option("--recursive", action = "store_true", default = FALSE,
              help = "search --dir recursively (cohorts in subdirectories)"),
  optparse::make_option("--pattern", type = "character", default = NULL,
              help = paste("regex for filenames under --dir [default '[.]fcs$'].",
                           "Use '[.]fcs( copy)?$' for duplicated archives")),
  optparse::make_option("--exclude", type = "character", default = NULL,
              help = paste("regex; matching paths are dropped after discovery.",
                           "Default drops single-stain compensation controls,",
                           "which are instrument setup files and not samples.",
                           "Pass '' to keep everything")),
  optparse::make_option("--max-events-per-file", type = "integer", default = 0L,
              dest = "max_events_per_file",
              help = paste("read at most N events per file, sampled evenly through",
                           "the acquisition (0 = all). Bounds peak memory on large",
                           "cohorts; gate thresholds are derived from the subsample")),
  optparse::make_option("--outdir", type = "character", default = "results",
              help = "output directory [%default]"),
  optparse::make_option("--sample-map", type = "character", default = NULL,
              dest = "sample_map", help = "CSV mapping file -> sample/patient (see README)"),
  optparse::make_option("--patient-table", type = "character", default = NULL,
              dest = "patient_table", help = "patient metadata CSV (German or English)"),
  optparse::make_option("--absolute-counts", type = "character", default = NULL,
              dest = "absolute_counts",
              help = paste("externally measured absolute cell counts: .xlsx, .csv or",
                           ".tsv, wide format (row 1 = population names, column 1 =",
                           "sample id/filename, matched against --sample-map). Requires",
                           "--sample-map. Produces absolute_counts_qc.png (inspect",
                           "first), absolute_counts.png/.csv and, with --group-column,",
                           "absolute_counts_stats.csv. See system.file(",
                           "'examples', package = 'cyRAVEN')")),
  optparse::make_option("--reference-date", type = "character", default = NULL,
              dest = "reference_date",
              help = paste("YYYY-MM-DD used to derive age from date of birth and",
                           "to place 2-digit years [default: today]. Set this to a",
                           "fixed study date to make ages reproducible across runs.")),
  optparse::make_option("--group-column", type = "character", default = NULL,
              dest = "group_column",
              help = paste("patient-table column defining the comparison groups",
                           "[default: the resolved 'cohort' column, if present].",
                           "Enables the between-group abundance figure and tests.")),
  optparse::make_option("--reference-group", type = "character", default = NULL,
              dest = "reference_group",
              help = paste("group every other group is tested against, e.g.",
                           "'controls' [default: first alphabetically].",
                           "Drawn unfilled and leftmost in each panel.")),
  optparse::make_option("--p-adjust-display", type = "character", default = "raw",
              dest = "p_adjust_display",
              help = paste("which p-value the figure brackets show: 'raw' or 'BH'",
                           "[default raw]. Both are always written to the CSV.")),
  optparse::make_option("--config", type = "character", default = NULL,
              help = "YAML config of threshold/population overrides"),
  optparse::make_option("--write-config", type = "character", default = NULL,
              dest = "write_config",
              help = "derive everything, write this YAML config, and exit"),
  optparse::make_option("--write-sample-map", type = "character", default = NULL,
              dest = "write_sample_map", help = "write a sample-map template and exit"),
  optparse::make_option("--cells-per-sample", type = "integer", default = 20000L,
              dest = "cells_per_sample", help = "cells drawn per sample [%default]"),
  optparse::make_option("--max-cells", type = "integer", default = 200000L,
              dest = "max_cells", help = "overall cell ceiling per embedding [%default]"),
  optparse::make_option("--cofactor", type = "double", default = NULL,
              help = "fixed asinh cofactor (default: derive per panel)"),
  optparse::make_option("--n-neighbors", type = "integer", default = 30L,
              dest = "n_neighbors", help = "UMAP n_neighbors [%default]"),
  optparse::make_option("--min-dist", type = "double", default = 0.3,
              dest = "min_dist", help = "UMAP min_dist [%default]"),
  optparse::make_option("--singlet-mad-k", type = "double", default = 3,
              dest = "singlet_mad_k", help = "singlet band width in MADs [%default]"),
  optparse::make_option("--min-cd45-pct", type = "double", default = 5,
              dest = "min_cd45_pct", help = "staining QC floor, %% CD45+ of live [%default]"),
  optparse::make_option("--include-qc-failed", action = "store_true", default = FALSE,
              dest = "include_qc_failed",
              help = paste("include declared samples (not controls) that fail staining",
                           "QC (no separable CD45+ mode, or below --min-cd45-pct) instead",
                           "of excluding them. Their pct_of_cd45_pos carries no real",
                           "evidence, see staining_qc.csv's verdict for which samples",
                           "were forced in, and treat their numbers with caution")),
  optparse::make_option("--viability-marker", type = "character", default = NULL,
              dest = "viability_marker", help = "viability dye name (default: auto-detect)"),
  optparse::make_option("--umap-markers", type = "character", default = NULL,
              dest = "umap_markers",
              help = paste("comma-separated lineage markers preferred as UMAP features",
                           "(default: a built-in panel-agnostic lineage list; markers",
                           "absent from a given panel are skipped). Override this for",
                           "a panel where that default list has fewer than 2 matches",
                           "instead of relying on the automatic all-eligible-markers",
                           "fallback")),
  optparse::make_option("--umap-markers-all", action = "store_true", default = FALSE,
              dest = "umap_markers_all",
              help = "use every eligible marker as a UMAP feature instead of preferring lineage markers"),
  optparse::make_option("--threads", type = "integer", default = 0L,
              help = "UMAP threads, 0 = all cores [%default]"),
  optparse::make_option("--seed", type = "integer", default = 42L, help = "RNG seed [%default]"),
  optparse::make_option("--keep-exprs", action = "store_true", default = FALSE,
              dest = "keep_exprs",
              help = "include raw expression matrices in the saved session (large)"),
  optparse::make_option("--no-session", action = "store_true", default = FALSE,
              dest = "no_session", help = "skip writing session_state.RData"),

  # ---- FlowJo hand-off (additive; changes nothing about the default run) ----
  # WHY a flag rather than always-on: the export writes one FCS per sample plus
  # a concatenated copy of every cell, which roughly doubles the run's output
  # footprint. Runs that only want the figures should not pay that.
  optparse::make_option("--flowjo-export", action = "store_true", default = FALSE,
              dest = "flowjo_export",
              help = paste("ALSO write UMAP-annotated FCS files for interactive",
                           "inspection in the FlowJo GUI. Every default output is",
                           "still produced; this adds to them. FlowJo's command line",
                           "cannot compute a UMAP (it is a plugin), so the embedding",
                           "is carried in as extra FCS parameters, see",
                           "docs/flowjo_interactive.md")),
  optparse::make_option("--flowjo-outdir", type = "character", default = NULL,
              dest = "flowjo_outdir",
              help = "directory for --flowjo-export output [default <outdir>/flowjo]"),
  optparse::make_option("--flowjo-no-concat", action = "store_true", default = FALSE,
              dest = "flowjo_no_concat",
              help = paste("with --flowjo-export, skip the concatenated",
                           "_ALL_SAMPLES.fcs and write per-sample files only")),
  optparse::make_option("--flowjo-no-groups", action = "store_true", default = FALSE,
              dest = "flowjo_no_groups",
              help = paste("with --flowjo-export, skip the per-group",
                           "_GROUP_<cohort>.fcs files")),
  # WHY "Other CD45+" IS HIDDEN BY DEFAULT: it is not a population, it is the
  # leftover -- every CD45+ cell that matched no definition in the spec. Here it
  # is the single largest label (36% of cells), so on a UMAP it blankets the
  # real populations and dominates every legend while meaning nothing
  # biologically. Hiding it is a DISPLAY choice only: those cells are still
  # gated, still counted, still in cells_umap.csv, gate_counts.csv and every
  # frequency table. --other puts them back on the figures.
  optparse::make_option("--other", action = "store_true", default = FALSE,
              dest = "include_other",
              help = paste("also draw the 'Other CD45+' catch-all on the UMAP",
                           "figures [hidden by default; never affects tables]")),

  # =========================================================================
  # EXTENDED ANALYSES
  # =========================================================================
  # DESIGN RULE FOR THESE FLAGS: anything that only ADDS an output is ON by
  # default, because an analysis you have to know to ask for is one nobody runs.
  # Anything that CHANGES an existing number is OFF by default or has an escape
  # hatch, because a silent change to a published figure is the worst outcome
  # this pipeline can produce.

  optparse::make_option("--no-differential-state", action = "store_true", default = FALSE,
              dest = "no_differential_state",
              help = paste("skip the sample-level differential-state tests",
                           "(marker_state_stats.csv + marker_state.png)")),
  optparse::make_option("--no-heatmaps", action = "store_true", default = FALSE,
              dest = "no_heatmaps",
              help = paste("skip the population phenotype heatmap and the cohort",
                           "composition heatmap")),
  optparse::make_option("--no-compositional", action = "store_true", default = FALSE,
              dest = "no_compositional",
              help = paste("skip the CLR (compositional) re-test of population",
                           "frequencies and its concordance table")),
  optparse::make_option("--no-confounding", action = "store_true", default = FALSE,
              dest = "no_confounding",
              help = paste("skip the age/sex confounding diagnostic")),
  optparse::make_option("--covariates", type = "character", default = "age,sex",
              dest = "covariates",
              help = paste("comma-separated patient-table columns to screen as",
                           "confounders [default %default]")),
  optparse::make_option("--rank-ancova", action = "store_true", default = FALSE,
              dest = "rank_ancova",
              help = paste("additionally fit an EXPLORATORY covariate-adjusted",
                           "rank ANCOVA. Off by default: at single-digit n per",
                           "cohort the adjustment is usually extrapolation,",
                           "read confounding_diagnostics.csv first")),
  optparse::make_option("--paired-column", type = "character", default = NULL,
              dest = "paired_column",
              help = paste("patient-table or sample-map column identifying the",
                           "PAIRING UNIT (e.g. donor) for a paired design.",
                           "Requires --condition-column. Pairing cannot be",
                           "inferred, so nothing paired runs without this")),
  optparse::make_option("--condition-column", type = "character", default = NULL,
              dest = "condition_column",
              help = paste("column naming the condition within a pair",
                           "(e.g. timepoint), used with --paired-column")),
  optparse::make_option("--batch-column", type = "character", default = NULL,
              dest = "batch_column",
              help = paste("column identifying the acquisition batch (run date,",
                           "operator, instrument). Enables the batch-mixing",
                           "diagnostic. Measuring the effect is separate from",
                           "removing it: add --correct-batch to correct, which",
                           "is refused when batch and group overlap too far")),
  optparse::make_option("--no-threshold-drift", action = "store_true", default = FALSE,
              dest = "no_threshold_drift",
              help = paste("skip the check for gating thresholds that differ",
                           "systematically between cohorts")),
  optparse::make_option("--cluster", action = "store_true", default = FALSE,
              dest = "unsupervised",
              help = paste("also run unsupervised clustering (SOM + consensus",
                           "metaclustering) and cross-check it against the gate",
                           "spec. Off by default because it adds runtime; it is",
                           "the only output that can find a population the",
                           "config does not describe")),
  optparse::make_option("--cluster-k", type = "integer", default = 12L, dest = "cluster_k",
              help = "metaclusters for --cluster [default %default]"),
  optparse::make_option("--cluster-grid", type = "integer", default = 10L,
              dest = "cluster_grid",
              help = "SOM grid side for --cluster; grid^2 nodes [default %default]"),
  optparse::make_option("--transform", type = "character", default = "arcsinh",
              help = paste("intensity transform: 'arcsinh' (default, cofactor",
                           "derived from the data), 'logicle' (the flow",
                           "cytometry convention -- linear near zero so",
                           "compensated negatives stay on scale, parameters",
                           "pooled per panel), or 'none' [default %default]")),
  optparse::make_option("--logicle-m", type = "double", default = 4.5,
              dest = "logicle_m",
              help = "decades on the logicle display scale [default %default]"),
  optparse::make_option("--correct-batch", action = "store_true", default = FALSE,
              dest = "correct_batch",
              help = paste("align each marker's distribution across batches",
                           "before embedding. Requires --batch-column. REFUSES",
                           "when batch is confounded with the study group,",
                           "because there removing the batch and removing the",
                           "finding are the same operation")),
  optparse::make_option("--force-batch-correction", action = "store_true",
              default = FALSE, dest = "force_batch_correction",
              help = paste("correct even when batch and group are confounded.",
                           "The decision is recorded in the run manifest")),
  optparse::make_option("--batch-max-cramers-v", type = "double", default = 0.6,
              dest = "batch_max_cramers_v",
              help = paste("refuse batch correction above this Cramer's V",
                           "between batch and group [default %default]")),
  optparse::make_option("--explain-clusters", action = "store_true", default = FALSE,
              dest = "explain_clusters",
              help = paste("for each cluster the population spec does not",
                           "describe, LEARN a two-marker gating strategy that",
                           "selects it, and report precision/recall/F1 on cells",
                           "held out of the fit. Turns \"cluster 7 matches",
                           "nothing\" into a gate a cytometrist can draw.",
                           "Requires --cluster. Descriptive only: proposes gates,",
                           "never changes the scored populations")),
  optparse::make_option("--explain-max-clusters", type = "integer", default = 4L,
              dest = "explain_max_clusters",
              help = paste("ceiling on how many gating strategies to derive;",
                           "the largest qualifying clusters are taken and the",
                           "number skipped is logged [default %default]")),
  optparse::make_option("--explain-max-depth", type = "integer", default = 4L,
              dest = "explain_max_depth",
              help = paste("maximum gates in a learned strategy. More levels",
                           "raise precision and cost recall; the search stops",
                           "early when another gate does not earn its place",
                           "[default %default]")),
  optparse::make_option("--auto-subcluster-k", action = "store_true", default = FALSE,
              dest = "auto_subcluster_k",
              help = paste("choose the subcluster count per population by mean",
                           "silhouette instead of the fixed k. CHANGES the",
                           "subcluster lettering in the multigraph overlay, so",
                           "it is opt-in")),
  optparse::make_option("--save-umap-model", type = "character", default = NULL,
              dest = "save_umap_model",
              help = paste("write the trained UMAP model to this path so a later",
                           "batch can be PROJECTED into the same embedding",
                           "instead of moving every coordinate")),
  optparse::make_option("--umap-model", type = "character", default = NULL,
              dest = "umap_model",
              help = paste("project this run's cells into a previously saved",
                           "model instead of embedding fresh. Right for adding",
                           "more of the same kind of sample; wrong for a new",
                           "panel or a new cell type, retrain instead")),
  # DEST DELIBERATELY DOES NOT START WITH "cofactor".
  #
  # `$` on a list partial-matches, and optparse OMITS options whose default is
  # NULL from the list it returns -- so with a dest of "cofactor_from_first_sample"
  # the existing `opt$cofactor` (default NULL, therefore absent) silently
  # partial-matched THIS flag and returned FALSE. FALSE is not NULL, so the
  # pipeline took it as a user-supplied cofactor of 0, every marker became
  # asinh(x/0) = Inf, every threshold came back Inf, and all 25 samples failed
  # staining QC with a message about "no separable CD45+ mode" that pointed
  # nowhere near the cause. The baseline hit the same trap once before -- see the
  # opt[["umap_markers", exact = TRUE]] comment in STEP 6 -- and the call site
  # below now uses exact indexing as well. Two guards, because this failure is
  # silent, total, and blames the data.
  optparse::make_option("--cofactor-from-first-sample", action = "store_true", default = FALSE,
              dest = "first_sample_cofactor",
              help = paste("derive the arcsinh cofactor from the first file only,",
                           "as before. The default now takes the median of up to",
                           "8 per-sample derivations, which changes derived",
                           "numbers when samples disagree"))
)

#' Directory this package is installed in
#'
#' Recorded in the run manifest so a results folder can be traced back to the
#' code that produced it. As a script this resolved the `--file=` argument;
#' in a package the installed location is both simpler and more truthful.
#' @keywords internal
script_dir <- function() {
  p <- tryCatch(find.package("cyRAVEN"), error = function(e) NULL)
  p %||% getwd()
}

#' Resolve the input file list, and REPORT what was skipped
#'
#' WHY the near-miss warning: real archives contain files whose names end in
#' something other than a clean ".fcs" -- "sample.fcs copy" from a filesystem
#' duplication, "sample.FCS", "sample.fcs.bak". A strict pattern silently drops
#' them, and silently analysing 19 of 25 files is far worse than failing: the
#' output looks complete and nothing in it reveals the missing cohort. So any
#' file whose name CONTAINS ".fcs" but did not match the pattern is listed as a
#' warning with the exact flag needed to include it.
#' @param opt Named list of options. See build_option_list() for the full set.
#' @keywords internal
resolve_input_files <- function(opt) {
  pat <- opt$pattern %||% "[.]fcs$"
  fs <- character(0)
  if (!is.null(opt$files)) {
    parts <- trimws(strsplit(opt$files, ",")[[1]])
    for (p in parts) fs <- c(fs, if (grepl("[*?]", p)) Sys.glob(p) else p)
  }
  if (!is.null(opt$dir)) {
    rec <- isTRUE(opt$recursive)
    fs <- c(fs, list.files(opt$dir, pattern = pat, ignore.case = TRUE,
                           full.names = TRUE, recursive = rec))
    # Anything that looks like an FCS file but did not match the pattern.
    all_f <- list.files(opt$dir, full.names = TRUE, recursive = rec)
    near  <- setdiff(all_f[grepl("[.]fcs", all_f, ignore.case = TRUE)], fs)
    near  <- near[!dir.exists(near)]
    if (length(near)) {
      log_msg("WARNING: ", length(near), " file(s) contain '.fcs' but do not match ",
              "--pattern '", pat, "' and were SKIPPED:")
      for (n in utils::head(near, 10)) log_msg("    ", basename(n))
      if (length(near) > 10) log_msg("    ... and ", length(near) - 10, " more")
      log_msg("  to include them, pass --pattern '[.]fcs( copy)?$' (or a pattern of your own)")
    }
    if (!rec && !length(fs) && length(list.dirs(opt$dir, recursive = FALSE)))
      log_msg("NOTE: no files matched at the top level of --dir, but it contains ",
              "subdirectories, pass --recursive to search them")
  }
  fs <- unique(fs[file.exists(fs)])

  # Drop instrument setup acquisitions.
  #
  # WHY a default rather than leaving it to the user: single-stain compensation
  # controls sit in the same archive as the samples and are FCS files like any
  # other, but each contains beads or cells carrying ONE fluorophore. They have no
  # CD45 population, no viable-cell population, and no marker set in common with
  # the panel -- so they form their own spurious "panel" group, fail staining QC,
  # and clutter every figure. They are also the reason a permissive --pattern is
  # dangerous. The default is a name-based heuristic, so it is announced in the
  # log and fully overridable; passing --exclude '' disables it.
  ex_pat <- if (is.null(opt$exclude)) "compensation|single[ _-]?stain|comp[ _-]?control"
            else opt$exclude
  if (nzchar(ex_pat)) {
    drop <- grepl(ex_pat, fs, ignore.case = TRUE)
    if (any(drop)) {
      log_msg("excluded ", sum(drop), " file(s) matching --exclude '", ex_pat, "':")
      for (n in utils::head(basename(fs[drop]), 5)) log_msg("    ", n)
      if (sum(drop) > 5) log_msg("    ... and ", sum(drop) - 5, " more")
      fs <- fs[!drop]
    }
  }
  if (!length(fs)) stop("no FCS files found, pass --files or --dir")
  fs
}
