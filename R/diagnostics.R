# =============================================================================
# cyRAVEN -- run diagnostics: what could make this result set wrong?
#   * batch / acquisition-date effect: diagnosed, not silently corrected
#   * threshold drift: does the GATE move with cohort?
#   * run manifest: what produced this results folder
#
# Every function here is read-only with respect to the analysis: it reports on a
# result set without altering it, so a diagnostic can never quietly change the
# numbers it is diagnosing.
# =============================================================================


# =============================================================================
# BATCH EFFECT: diagnosed, deliberately not corrected
# =============================================================================
#
# THE GAP THIS CLOSES. The comparable toolboxes in this field ship batch
# CORRECTORS -- Harmony on the intensities or the PCA, and CytoNorm's
# FlowSOM-based quantile splines against reference samples. This package has
# neither, and -- more to the point -- until this function it had no way to tell
# whether it needed one. An island on the UMAP that is an acquisition batch is
# drawn identically to an island that is a phenotype.
#
# WHY THIS IMPLEMENTS THE DIAGNOSTIC AND NOT THE CORRECTION. Correction is not
# free and is not neutral. Harmony moves cells so that batch labels mix; when the
# batch is confounded with the biological group -- samples of one cohort acquired
# on one set of days, which is the normal way a clinical cohort is collected --
# "removing the batch" and "removing the effect you are looking for" are the same
# operation, and the algorithm cannot distinguish them. Running a corrector under
# that confounding produces a clean-looking UMAP with the finding deleted.
#
# So the honest order is: measure the batch effect, measure the CONFOUNDING
# between batch and cohort, and report both. If batch and cohort are separable
# and mixing is poor, correction is worth considering and the numbers here say so
# explicitly. If they are confounded, no correction is safe at any setting, and
# that is the single most important thing this diagnostic can tell you.
#
# THE METRIC. Local Inverse Simpson's Index (Korsunsky et al. 2019, the metric
# Harmony is evaluated with). For each cell, look at its k nearest neighbours and
# compute 1 / sum(p_b^2) over the batch proportions p_b among them. It runs from
# 1 (every neighbour shares the cell's batch -- no mixing) to the effective number
# of batches (neighbourhoods look like the dataset as a whole -- full mixing).
#
# WHY A PERMUTATION NULL AND NOT A THRESHOLD ON THE RAW SCORE: the achievable
# iLISI depends on how many batches there are and how unevenly sized they are, so
# "1.7" means nothing on its own. Shuffling the batch labels across cells and
# recomputing gives the score this dataset would produce with NO batch structure
# at all, which is the only reference that makes the observed value readable.

#' Local Inverse Simpson's Index of a label over a neighbourhood graph
#'
#' @param coords numeric matrix, cells x dimensions (the embedding, or scaled
#'   marker space)
#' @param labels character/factor vector, one per row of `coords`
#' @param k neighbourhood size
#' @return numeric vector of per-cell iLISI
#' @export
lisi_score <- function(coords, labels, k = 30L) {
  n <- nrow(coords)
  if (n < 3L) return(rep(NA_real_, n))
  k <- max(2L, min(as.integer(k), n - 1L))
  lab <- factor(labels)
  L <- nlevels(lab)
  if (L < 2L) return(rep(1, n))
  code <- as.integer(lab)
  out <- numeric(n)
  # Blocked distance computation: the full n x n matrix is what makes a naive
  # implementation unusable, so distances are taken a block of rows at a time and
  # discarded immediately. Peak memory is block x n, not n x n.
  block <- max(1L, min(512L, n))
  sq <- rowSums(coords^2)
  for (start in seq(1L, n, by = block)) {
    idx <- start:min(n, start + block - 1L)
    d2 <- outer(sq[idx], sq, "+") - 2 * (coords[idx, , drop = FALSE] %*% t(coords))
    for (r in seq_along(idx)) {
      dr <- d2[r, ]
      dr[idx[r]] <- Inf                     # a cell is not its own neighbour
      # PARTIAL selection, not a full order(). Only the IDENTITY of the k nearest
      # neighbours matters -- their ordering among themselves never enters the
      # Simpson index below -- so sorting all n distances is wasted work: O(n) for
      # a partial select against O(n log n) for a full sort. At 4,000 cells that
      # is 4,000 complete sorts of 4,000 elements per LISI evaluation, and the
      # permutation null runs twenty more evaluations on top.
      #
      # Two steps because sort.int() cannot return indices and take `partial` at
      # the same time: find the k-th smallest distance, then resolve which
      # elements reach it. `which()` yields candidates in index order and
      # order() is stable, so ties at the boundary break by index -- identical to
      # the order(dr)[seq_len(k)] this replaces.
      kth <- sort.int(dr, partial = k)[k]
      nb <- which(dr <= kth)
      if (length(nb) > k) nb <- nb[order(dr[nb])[seq_len(k)]]
      p <- tabulate(code[nb], nbins = L) / k
      out[idx[r]] <- 1 / sum(p^2)
    }
  }
  out
}

#' Batch-mixing diagnostic with a permutation null and a confounding check
#'
#' @param cells embedding cell table with umap_1/umap_2
#' @param batch_col column naming the batch (acquisition date, run, operator)
#' @param group_col biological grouping, used for the confounding check
#' @param k neighbourhood size
#' @param max_cells subsample ceiling; LISI is a local statistic and converges
#'   quickly, so a few thousand cells is ample and keeps this from dominating
#'   the runtime of the whole pipeline
#' @param n_perm permutations for the null
#' @param seed Random seed. The RNG stream is restored afterwards, so this cannot perturb later sampling. Default `42L`.
#' @return list(summary = data.frame, per_cell = numeric, confounding = data.frame)
#' @export
batch_mixing_report <- function(cells, batch_col, group_col = "cohort",
                                k = 30L, max_cells = 4000L, n_perm = 20L,
                                seed = 42L) {
  if (is.null(cells) || !nrow(cells)) return(NULL)
  if (!batch_col %in% names(cells)) return(NULL)
  if (!all(c("umap_1", "umap_2") %in% names(cells))) return(NULL)
  b <- as.character(cells[[batch_col]])
  ok <- !is.na(b) & nzchar(trimws(b)) &
        is.finite(cells$umap_1) & is.finite(cells$umap_2)
  if (sum(ok) < 100L) return(NULL)
  if (length(unique(b[ok])) < 2L) {
    log_msg("  batch diagnostic: only one level of '", batch_col, "' \u2014 skipped")
    return(NULL)
  }

  # RNG hygiene, for the same reason subcluster_by_reference() and fig_gating_qc()
  # do it: this runs mid-pipeline and must not shift which cells later steps draw.
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  idx <- which(ok)
  if (length(idx) > max_cells) idx <- sort(sample(idx, max_cells))
  co <- as.matrix(cells[idx, c("umap_1", "umap_2"), drop = FALSE])
  bl <- b[idx]

  obs <- lisi_score(co, bl, k = k)
  # Null: the same cells and the same neighbourhood graph, with batch labels
  # shuffled. Anything the geometry alone can produce shows up here.
  perm <- vapply(seq_len(max(1L, n_perm)), function(i)
    stats::median(lisi_score(co, sample(bl), k = k), na.rm = TRUE), numeric(1))

  med <- stats::median(obs, na.rm = TRUE)
  null_med <- stats::median(perm)
  # One-sided empirical p: how often does a batch-free arrangement mix as POORLY
  # as what was observed? Small p means the observed mixing is worse than chance,
  # i.e. a real batch effect. (+1 in numerator and denominator is the standard
  # correction that keeps an empirical p from ever being exactly 0.)
  p_emp <- (sum(perm <= med) + 1) / (length(perm) + 1)
  n_batch <- length(unique(bl))

  # ---- confounding between batch and biological group ----------------------
  conf <- NULL
  if (group_col %in% names(cells)) {
    g <- as.character(cells[[group_col]])[idx]
    keep <- !is.na(g) & nzchar(trimws(g))
    if (sum(keep) > 10L && length(unique(g[keep])) > 1L) {
      tb <- table(bl[keep], g[keep])
      # Cramer's V: association between two categoricals, scaled to \code{[0,1]} so it
      # does not depend on table size the way chi-squared does. 1 means each batch
      # contains exactly one group -- total confounding, no correction is safe.
      chi <- suppressWarnings(stats::chisq.test(tb)$statistic)
      nn <- sum(tb)
      v <- sqrt(as.numeric(chi) / (nn * (min(dim(tb)) - 1)))
      conf <- data.frame(
        batch_column = batch_col, group_column = group_col,
        n_batches = nrow(tb), n_groups = ncol(tb),
        cramers_v = round(v, 3),
        verdict = if (!is.finite(v)) "not estimable"
                  else if (v >= 0.8) "SEVERE - batch and group are near-identical; NO batch correction is safe"
                  else if (v >= 0.5) "substantial - batch and group overlap; correction would remove real signal"
                  else if (v >= 0.3) "moderate - interpret batch-adjacent findings with care"
                  else "low - batch and group are largely separable",
        stringsAsFactors = FALSE)
    }
  }

  summary <- data.frame(
    batch_column = batch_col, n_cells_scored = length(idx), k = k,
    n_batches = n_batch,
    ilisi_observed_median = round(med, 3),
    ilisi_null_median = round(null_med, 3),
    ilisi_max_possible = n_batch,
    mixing_fraction = round((med - 1) / max(1e-9, null_med - 1), 3),
    n_permutations = length(perm), p_empirical = round(p_emp, 4),
    verdict = if (p_emp > 0.05) "no detectable batch structure in the embedding"
              else if (med / null_med > 0.9) "statistically detectable but small batch structure"
              else "SUBSTANTIAL batch structure - cells sit next to same-batch cells",
    stringsAsFactors = FALSE)
  list(summary = summary, per_cell = obs, cell_index = idx, confounding = conf)
}

#' Batch diagnostic figure: the embedding coloured by batch, plus the LISI
#' distribution against its permutation null
#'
#' WHY BOTH PANELS: the scatter is what a reader will look at anyway and it shows
#' WHERE any batch structure sits; the LISI distribution is what says whether the
#' pattern the eye finds is more than chance. Neither alone is enough -- a
#' scattered-looking plot can still be significantly structured at these n, and a
#' significant score with no visible pattern is usually one small batch.
#' @param cells Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`, `sample_id`, `population_label` and one column per marker.
#' @param report The report.
#' @param outfile Path to write the figure to.
#' @param batch_col Column naming the acquisition batch.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_batch_diagnostic <- function(cells, report, outfile, batch_col,
                                 panel_label = "", dpi = 200, colors = fcs_colors()) {
  if (is.null(report)) return(invisible(NULL))
  d <- cells[report$cell_index, , drop = FALSE]
  d$.batch <- as.character(cells[[batch_col]])[report$cell_index]
  d$.lisi  <- report$per_cell

  lv <- sort(unique(d$.batch))
  cols <- population_colours(lv, colors = colors)
  aes_pt <- auto_point_aes(nrow(d))

  p1 <- ggplot(d, aes(umap_1, umap_2, colour = .batch)) +
    geom_point(size = aes_pt$size, alpha = aes_pt$alpha, stroke = 0) +
    scale_colour_manual(values = cols, name = batch_col) +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = paste0("Embedding by ", batch_col)) +
    theme_cyto(colors = colors) +
    theme(legend.position = "right", legend.key.size = unit(9, "pt"),
          legend.text = element_text(size = 7))

  s <- report$summary
  p2 <- ggplot(data.frame(lisi = d$.lisi), aes(lisi)) +
    geom_histogram(bins = 40, fill = colors$heatmap_low, colour = NA, alpha = 0.85) +
    geom_vline(xintercept = s$ilisi_observed_median, linewidth = 0.7,
               colour = colors$bracket) +
    geom_vline(xintercept = s$ilisi_null_median, linewidth = 0.7, linetype = "22",
               colour = colors$gate_highlight) +
    labs(x = paste0("iLISI (k = ", s$k, ")"), y = "cells",
         title = "Local batch mixing",
         subtitle = paste0("solid = observed median ", s$ilisi_observed_median,
                           "; dashed = permutation null ", s$ilisi_null_median,
                           "\n1 = no mixing, ", s$n_batches, " = full mixing")) +
    theme_cyto(colors = colors)

  cap <- paste0(
    "Empirical p = ", s$p_empirical, " over ", s$n_permutations,
    " label permutations (one-sided: is mixing WORSE than chance?). Verdict: ",
    s$verdict, ".")
  if (!is.null(report$confounding))
    cap <- paste0(cap, "\nBatch vs group confounding: Cramer's V = ",
                  report$confounding$cramers_v, " \u2014 ", report$confounding$verdict, ".")
  cap <- paste0(cap,
    "\nDIAGNOSTIC ONLY: no batch correction is applied by this pipeline. ",
    "Computed on the 2-D embedding, so it measures batch structure as PLOTTED.")

  fig <- patchwork::wrap_plots(list(p1, p2), ncol = 2, widths = c(1.35, 1)) +
    patchwork::plot_annotation(
      title = paste0("Batch-effect diagnostic",
                     if (nzchar(panel_label)) paste0(" \u2014 ", panel_label) else ""),
      caption = cap,
      theme = ggplot2::theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.caption = element_text(size = 7, hjust = 0, colour = colors$caption_text)))
  safe_ggsave(outfile, plot = fig, width = 11.5, height = 4.6, dpi = dpi,
              limitsize = FALSE)
  log_msg("[fig] wrote ", outfile, " (iLISI ", s$ilisi_observed_median,
          " vs null ", s$ilisi_null_median, ", p = ", s$p_empirical, ")")
  invisible(fig)
}


# =============================================================================
# THRESHOLD DRIFT: does the GATE move with cohort?
# =============================================================================
#
# THE GAP THIS CLOSES, AND WHY IT IS SPECIFIC TO THIS PACKAGE. A clustering-first
# workflow has no equivalent of this check because it never gates: one model is
# fitted to the pooled data, so every cell is judged by the same rule. This
# package derives a threshold PER SAMPLE from that sample's own density valley --
# which is what a human gater does, adapts correctly to staining variation, and
# is the right default.
#
# It also means "CD4-positive" is not literally the same predicate in every
# sample. If the derived CD4 threshold happens to sit systematically higher in
# one cohort, that cohort will show fewer CD4 T cells FOR THAT REASON ALONE, and
# the abundance test will report it as a biological difference with a small
# p-value. Nothing downstream can detect this, because by then the thresholds
# have been applied and discarded.
#
# The check is cheap and the pipeline already writes the data it needs
# (thresholds_used.csv, one row per sample x marker). If a marker's threshold
# does not differ across cohorts, per-sample gating cost nothing and the
# populations built on it are safe. If it does, that marker is flagged and every
# population whose definition uses it is named, so a reader knows exactly which
# results to treat as provisional.

#' Test each marker's per-sample threshold for a cohort difference
#'
#' @param thr thresholds_used table (sample_id, marker, threshold, source)
#' @param group_of named vector sample_id -> cohort
#' @param spec population spec, used to name which populations each flagged
#'   marker feeds into
#' @param min_n Minimum samples per group before a test is attempted. Default `3L`.
#' @export
stats_threshold_drift <- function(thr, group_of, spec = NULL, min_n = 3L) {
  if (is.null(thr) || !nrow(thr)) return(NULL)
  need <- c("sample_id", "marker", "threshold")
  if (!all(need %in% names(thr))) return(NULL)
  d <- thr[is.finite(thr$threshold), , drop = FALSE]
  d <- qc_pass_rows(d)
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  groups <- sort(unique(d$group))
  if (length(groups) < 2L) return(NULL)

  # Which populations would be affected if this marker's cut moves.
  uses <- list()
  if (!is.null(spec) && length(spec)) {
    for (pn in names(spec)) {
      r <- spec[[pn]]
      if (!is.list(r)) next
      mk <- setdiff(names(r), "any_of")
      if (!is.null(r$any_of) && is.list(r$any_of)) mk <- c(mk, names(r$any_of))
      for (m in unique(mk)) uses[[m]] <- c(uses[[m]], pn)
    }
  }

  rows <- list()
  for (mk in sort(unique(d$marker))) {
    dm <- d[d$marker == mk, , drop = FALSE]
    n_by <- table(dm$group)
    if (length(n_by) < 2L || any(n_by < min_n)) next
    kw <- try(stats::kruskal.test(dm$threshold, factor(dm$group)), silent = TRUE)
    if (inherits(kw, "try-error")) next
    med <- vapply(split(dm$threshold, dm$group), stats::median, numeric(1))
    # Spread expressed against the WITHIN-cohort variability, because a threshold
    # gap only matters relative to how much the cut wanders anyway. A gap of 0.3
    # asinh units is nothing if samples within a cohort already vary by 0.5, and
    # is decisive if they vary by 0.02.
    within_sd <- stats::median(vapply(split(dm$threshold, dm$group),
                                      function(v) stats::sd(v), numeric(1)),
                               na.rm = TRUE)
    gap <- max(med) - min(med)
    rows[[length(rows) + 1L]] <- data.frame(
      marker = mk, n_samples = nrow(dm), n_groups = length(n_by),
      median_gap_between_groups = round(gap, 4),
      typical_within_group_sd = round(within_sd, 4),
      gap_over_sd = round(gap / max(within_sd, 1e-9), 2),
      medians_by_group = paste(sprintf("%s=%.3f", names(med), med), collapse = "; "),
      dominant_source = names(sort(table(dm$source %||% "unknown"),
                                   decreasing = TRUE))[1] %||% NA_character_,
      test = "Kruskal-Wallis", p_value = kw$p.value,
      populations_using_this_marker =
        paste(sort(unique(uses[[mk]] %||% character(0))), collapse = "; "),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj_BH <- stats::p.adjust(out$p_value, method = "BH")
  out$drift_flag <- ifelse(
    !is.na(out$p_adj_BH) & out$p_adj_BH < 0.05 & out$gap_over_sd >= 1,
    "FLAGGED - threshold differs by cohort; abundance differences for the listed populations are partly definitional",
    ifelse(!is.na(out$p_value) & out$p_value < 0.05, "watch - nominal difference only", ""))
  out[order(out$p_value), , drop = FALSE]
}

#' Threshold-drift figure: per-sample cut for every marker, coloured by cohort
#' @param thr The thr.
#' @param outfile Path to write the figure to.
#' @param group_of Named character vector mapping sample_id to group label.
#' @param stats Statistics table from the matching stats_ function, used to annotate the figure.
#' @param reference The group every other group is compared against.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_threshold_drift <- function(thr, outfile, group_of, stats = NULL,
                                reference = NULL, panel_label = "", dpi = 200,
                                colors = fcs_colors()) {
  if (is.null(thr) || !nrow(thr)) return(invisible(NULL))
  d <- thr[is.finite(thr$threshold), , drop = FALSE]
  d <- qc_pass_rows(d)
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  glev <- sort(unique(d$group))
  if (!is.null(reference) && reference %in% glev)
    glev <- c(reference, setdiff(glev, reference))
  d$group <- factor(d$group, levels = glev)
  fills <- semantic_colours(glev, reference = reference, colors = colors)

  # Flagged markers first and marked in the facet label, so the eye lands on the
  # ones that matter instead of scanning an alphabetical grid.
  flagged <- character(0)
  if (!is.null(stats))
    flagged <- stats$marker[grepl("^FLAGGED", stats$drift_flag %||% "")]
  d$facet <- ifelse(d$marker %in% flagged, paste0(d$marker, "  [DRIFT]"), d$marker)
  ord <- unique(c(sort(unique(d$facet[d$marker %in% flagged])),
                  sort(unique(d$facet[!d$marker %in% flagged]))))
  d$facet <- factor(d$facet, levels = ord)

  fig <- ggplot(d, aes(group, threshold, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3,
                 colour = colors$bar_outline, alpha = 0.55) +
    geom_point(position = position_jitter(width = 0.14, height = 0), size = 1.1,
               colour = colors$bracket) +
    scale_fill_manual(values = fills, name = NULL) +
    facet_wrap(~ facet, scales = "free_y") +
    labs(x = NULL, y = "derived threshold (arcsinh units)",
         title = paste0("Gating-threshold drift across cohorts",
                        if (nzchar(panel_label)) paste0(" \u2014 ", panel_label) else ""),
         subtitle = paste0(
           "Each point is one sample's own derived cut for that marker. ",
           "Thresholds are derived PER SAMPLE by design, so some spread is expected ",
           "and healthy."),
         caption = paste0(
           "A marker labelled [DRIFT] has thresholds that differ systematically by cohort ",
           "(Kruskal-Wallis, BH-adjusted p < 0.05, gap exceeding within-cohort SD).\n",
           "For such a marker, part of any abundance difference in the populations it ",
           "defines is a difference in the DEFINITION \u2014 see threshold_drift_stats.csv ",
           "for which populations are affected.")) +
    theme_cyto(colors = colors) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "right", legend.key.size = unit(9, "pt"),
          legend.text = element_text(size = 7),
          strip.text = element_text(size = 7.5, face = "bold"),
          plot.caption = element_text(size = 7, hjust = 0, colour = colors$caption_text))

  nf <- length(levels(d$facet))
  nc <- max(1L, ceiling(sqrt(nf * 1.4)))
  nr <- ceiling(nf / nc)
  safe_ggsave(outfile, plot = fig, width = 1.55 * nc + 1.6, height = 1.5 * nr + 1.8,
              dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile, " (", nf, " markers, ",
          length(flagged), " flagged for drift)")
  invisible(fig)
}


# =============================================================================
# RUN MANIFEST: what produced this results folder
# =============================================================================
#
# THE GAP THIS CLOSES. The baseline captures sessionInfo() only inside
# save_session(), so a default run leaves no record of the R version, the package
# versions, the command line, or which input files were read. A methods section
# needs all four, and a results folder that cannot say what produced it cannot be
# reproduced or defended a year later.
#
# WHY IT IS WRITTEN FIRST AND UPDATED AT THE END, rather than written once when
# the run finishes: a run that CRASHES is exactly when you most want to know what
# it was doing. The first write happens before the expensive steps and records
# the invocation; the final write adds the outcome. A manifest whose status still
# reads "running" is itself information.

#' Write the run manifest
#'
#' @param path output file
#' @param opt parsed options
#' @param files input FCS paths (may be NULL on the first, pre-resolution write)
#' @param extra named list of run-derived values to record (cofactors, seeds,
#'   embedding parameters, sample counts)
#' @param status "running" | "completed" | "failed"
#' @param started Run start time.
#' @export
write_run_manifest <- function(path, opt = NULL, files = NULL, extra = NULL,
                               status = c("running", "completed", "failed"),
                               started = NULL) {
  status <- match.arg(status)
  ln <- c(
    "# cyRAVEN run manifest",
    "# Everything needed to say what produced this results folder.",
    paste0("status: ", status),
    paste0("started_utc: ", format(started %||% Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("written_utc: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "",
    "## environment",
    paste0("r_version: ", R.version.string),
    paste0("platform: ", R.version$platform),
    paste0("hostname: ", tryCatch(Sys.info()[["nodename"]], error = function(e) "unknown")),
    "")

  # Package versions for everything actually loaded. Recorded by querying the
  # namespace rather than the declared manifest, so what is written is what RAN.
  # Dependencies read from the package's OWN DESCRIPTION rather than from a
  # hand-maintained list. As a script there was no other source of truth and the
  # list had to be declared; as a package, DESCRIPTION already is that source,
  # and a second copy could only ever drift from it.
  deps <- tryCatch({
    d <- utils::packageDescription("cyRAVEN")
    fields <- unlist(d[c("Depends", "Imports", "Suggests")], use.names = FALSE)
    p <- unlist(strsplit(paste(stats::na.omit(fields), collapse = ","), ","))
    p <- trimws(sub("\\(.*\\)", "", p))
    setdiff(p[nzchar(p)], "R")
  }, error = function(e) character(0))
  pk <- unique(c(deps, "cyRAVEN"))
  ln <- c(ln, "## package versions (installed; '-' = not installed)")
  for (p in sort(pk)) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "-")
    ln <- c(ln, paste0("  ", p, ": ", v))
  }
  ln <- c(ln, "")

  # Pipeline provenance: the git commit if this is a checkout, so a results
  # folder can be traced to the exact code that made it.
  sha <- tryCatch({
    d <- script_dir()
    o <- suppressWarnings(system2("git", c("-C", shQuote(d), "rev-parse", "--short", "HEAD"),
                                  stdout = TRUE, stderr = FALSE))
    if (length(o) && nzchar(o[1])) o[1] else NA_character_
  }, error = function(e) NA_character_)
  dirty <- tryCatch({
    d <- script_dir()
    o <- suppressWarnings(system2("git", c("-C", shQuote(d), "status", "--porcelain"),
                                  stdout = TRUE, stderr = FALSE))
    length(o) > 0
  }, error = function(e) NA)
  ln <- c(ln, "## pipeline",
          paste0("script_dir: ", tryCatch(script_dir(), error = function(e) "unknown")),
          paste0("git_commit: ", sha %||% "not a git checkout"),
          paste0("git_uncommitted_changes: ", if (isTRUE(dirty)) "YES - results may not be reproducible from this commit" else "no"),
          "")

  ln <- c(ln, "## invocation",
          paste0("command: ", paste(commandArgs(trailingOnly = FALSE), collapse = " ")),
          "")
  if (!is.null(opt)) {
    ln <- c(ln, "## options")
    for (nm in sort(names(opt))) {
      v <- opt[[nm]]
      if (is.null(v) || (is.list(v) && !length(v))) next
      vs <- tryCatch(paste(utils::head(as.character(unlist(v)), 20), collapse = ", "),
                     error = function(e) "<unprintable>")
      ln <- c(ln, paste0("  ", nm, ": ", vs))
    }
    ln <- c(ln, "")
  }
  if (!is.null(files) && length(files)) {
    ln <- c(ln, paste0("## input files (", length(files), ")"))
    for (f in files) {
      sz <- tryCatch(file.size(f), error = function(e) NA)
      mt <- tryCatch(format(file.mtime(f), tz = "UTC"), error = function(e) NA)
      ln <- c(ln, paste0("  ", basename(f), "  bytes=", sz, "  mtime_utc=", mt))
    }
    ln <- c(ln, "")
  }
  if (!is.null(extra) && length(extra)) {
    ln <- c(ln, "## run-derived values")
    for (nm in names(extra)) {
      v <- extra[[nm]]
      vs <- tryCatch(paste(as.character(unlist(v)), collapse = ", "),
                     error = function(e) "<unprintable>")
      ln <- c(ln, paste0("  ", nm, ": ", vs))
    }
    ln <- c(ln, "")
  }
  ln <- c(ln, "## sessionInfo()",
          paste0("  ", tryCatch(utils::capture.output(utils::sessionInfo()),
                                error = function(e) "unavailable")))
  writeLines(ln, path)
  invisible(path)
}
