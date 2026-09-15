# =============================================================================
# SPLIT EVERY FIGURE BY TIMEPOINT
# =============================================================================
#
# WHAT THIS IS FOR. On a repeated-measures design the pooled figure is the wrong
# default for almost everything: a boxplot over all visits, a UMAP with every
# visit's cells superimposed, a heatmap whose columns mix d0 and d7, all answer a
# question about the cohort-as-a-whole that nobody asked, and hide the one the
# study exists to ask. Individual by-timepoint twins were being added one figure
# at a time, which is both laborious and unreliable -- a run writes around a
# hundred figures, and any figure added later starts out pooled again.
#
# This is the general mechanism instead: once the tables and the embedding are
# computed, re-emit the figure set once per timepoint into
# <outdir>/by_timepoint/<level>/, from the same artefacts, with the rows
# restricted to the samples drawn at that visit.
#
# WHY A SECOND PASS OVER THE SAME ARTEFACTS, AND NOT A PIPELINE RUN PER VISIT.
# Running the whole pipeline per timepoint would recompute the embedding from
# each visit's cells alone, and three embeddings computed separately are NOT
# comparable -- a position in one has no relationship to the same position in
# another, so "the d7 island moved" would be an artefact of running UMAP three
# times. This package computes one embedding across all samples precisely so
# positions can be compared, and the split has to preserve that. Gating and the
# thresholds are likewise left alone: re-deriving them per visit would make
# "CD4-positive" a different predicate at d0 and d7, which is the confound the
# threshold-drift check exists to detect.
#
# So: SHARED embedding, SHARED thresholds, SHARED statistics. Only the rows
# drawn are restricted. Every figure below is the same figure the pooled run
# writes, over one visit's samples.
#
# WHAT IS DELIBERATELY NOT SPLIT.
#   - The timepoint_* figures. They are already resolved by visit -- splitting a
#     per-patient trajectory by visit leaves one point per patient and destroys
#     the only thing it shows.
#   - clinical_timepoint / clinical_effects_timepoint. The timepoint IS the
#     variable under test there; holding it constant empties the figure.
#   - Between-group statistics. They are not recomputed per visit, and the
#     figures carry no brackets here: with 9/6/5 samples split further by study
#     group, every comparison is below --min-group-n, and drawing a test that
#     was computed on the pooled cohort next to one visit's data would attribute
#     a cohort-level result to a subset of it.

#' The timepoint levels present in a sample sheet, in collection order
#'
#' @param smap Sample map.
#' @param min_samples Levels with fewer samples than this are dropped: a figure
#'   set for a single acquisition is a set of one-column boxplots and one-sample
#'   UMAPs, which is noise rather than a view.
#' @return Character vector, or `NULL` when there is nothing to split by.
#' @keywords internal
timepoint_levels <- function(smap, min_samples = 2L) {
  if (is.null(smap) || !all(c("sample_id", "timepoint") %in% names(smap)))
    return(NULL)
  tv <- as.character(smap$timepoint)
  ok <- !is.na(tv) & nzchar(trimws(tv)) & !is.na(smap$sample_id)
  if (!any(ok)) return(NULL)
  tab <- table(tv[ok])
  keep <- names(tab)[tab >= min_samples]
  if (length(keep) < 2L) return(NULL)
  # Collection order, not alphabetical: d10 follows d7 rather than sitting
  # between d0 and d3.
  levels(natural_cluster_factor(keep))
}

#' The sample ids drawn at one timepoint
#' @param smap Sample map.
#' @param level One timepoint level.
#' @return Character vector of sample ids.
#' @keywords internal
samples_at_timepoint <- function(smap, level) {
  tv <- as.character(smap$timepoint)
  as.character(smap$sample_id[!is.na(tv) & tv == level])
}

#' Restrict any sample-keyed table to a set of samples
#'
#' Returns NULL rather than an empty frame when nothing is left, because every
#' figure function in this package treats NULL as "nothing to draw" and an empty
#' frame as "draw an empty figure".
#' @param x A data.frame carrying `sample_id`, or NULL.
#' @param ids Sample ids to keep.
#' @return The subset, or NULL.
#' @keywords internal
subset_by_sample <- function(x, ids) {
  if (is.null(x) || !nrow(x) || !"sample_id" %in% names(x)) return(NULL)
  out <- x[as.character(x$sample_id) %in% ids, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out
}

#' Write the figure set for one timepoint
#'
#' Every call is guarded individually: a figure that cannot be drawn from one
#' visit's samples -- a population with no cells at that visit, a marker absent
#' from those acquisitions -- must not take the rest of the set down with it.
#'
#' @param ctx Named list of the shared artefacts: `freq`, `mfi`, `fx`, `rt`,
#'   `cells`, `tc`, `ufreq`, `markers`, `covariates`, `feature_cols`,
#'   `panel_label`, `colors`.
#' @param dir Destination directory; created if absent.
#' @param label Timepoint label, added to every figure title.
#' @return Character vector of the files written.
#' @keywords internal
emit_timepoint_figures <- function(ctx, dir, label) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  written <- character(0)
  plab <- if (nzchar(ctx$panel_label %||% ""))
    paste0(ctx$panel_label, " - ", label) else label
  cols <- ctx$colors %||% fcs_colors()

  # One guarded emitter. The figure functions return their path or NULL, but
  # several return a ggplot object instead, so existence on disk is what counts
  # as written rather than the return value.
  emit <- function(name, expr) {
    f <- file.path(dir, name)
    ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
      log_msg("  WARNING ", label, "/", name, " failed: ", conditionMessage(e),
              ", every other figure is unaffected")
      FALSE
    })
    if (ok && file.exists(f)) written <<- c(written, f)
    invisible(NULL)
  }

  # ---- per-sample abundance and state --------------------------------------
  if (!is.null(ctx$freq)) {
    emit("population_frequencies.png",
         fig_population_frequencies(ctx$freq,
           file.path(dir, "population_frequencies.png"),
           panel_label = plab, colors = cols))
    if ("cells_absolute" %in% names(ctx$freq) &&
        any(is.finite(ctx$freq$cells_absolute)))
      emit("absolute_vs_share.png",
           fig_absolute_vs_share(ctx$freq,
             file.path(dir, "absolute_vs_share.png"), panel_label = plab))
  }
  if (!is.null(ctx$mfi)) {
    emit("population_marker_heatmap.png",
         fig_population_marker_heatmap(ctx$mfi,
           file.path(dir, "population_marker_heatmap.png"),
           panel_label = plab, colors = cols))
    emit("marker_state.png",
         fig_marker_state(ctx$mfi, file.path(dir, "marker_state.png"),
           panel_label = plab, colors = cols))
  }
  if (!is.null(ctx$fx))
    emit("functional_markers.png",
         fig_functional_markers(ctx$fx,
           file.path(dir, "functional_markers.png"),
           panel_label = plab, colors = cols))
  if (!is.null(ctx$rt))
    emit("population_ratios.png",
         fig_population_ratios(ctx$rt,
           file.path(dir, "population_ratios.png"),
           panel_label = plab, colors = cols))
  if (!is.null(ctx$tc))
    emit("total_counts_qc.png",
         fig_total_counts_qc(ctx$tc, file.path(dir, "total_counts_qc.png"),
           panel_label = plab))
  if (!is.null(ctx$ufreq)) {
    emit("frequency_uncertainty.png",
         fig_frequency_uncertainty(ctx$ufreq,
           file.path(dir, "frequency_uncertainty.png"), colors = cols))
    emit("detection_limits.png",
         fig_detection_limits(ctx$ufreq,
           file.path(dir, "detection_limits.png"), colors = cols))
  }

  # ---- the shared embedding, one visit's cells ------------------------------
  # THE COORDINATES ARE NOT RECOMPUTED. These are the same positions as the
  # pooled figures, with the other visits' points not drawn, so an island in
  # one visit's panel is the same island in another's.
  if (!is.null(ctx$cells) && nrow(ctx$cells)) {
    emit("umap_overview.png",
         fig_umap_overview(ctx$cells, file.path(dir, "umap_overview.png"),
           panel_label = plab, covariates = ctx$covariates,
           feature_cols = ctx$feature_cols %||% character(0), colors = cols))
    emit("umap_density.png",
         fig_density_by_sample(ctx$cells, file.path(dir, "umap_density.png"),
           panel_label = plab, colors = cols))
    if (length(ctx$markers))
      emit("umap_markers.png",
           fig_marker_grid(ctx$cells, ctx$markers,
             file.path(dir, "umap_markers.png"), panel_label = plab,
             colors = cols))
    if ("population_label" %in% names(ctx$cells))
      emit("cohort_composition_heatmap.png",
           fig_cohort_confusion(ctx$cells,
             file.path(dir, "cohort_composition_heatmap.png"),
             group_col = "sample_id", panel_label = plab, colors = cols))
  }
  written
}

#' Re-emit the whole figure set once per timepoint
#'
#' @param ctx Shared artefacts; see [emit_timepoint_figures()]. Must also carry
#'   `smap`.
#' @param outdir Run output directory; figures go to `by_timepoint/<level>/`.
#' @return Character vector of files written, invisibly.
#' @keywords internal
split_figures_by_timepoint <- function(ctx, outdir) {
  lv <- timepoint_levels(ctx$smap)
  if (is.null(lv)) {
    log_msg("  --split-by-timepoint: the sample sheet carries no timepoint ",
            "with two or more levels of at least two samples; nothing to split")
    return(invisible(character(0)))
  }
  root <- file.path(outdir, "by_timepoint")
  written <- character(0)
  for (tp in lv) {
    ids <- samples_at_timepoint(ctx$smap, tp)
    sub <- ctx
    for (nm in c("freq", "mfi", "fx", "rt", "tc", "ufreq"))
      sub[[nm]] <- subset_by_sample(ctx[[nm]], ids)
    sub$cells <- if (!is.null(ctx$cells) && "sample_id" %in% names(ctx$cells))
      ctx$cells[as.character(ctx$cells$sample_id) %in% ids, , drop = FALSE]
      else NULL
    w <- emit_timepoint_figures(sub, file.path(root, split_dir_name(tp)), tp)
    written <- c(written, w)
    log_msg("  ", tp, ": ", length(ids), " sample(s), ", length(w),
            " figure(s) -> by_timepoint/", split_dir_name(tp))
  }
  log_msg("--split-by-timepoint: ", length(written), " figure(s) over ",
          length(lv), " timepoint(s). Same embedding, same thresholds and the ",
          "same statistics as the pooled figures -- only the rows drawn differ, ",
          "so positions and cuts are comparable between visits.")
  invisible(written)
}

#' A filesystem-safe directory name for a split level
#' @param x Character.
#' @keywords internal
split_dir_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
