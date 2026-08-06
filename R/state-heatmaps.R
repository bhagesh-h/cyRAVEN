# =============================================================================
# Differential state, and the phenotype / composition heatmaps
# =============================================================================
#
# WHAT THIS FILE IS: the analyses that read the scored populations back and ask
# what is INSIDE them -- whether a marker sits at a different level between
# cohorts, whether each population carries the phenotype its label claims, and
# which cohorts each population is actually drawn from. Nothing here gates,
# scores or embeds; every function below consumes tables the earlier stages
# produced and is therefore safe to run, re-run or omit without changing them.
#
# WHAT DRIVES THE CONTENTS: a function-by-function comparison against the
# established general-purpose cytometry analysis toolboxes this package is
# closest in scope to. Their algorithm chain is
#   read -> PCA -> UMAP -> Phenograph/FlowSOM -> metaclustering
#   -> batch harmonisation -> per-sample aggregation -> diffcyt DA/DS
# and the gaps identified are implemented below in order of scientific impact.
#
# HOUSE RULES FOLLOWED HERE, because they are what the package is built on:
#   - every function carries WHAT it does and WHY it is done that way
#   - a statistic that cannot be computed honestly returns NULL, never a
#     plausible-looking number
#   - replicates are SAMPLES, never cells; anything cell-level is labelled
#     descriptive and carries no p-value
#   - figures reuse the package's theme/colour/layout machinery so a reader
#     does not have to learn a second visual language
# =============================================================================


# =============================================================================
# DIFFERENTIAL STATE: does a marker sit at a different level inside the
# same population, between cohorts?
# =============================================================================
#
# THE GAP THIS CLOSES. Abundance testing can tell you a population changed SIZE
# (group_comparison.png, Wilcoxon on per-sample abundance -- correct, sample-level
# statistics). The overlay figures can SHOW you a marker shifted inside a
# population (umap_multigraph_overlay.png) and RANK those shifts
# (subcluster_marker_shifts.csv). What none of that can do is TEST one, because
# the ranking is computed over pooled cells: n is the number of cells, so one
# donor contributing 4,000 cells can produce a large Cliff's delta on their own.
# That is pseudoreplication. It is a common enough mistake that the widely used
# cell-level Wilcoxon helpers in this field ship with a warning telling you not
# to interpret their p-values.
#
# THE FIX IS THE ESTABLISHED ONE. diffcyt's differential-state pathway
# (Weber et al. 2019, Commun Biol 2:183) aggregates each marker to ONE VALUE PER
# SAMPLE PER POPULATION first (calcMedians), then models those sample-level
# values. n becomes the number of donors, which is what the design actually
# provides. This pipeline already computes that table -- population_marker_mfi.csv
# is exactly per sample x population x marker medians -- so the aggregation step
# is already done and only the testing was missing.
#
# WHY NOT limma/LMM AS diffcyt USES: limma's moderated t assumes approximate
# normality of the per-sample values and borrows variance across markers via an
# empirical-Bayes prior. That is a good trade at the marker counts of a mass
# cytometry panel (30-40) and the sample counts of a designed experiment. At the
# sizes here -- a dozen markers, single-digit donors per cohort -- the prior is
# estimated from too few markers to help, and the normality assumption is doing
# real work rather than being a formality. The baseline already made this call
# for abundance (see stats_group_comparison) and settled on Wilcoxon rank-sum
# plus Kruskal-Wallis, which assume only independence and ordinal comparability.
# Using the SAME tests here means an abundance result and a state result are
# directly comparable and a reader learns one set of caveats, not two.

#' Composite panel key for population x marker tables
#'
#' WHY: fig_group_comparison() panels on `population`. A marker only means
#' something inside the population it was measured in -- CD16 in monocytes and
#' CD16 in NK cells are different biology -- so the panel unit is the PAIR.
#' Mirrors fx_panel_key() in the baseline, which solved the same problem for
#' functional blocks.
#' @param mfi Per-sample x population x marker summary table.
#' @keywords internal
mfi_panel_key <- function(mfi) {
  mfi$population <- paste0(mfi$population, ", ", mfi$marker)
  mfi
}

#' Sample-level differential-state test: population x marker x group
#'
#' WHAT: for every (population, marker) pair, compares each non-reference cohort
#' against the reference using the per-sample values already in
#' population_marker_mfi.csv. One row per pair per comparison.
#'
#' WHY BOTH MEASURES ARE TESTED AND REPORTED SEPARATELY:
#'   `median_asinh`  the population's median intensity for that marker. Moves
#'                   when the whole population shifts up or down -- the classic
#'                   "activation marker is brighter in patients" result.
#'   `pct_positive`  the fraction of that population above the sample's own
#'                   gating threshold. Moves when a SUBSET of the population
#'                   turns positive while the rest does not.
#' A bimodal shift (30% of cells go bright, 70% unchanged) moves pct_positive
#' sharply and the median barely at all; a uniform shift does the opposite.
#' Reporting only one measure makes the other kind of change invisible, so both
#' are computed and the measure is a column, not an assumption.
#'
#' WHY min_cells: the per-sample median of a marker over 12 cells is noise
#' shaped like a measurement. Rows below the floor are dropped BEFORE testing
#' rather than down-weighted, because a rank test has no weighting mechanism and
#' silently treats a 12-cell median as equal evidence to a 4,000-cell one.
#'
#' WHY THE FDR FAMILY IS THE WHOLE TABLE: this tests every marker in every
#' population, which at a 12-marker panel and 12 populations is ~144 pairs x the
#' number of comparisons. That is a screen, and its multiplicity has to be
#' accounted for across the screen, not within each population.
#'
#' @param mfi population_marker_mfi table (sample_id, population, marker,
#'   n_cells, median_asinh, pct_positive, and -- after the deliverable change
#'   documented in the README -- is_control and qc_status)
#' @param group_of named character vector sample_id -> group label
#' @param reference group every other group is compared against
#' @param measures which per-sample quantities to test
#' @param min_cells per-sample floor on the population's cell count
#' @param min_n minimum samples per group for a test to be attempted
#' @return data.frame, one row per population x marker x comparison x measure,
#'   or NULL when nothing is testable
#' @export
stats_marker_state <- function(mfi, group_of, reference = NULL,
                               measures = c("median_asinh", "pct_positive"),
                               min_cells = 20L, min_n = 3L) {
  if (is.null(mfi) || !nrow(mfi)) return(NULL)
  need <- c("sample_id", "population", "marker", "n_cells")
  if (!all(need %in% names(mfi))) {
    log_msg("NOTE differential state skipped: population_marker_mfi is missing ",
            paste(setdiff(need, names(mfi)), collapse = ", "))
    return(NULL)
  }
  measures <- intersect(measures, names(mfi))
  if (!length(measures)) return(NULL)

  # Controls and staining-QC failures out first, for exactly the reason
  # qc_pass_rows() exists: their marker medians are real numbers describing a
  # tube that was not stained, or was stained and failed.
  d0 <- qc_pass_rows(mfi)
  dropped <- nrow(mfi) - nrow(d0)
  d0 <- d0[is.finite(d0$n_cells) & d0$n_cells >= min_cells, , drop = FALSE]
  if (!nrow(d0)) {
    log_msg("NOTE differential state skipped: no population x marker row has ",
            min_cells, "+ cells in a QC-passing sample")
    return(NULL)
  }
  d0$group <- unname(group_of[d0$sample_id])
  d0 <- d0[!is.na(d0$group), , drop = FALSE]
  if (!nrow(d0)) return(NULL)
  groups <- sort(unique(d0$group))
  if (length(groups) < 2L) return(NULL)
  if (is.null(reference) || !reference %in% groups) reference <- groups[1]
  others <- setdiff(groups, reference)

  rows <- list()
  for (meas in measures) {
    d <- d0[is.finite(d0[[meas]]), , drop = FALSE]
    if (!nrow(d)) next
    # A control character as the key separator, not a printable string: a
    # population name that happens to contain the separator would otherwise
    # merge two different population x marker pairs into one test.
    d$pair <- paste(d$population, d$marker, sep = "\r")
    for (pk in sort(unique(d$pair))) {
      dp <- d[d$pair == pk, , drop = FALSE]
      vals <- split(dp[[meas]], dp$group)
      n_by <- vapply(groups, function(g) length(vals[[g]] %||% numeric(0)), integer(1))

      # Omnibus first with 3+ groups, so a reader sees whether ANY cohort differs
      # before reading pairwise columns -- same order of presentation the
      # abundance table uses.
      p_omni <- NA_real_
      if (length(groups) > 2L && all(n_by >= min_n)) {
        kw <- try(stats::kruskal.test(dp[[meas]], factor(dp$group)), silent = TRUE)
        if (!inherits(kw, "try-error")) p_omni <- kw$p.value
      }
      for (g in others) {
        a <- vals[[reference]]; b <- vals[[g]]
        if (is.null(a) || is.null(b) || length(a) < min_n || length(b) < min_n) next
        wt <- try(stats::wilcox.test(a, b, exact = FALSE), silent = TRUE)
        p <- if (inherits(wt, "try-error")) NA_real_ else wt$p.value
        # Cliff's delta over SAMPLES (not cells). Bounded \co-1 to 1, rank-based, and
        # at single-digit n far more stable than a ratio of means.
        cliff <- mean(outer(b, a, ">")) - mean(outer(b, a, "<"))
        rows[[length(rows) + 1L]] <- data.frame(
          population = dp$population[1], marker = dp$marker[1], measure = meas,
          reference_group = reference, comparison_group = g,
          n_reference = length(a), n_comparison = length(b),
          median_reference  = stats::median(a),
          median_comparison = stats::median(b),
          # DIFFERENCE, not ratio. median_asinh is on an arcsinh scale that runs
          # negative, where a ratio is uninterpretable and changes sign for no
          # biological reason. pct_positive gets the same treatment so one column
          # means one thing throughout the table.
          delta_median = stats::median(b) - stats::median(a),
          cells_reference  = sum(dp$n_cells[dp$group == reference]),
          cells_comparison = sum(dp$n_cells[dp$group == g]),
          cliffs_delta = round(cliff, 3),
          test = "Wilcoxon rank-sum (samples)", p_value = p,
          p_omnibus_kruskal = p_omni, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) {
    log_msg("NOTE differential state produced no testable pair (need ", min_n,
            "+ samples per group with ", min_cells, "+ cells)")
    return(NULL)
  }
  out <- do.call(rbind, rows)
  # FDR WITHIN each measure, not across both: median_asinh and pct_positive are
  # two views of the same cells and are strongly correlated, so pooling them
  # would inflate the family size with tests that are not independent questions
  # and over-penalise every hit.
  out$p_adj_BH <- NA_real_
  for (meas in unique(out$measure)) {
    i <- out$measure == meas
    out$p_adj_BH[i] <- stats::p.adjust(out$p_value[i], method = "BH")
  }
  out$significant_raw <- !is.na(out$p_value) & out$p_value < 0.05
  out$significant_BH  <- !is.na(out$p_adj_BH) & out$p_adj_BH < 0.05
  if (dropped > 0)
    log_msg("  differential state: ", dropped,
            " control/QC-failed row(s) excluded before testing")
  out[order(out$measure, out$p_value), , drop = FALSE]
}

#' Differential-state figure, drawn by the baseline's own comparison figure
#'
#' WHY IT DELEGATES RATHER THAN DRAWING: fig_group_comparison() already renders
#' exactly this shape -- one panel per unit, bar = mean, rule = median, whiskers =
#' SD, points = samples, Wilcoxon brackets. Reusing it means the DS figure cannot
#' drift from the abundance figure's layout, and a reader who has learned to read
#' one has learned to read the other. Exactly the pattern
#' fig_functional_markers() already uses.
#'
#' WHY pct_positive AND NOT median_asinh IS PLOTTED: fig_group_comparison() fixes
#' the y axis at c(0, headroom) -- deliberately, since every quantity it was built
#' for (percentages, concentrations, ratios) is non-negative. Arcsinh medians run
#' negative for dim markers, and handing them to that scale would silently clip
#' the bar to zero and draw a confident-looking panel of a value that is not
#' there. pct_positive is bounded 0 to 100, is a genuine differential-state
#' measure, and is what a bimodal shift moves most. The median_asinh results are
#' not discarded -- they carry the full test in the CSV and are the quantity the
#' state heatmap in fig_population_marker_heatmap() draws, where a diverging scale handles sign correctly.
#' @param mfi Per-sample x population x marker summary table.
#' @param outfile Path to write the figure to.
#' @param group_of Named character vector mapping sample_id to group label.
#' @param stats Statistics table from the matching stats_ function, used to annotate the figure.
#' @param reference The group every other group is compared against.
#' @param p_source Which p-value the figure annotates: raw or BH-adjusted. Default `c("raw", "BH")`.
#' @param min_cells Minimum cells a sample must contribute before it is used. Default `20L`.
#' @param ncol Number of panel columns; NULL computes one that keeps the canvas roughly square.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_marker_state <- function(mfi, outfile, group_of = NULL, stats = NULL,
                             reference = NULL, p_source = c("raw", "BH"),
                             min_cells = 20L, ncol = NULL, panel_label = "",
                             dpi = 300, colors = fcs_colors()) {
  p_source <- match.arg(p_source)
  if (is.null(mfi) || !nrow(mfi) || !"pct_positive" %in% names(mfi))
    return(invisible(NULL))
  d <- mfi[is.finite(mfi$n_cells) & mfi$n_cells >= min_cells, , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  d <- mfi_panel_key(d)

  # COLUMN COUNT IS COMPUTED, NOT FIXED. This figure has one panel per
  # population x marker pair, so a 12-population 14-marker panel asks for ~168 of
  # them -- an order of magnitude more than the abundance figure this delegates
  # to, which was tuned for about a dozen. At a fixed 4 columns that came out 65
  # INCHES TALL and 12 wide: a ribbon nobody can read, print or place in a
  # document, and one that wastes most of its own canvas.
  #
  # Solving for a roughly square canvas instead: panels are 2.75in wide by 2.15in
  # tall, so the column count that squares the figure is sqrt(n * 2.15/2.75).
  # The 168-panel case becomes ~33 x 31in -- still a large figure, because 168
  # panels IS large, but a balanced one that scales down to a page.
  #
  # NO SILENT CAP: every pair is drawn. Truncating to the top-N by p-value would
  # make the figure's contents depend on the test result, which is exactly the
  # kind of selection this pipeline avoids elsewhere. If it is too big to read,
  # marker_state_stats.csv is the sorted artefact to read instead.
  if (is.null(ncol)) {
    npair <- length(unique(d$population))
    ncol <- max(4L, as.integer(ceiling(sqrt(npair * 2.15 / 2.75))))
  }
  st <- NULL
  if (!is.null(stats)) {
    st <- stats[stats$measure == "pct_positive", , drop = FALSE]
    if (nrow(st)) st$population <- paste0(st$population, ", ", st$marker)
    else st <- NULL
  }
  if (is.null(group_of)) {
    ids <- unique(d$sample_id)
    group_of <- setNames(rep("All samples", length(ids)), ids)
    st <- NULL
  }
  fig_group_comparison(
    d, outfile, group_of = group_of, stats = st, reference = reference,
    p_source = p_source, ncol = ncol, panel_label = panel_label, dpi = dpi,
    value_col = "pct_positive", value_label = "% positive",
    value_caveat = paste(
      "Differential state: every marker in every population, tested on",
      "PER-SAMPLE values (n = donors, not cells). Positivity is the fraction of",
      "the population above that sample's own gating threshold;",
      sprintf("populations with fewer than %d cells in a sample are omitted.", min_cells),
      "Intensity shifts (median arcsinh) are tested in marker_state_stats.csv",
      "and drawn in population_marker_heatmap.png."),
    title_noun = "Marker state", colors = colors)
}


# =============================================================================
# PHENOTYPE HEATMAPS: are the populations what their labels claim, and
# which cohort fills each one?
# =============================================================================
#
# THE GAP THIS CLOSES. A marker-by-population heatmap -- markers in rows,
# populations in columns, z-scored across populations -- is the standard
# annotation figure in this field. The gating stage already produces the data for
# it (population_marker_mfi.csv); what was missing was the figure.
#
# WHY IT MATTERS MORE HERE THAN IN A CLUSTERING-FIRST WORKFLOW. There, populations
# come from unsupervised clustering and the heatmap's job is to NAME them. Here
# they come
# from Boolean threshold gates that were written down in advance, so the heatmap's
# job is to CHECK them: if the column labelled "CD4 T cells" does not show high
# CD4 and low CD8, the gate is wrong, and that is visible in one glance instead
# of being inferred from a frequency that looks surprising. A pipeline that
# assigns labels from thresholds needs a figure that can contradict the labels.
#
# WHY z-SCORE ACROSS POPULATIONS (each marker row centred and scaled over the
# population columns): markers differ by orders of magnitude in raw brightness,
# so an unscaled heatmap is a picture of which fluorophore is brightest. Scaling
# within the row asks the only question worth asking -- for THIS marker, which
# populations are high and which are low -- which is what makes the gate check
# work. The unscaled medians are written alongside in the CSV.

#' Hard-wrap a title/subtitle/caption to the figure's actual width
#'
#' WHY THIS IS NEEDED AT ALL: ggplot2 does not wrap plot titles, subtitles or
#' captions. A long one is drawn as a single line and CLIPPED at the device edge,
#' silently -- the text is simply gone, with no warning and no visual cue that
#' anything is missing. On these heatmaps that lost the last third of both
#' subtitles, including the sentence explaining what an even split is.
#'
#' It is worse than an ordinary layout bug because the truncation lands
#' mid-sentence and reads as if the sentence ended there. A reader has no way to
#' tell a clipped subtitle from a badly written one.
#'
#' HOW THE WIDTH IS ESTIMATED: at a given point size the mean advance width of
#' the default sans face is close to 0.52 * size, so characters-per-inch is
#' 72 / (0.52 * pt). `margin_in` covers the plot margins and any inset. This is
#' an estimate, not a measurement -- strwrap is given a deliberately conservative
#' character count so a slightly wider-than-average string still fits rather than
#' spilling. Callers pair it with plot.title.position = "plot", which starts the
#' text at the figure edge instead of the panel edge and so makes the full width
#' available.
#'
#' @param txt the string; any newlines already in it are preserved and each
#'   resulting line is wrapped independently (captions rely on this)
#' @param width_in the figure width being passed to ggsave
#' @param pt point size the text will be rendered at
#' @param margin_in Margin to reserve, in inches. Default `0.35`.
#' @keywords internal
wrap_plot_text <- function(txt, width_in, pt = 9.5, margin_in = 0.35) {
  if (!length(txt) || !nzchar(txt)) return(txt)
  usable <- max(1, width_in - margin_in)
  n <- max(20L, as.integer(usable * 72 / (0.52 * pt)))
  parts <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  paste(vapply(parts, function(p)
    paste(strwrap(p, width = n), collapse = "\n"), character(1)),
    collapse = "\n")
}

#' Pick black or white text per tile, whichever the fill can actually be read
#' against
#'
#' WHY: a single fixed text colour cannot work on a sequential fill. The
#' composition heatmap used one grey for every label, so numbers sitting on the
#' dark end of viridis were effectively invisible -- grey on near-black -- while
#' the same grey was fine on the bright end. Which cells become unreadable
#' depends on the data, so it is not something a fixed colour can be chosen
#' around.
#'
#' HOW: compute the fill each value will actually receive, take its WCAG relative
#' luminance, and return whichever of white / near-black has the higher contrast
#' ratio against it. Comparing the two ratios rather than thresholding on
#' luminance means the answer stays correct if the palette is re-themed -- the
#' viridis option is configurable, and "cividis" or "magma" put their dark end in
#' a different place than "D" does.
#'
#' @param values numeric vector
#' @param limits the fill scale's limits; MUST match the scale, so the colour
#'   computed here is the colour ggplot draws
#' @param option viridis option name
#' @param dark the dark ink to use when the fill is light
#' @keywords internal
contrast_text_colour <- function(values, limits, option = "D", dark = "grey10") {
  n <- length(values)
  if (!n) return(character(0))
  if (!requireNamespace("viridisLite", quietly = TRUE))
    return(rep(dark, n))                       # no palette to reason about
  ramp <- viridisLite::viridis(256, option = option)
  frac <- (values - limits[1]) / max(1e-9, diff(limits))
  frac[!is.finite(frac)] <- 0
  fills <- ramp[1L + round(255 * pmin(1, pmax(0, frac)))]

  # WCAG 2.x relative luminance: linearise each sRGB channel, then weight.
  rel_lum <- function(cols) {
    m <- grDevices::col2rgb(cols) / 255
    lin <- ifelse(m <= 0.03928, m / 12.92, ((m + 0.055) / 1.055)^2.4)
    as.numeric(0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ])
  }
  ratio <- function(a, b) (pmax(a, b) + 0.05) / (pmin(a, b) + 0.05)
  L  <- rel_lum(fills)
  ifelse(ratio(L, rel_lum("white")) >= ratio(L, rel_lum(dark)), "white", dark)
}

#' Population x marker phenotype heatmap
#'
#' @param mfi population_marker_mfi table
#' @param scale_by "marker" (row z-score, the annotation/validation view),
#'   "none" (raw arcsinh medians, comparable across populations within a marker
#'   but not between markers)
#' @param min_cells per-sample floor before a sample contributes to a cell
#' @param annotate_expected optional named list population -> character vector of
#'   markers the gate definition REQUIRES to be positive. Drawn as an outline on
#'   those cells, turning the figure into an explicit gate audit.
#' @param outfile Path to write the figure to.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_population_marker_heatmap <- function(mfi, outfile, scale_by = c("marker", "none"),
                                          min_cells = 20L, annotate_expected = NULL,
                                          panel_label = "", dpi = 300,
                                          colors = fcs_colors()) {
  scale_by <- match.arg(scale_by)
  if (is.null(mfi) || !nrow(mfi)) return(invisible(NULL))
  d <- qc_pass_rows(mfi)
  d <- d[is.finite(d$n_cells) & d$n_cells >= min_cells &
         is.finite(d$median_asinh), , drop = FALSE]
  if (!nrow(d)) {
    log_msg("[fig] no QC-passing population x marker cell with ", min_cells,
            "+ cells, phenotype heatmap skipped")
    return(invisible(NULL))
  }

  # Aggregate over samples with the MEDIAN of per-sample medians, not the mean of
  # pooled cells: the latter weights a deeply-acquired sample more heavily, which
  # is the same pooling mistake the DS test above exists to avoid.
  agg <- stats::aggregate(median_asinh ~ population + marker, d, stats::median)
  n_s <- stats::aggregate(sample_id ~ population + marker, d,
                          function(x) length(unique(x)))
  names(n_s)[3] <- "n_samples"
  agg <- merge(agg, n_s, by = c("population", "marker"), all.x = TRUE)

  agg$value <- agg$median_asinh
  legend_lab <- "median arcsinh"
  if (scale_by == "marker") {
    agg$value <- unsplit(lapply(split(agg$median_asinh, agg$marker), function(v) {
      s <- stats::sd(v)
      # A marker with no variation between populations is not evidence of
      # anything; centring it and dividing by ~0 would manufacture extreme
      # z-scores out of rounding. Report it flat and let the reader see it is flat.
      if (!is.finite(s) || s < 1e-8) rep(0, length(v)) else (v - mean(v)) / s
    }), agg$marker)
    legend_lab <- "z-score\n(within marker)"
  }

  # Order populations by overall brightness so related lineages land near each
  # other, and markers alphabetically so a reader can find one by name.
  pop_ord <- stats::aggregate(value ~ population, agg, mean)
  agg$population <- factor(agg$population,
                           levels = pop_ord$population[order(-pop_ord$value)])
  agg$marker <- factor(agg$marker, levels = sort(unique(as.character(agg$marker))))

  lim <- max(abs(agg$value), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1

  fig <- ggplot(agg, aes(population, marker, fill = value)) +
    geom_tile(colour = colors$tile_border, linewidth = 0.3)

  # Gate audit: outline the cells the population's definition requires to be
  # positive. A required cell that is NOT bright is a gate that did not do what
  # the config says it does.
  if (!is.null(annotate_expected) && length(annotate_expected)) {
    exp_df <- do.call(rbind, lapply(names(annotate_expected), function(p) {
      mk <- intersect(annotate_expected[[p]], levels(agg$marker))
      if (!length(mk) || !p %in% levels(agg$population)) return(NULL)
      data.frame(population = p, marker = mk, stringsAsFactors = FALSE)
    }))
    if (!is.null(exp_df) && nrow(exp_df)) {
      exp_df$population <- factor(exp_df$population, levels = levels(agg$population))
      exp_df$marker <- factor(exp_df$marker, levels = levels(agg$marker))
      fig <- fig + geom_tile(data = exp_df, aes(population, marker),
                             fill = NA, colour = colors$gate_highlight,
                             linewidth = 0.8, inherit.aes = FALSE)
    }
  }

  # Canvas width first: the subtitle has to be wrapped to a width that is
  # already known (see wrap_plot_text).
  w <- max(6.0, 0.42 * length(levels(agg$population)) + 3.4)
  h <- max(4.0, 0.28 * length(levels(agg$marker)) + 2.6)

  fig <- fig +
    scale_fill_gradient2(low = colors$heatmap_low, mid = colors$heatmap_mid,
                         high = colors$heatmap_high,
                         midpoint = if (scale_by == "marker") 0 else stats::median(agg$value),
                         limits = if (scale_by == "marker") c(-lim, lim) else NULL,
                         name = legend_lab) +
    labs(x = NULL, y = NULL,
         title = paste0("Population phenotype",
                        if (nzchar(panel_label)) paste0(", ", panel_label) else ""),
         subtitle = wrap_plot_text(paste0(
           "Median of per-sample median arcsinh intensity",
           if (scale_by == "marker")
             ", z-scored within each marker row across populations" else "",
           ".",
           if (!is.null(annotate_expected) && length(annotate_expected))
             " Outlined cells are markers the gate definition requires to be POSITIVE." else ""),
           w, pt = 9.5),
         caption = wrap_plot_text(paste0(
           "Aggregated across QC-passing samples only; a population x marker cell needs ",
           min_cells, "+ cells in a sample for that sample to contribute.\n",
           "GATE AUDIT: a population whose defining markers are not bright here has a ",
           "threshold problem, not a biological one, check thresholds_used.csv."),
           w, pt = 7)) +
    theme_cyto(colors = colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          panel.grid = element_blank(),
          plot.title.position = "plot",
          plot.subtitle = element_text(lineheight = 1.15),
          plot.caption = element_text(size = 7, hjust = 0, lineheight = 1.2,
                                      colour = colors$caption_text),
          legend.key.height = unit(18, "pt"))

  h <- h + 0.16 * (lengths(regmatches(fig$labels$subtitle,
                                      gregexpr("\n", fig$labels$subtitle))) +
                   lengths(regmatches(fig$labels$caption,
                                      gregexpr("\n", fig$labels$caption))))
  safe_ggsave(outfile, plot = fig, width = w, height = h, dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile, " (", length(levels(agg$population)),
          " populations x ", length(levels(agg$marker)), " markers)")
  invisible(agg)
}

#' Cohort composition heatmap ("confusion matrix")
#'
#' WHAT: for each population, the share of its cells contributed by each cohort
#' AFTER equalising every cohort to the same notional cell count.
#'
#' WHY THE EQUALISATION IS THE WHOLE POINT: cohorts differ in sample count and in
#' acquisition depth, so raw contributions are dominated by whichever cohort has
#' the most cells and every population looks like it belongs to that cohort. The
#' normalisation -- the same device used by the cohort-confusion heatmaps
#' elsewhere in the field, which equalise to a fixed cell count per group --
#' removes that, so a row that is not near-uniform means a population is
#' genuinely cohort-skewed.
#'
#' WHY IT IS NOT A SUBSTITUTE FOR group_comparison.png: this pools cells across
#' donors, so it shows a pattern, not a tested effect. A population that looks
#' skewed here is a lead to check against the sample-level abundance test, which
#' is the figure that carries statistics.
#' @param cells Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`, `sample_id`, `population_label` and one column per marker.
#' @param outfile Path to write the figure to.
#' @param group_col Name of the column holding the biological grouping (cohort). Default `"cohort"`.
#' @param pop_col Name of the column holding the population label. Default `"population_label"`.
#' @param norm_to Notional cell count each group is normalised to before shares are taken. Default `1000`.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_cohort_confusion <- function(cells, outfile, group_col = "cohort",
                                 pop_col = "population_label", norm_to = 1000,
                                 panel_label = "", dpi = 300, colors = fcs_colors()) {
  if (is.null(cells) || !nrow(cells)) return(invisible(NULL))
  if (!all(c(group_col, pop_col) %in% names(cells))) return(invisible(NULL))
  g <- as.character(cells[[group_col]]); p <- as.character(cells[[pop_col]])
  ok <- !is.na(g) & nzchar(trimws(g)) & !is.na(p) & nzchar(trimws(p))
  if (!any(ok)) return(invisible(NULL))
  tb <- table(p[ok], g[ok])
  if (nrow(tb) < 1L || ncol(tb) < 2L) {
    log_msg("[fig] fewer than 2 cohorts with labelled cells, confusion heatmap skipped")
    return(invisible(NULL))
  }
  # Equalise cohorts, then read each population as a composition over cohorts.
  eq <- sweep(tb, 2, colSums(tb), "/") * norm_to
  pct <- sweep(eq, 1, rowSums(eq), "/") * 100

  d <- as.data.frame(as.table(pct), stringsAsFactors = FALSE)
  names(d) <- c("population", "group", "pct")
  d <- d[is.finite(d$pct), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))

  glev <- sort(unique(d$group))
  if (!is.null(reference_group()) && reference_group() %in% glev)
    glev <- c(reference_group(), setdiff(glev, reference_group()))
  d$group <- factor(d$group, levels = glev)
  # Order populations by how far they depart from an even split, so the most
  # cohort-skewed populations sit at the top where they will be read first.
  even <- 100 / length(glev)
  skew <- stats::aggregate(pct ~ population, d, function(v) max(abs(v - even)))
  d$population <- factor(d$population, levels = skew$population[order(skew$pct)])

  # Canvas size is decided BEFORE the plot is built, because the text has to be
  # wrapped to a width that is already known. Building first and sizing after is
  # what let the subtitle run off the edge.
  w <- max(5.0, 1.05 * length(glev) + 3.6)
  h <- max(3.6, 0.30 * length(levels(d$population)) + 2.4)

  # Fill limits pinned to the full 0-100%, not left to the data range. Two
  # reasons: a percentage scale that silently rescales to "0-87%" invites the
  # brightest tile to be read as 100, and -- the operative one here -- the text
  # colour below has to reproduce EXACTLY the fill ggplot will draw, which is
  # only possible if the limits are stated rather than inferred.
  vopt <- colors$count_viridis %||% "D"
  d$.txt_col <- contrast_text_colour(d$pct, limits = c(0, 100), option = vopt)

  fig <- ggplot(d, aes(group, population, fill = pct)) +
    geom_tile(colour = colors$tile_border, linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f", pct), colour = .txt_col),
              size = 2.6, fontface = "bold", show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_viridis_c(option = vopt, limits = c(0, 100),
                         name = "% of\npopulation") +
    labs(x = NULL, y = NULL,
         title = paste0("Cohort composition of each population",
                        if (nzchar(panel_label)) paste0(", ", panel_label) else ""),
         subtitle = wrap_plot_text(
           paste0("Each cohort normalised to ", norm_to,
                  " cells before the split is taken, so unequal sample ",
                  "numbers and acquisition depth cannot drive the pattern. ",
                  "An even split is ", sprintf("%.0f", even), "%."), w, pt = 9.5),
         caption = wrap_plot_text(paste0(
           "Descriptive: cells are pooled within a cohort, so this is a pattern, not a ",
           "tested effect.\nConfirm anything that looks skewed against the sample-level ",
           "test in group_comparison_stats.csv."), w, pt = 7)) +
    theme_cyto(colors = colors) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          # Titles anchored to the FIGURE edge, not the panel edge: the y-axis
          # labels here are population names and eat well over an inch, all of
          # which was width the subtitle could not use.
          plot.title.position = "plot",
          plot.subtitle = element_text(lineheight = 1.15),
          plot.caption = element_text(size = 7, hjust = 0, lineheight = 1.2,
                                      colour = colors$caption_text))

  # Wrapping adds lines, so the canvas has to grow to hold them or the figure
  # just crops in a different place.
  h <- h + 0.16 * (lengths(regmatches(fig$labels$subtitle,
                                      gregexpr("\n", fig$labels$subtitle))) +
                   lengths(regmatches(fig$labels$caption,
                                      gregexpr("\n", fig$labels$caption))))
  safe_ggsave(outfile, plot = fig, width = w, height = h, dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile, " (", length(levels(d$population)),
          " populations x ", length(glev), " cohorts)")
  invisible(d)
}

#' Markers each population's gate definition requires to be POSITIVE
#'
#' Feeds the heatmap's gate-audit outline. Reads the same `populations:` spec the
#' scoring uses, so the audit cannot drift from the definitions it audits.
#' `any_of` members are included: any one of them being bright satisfies the gate,
#' so all are legitimate places to look.
#' @param spec Population specification mapping population name to marker directions. See [default_population_spec()].
#' @export
expected_positive_markers <- function(spec) {
  if (!length(spec)) return(NULL)
  out <- lapply(spec, function(rule) {
    if (!is.list(rule)) return(character(0))
    direct <- names(rule)[vapply(rule, function(v)
      is.character(v) && length(v) == 1L && identical(v, "above"), logical(1))]
    anyof <- character(0)
    if (!is.null(rule$any_of) && is.list(rule$any_of))
      anyof <- names(rule$any_of)[vapply(rule$any_of, function(v)
        identical(v, "above"), logical(1))]
    setdiff(unique(c(direct, anyof)), "any_of")
  })
  out[vapply(out, length, integer(1)) > 0L]
}
