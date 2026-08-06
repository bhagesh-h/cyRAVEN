# 4. FIGURES
# =============================================================================

#' Every colour this script draws with, in one place
#'
#' WHY A CONFIG-DRIVEN DEFAULTS LIST RATHER THAN LITERALS SCATTERED THROUGH
#' EACH FIGURE FUNCTION: re-theming used to mean editing R code at a dozen
#' call sites, with no guarantee of catching every one. Every figure function
#' below takes a `colors = fcs_colors()` parameter (default arguments are resolved
#' at CALL time in R, so a run that overrides COLORS via --config, see
#' apply_color_config() in main(), changes what every later figure call
#' picks up with no other code change) and reads exclusively from this list --
#' no hex code or named colour is hardcoded anywhere past this point. Override
#' any subset from --config's `colors:` block; unset keys keep these defaults
#' (see write_config(), which round-trips the resolved list back out as YAML).
#' @export
default_colors <- function() list(
  # Discrete: populations, samples, groups. See pop_palette()'s own docstring
  # for how this specific list was derived and verified.
  # Bright, saturated, one clear hue each, walking the wheel red -> orange ->
  # yellow -> lime -> green -> cyan -> blue -> pink. No black, no violet, no
  # grey. Verified min pairwise CIE Lab distance 28.2 at full opacity, which is
  # well clear of the 12-18 band where colours start being confused; the closest
  # pair is the two greens (#00B34A / #00E0A0), so if a palette entry ever has
  # to be dropped, drop one of those. Violet is deliberately absent even though
  # it would open the blue->pink gap, which is why 9 is the practical ceiling
  # here rather than 12 -- beyond this pop_palette() falls back to hcl.colors().
  population_palette = c("#E8112D","#FF7A00","#FFC800","#8FD400","#00B34A",
                         "#00C2D1","#0072F0","#FF3D9E","#00E0A0"),
  # "Other CD45+" catch-all -- deliberately excluded from population_palette
  # (see population_colours()) so a real population can never be issued it.
  # Not grey any more: it is now hidden from the UMAPs unless --other is passed,
  # and when it IS shown it should read as a real, bright category.
  other_grey       = "#7A5C00",
  # Study/cohort curves in the multigraph overlay: reference first (green), then
  # the patient studies. Fixed so a colour means the same study in every panel.
  study_palette    = c("#00A651", "#E8112D", "#8E44E8"),
  # Sex, fixed everywhere it is drawn.
  sex_palette      = c(male = "#00A651", female = "#E8112D"),
  reference_fill   = "white",   # unfilled bar for the group_comparison reference
  na_fill          = "grey85",  # tile with no value (fig_absolute_counts_qc)

  # Gating diagnostics (recon_diagnostics.png): leukocyte gate box outline and
  # the CD45 threshold's dashed line -- same colour, they're the same kind of
  # "here is the derived cut" annotation. Also the polygon outline in
  # fig_gate_strategy(), which is the same annotation again: a derived boundary
  # drawn over the cells it was derived from.
  gate_highlight   = "#D62728",

  # Learned gate proposals (fig_gate_strategy). Two categories only -- the
  # population being explained, and everything it is being separated from -- so
  # this is deliberately NOT a palette. The non-targets are grey because they
  # are context: colouring both sides equally lets a dense background cloud hide
  # the few hundred cells the gate is actually about.
  gate_target      = "#E8112D",
  gate_nontarget   = "grey75",

  # Threshold review flag (gating_qc.png): a valley-derived cut vs. one that
  # needed a fallback and should be checked (see resolve_threshold()).
  threshold_ok     = "grey20",
  threshold_review = "firebrick",

  # Structural chrome: bar/tile outlines, error bars, significance brackets.
  bar_outline      = "grey20",
  bracket          = "black",
  tile_border      = "white",

  # Theme text and gridlines (theme_cyto()) and figure-level caption text.
  grid_major       = "grey90",
  axis_ticks       = "grey30",
  subtitle_text    = "grey35",
  caption_text     = "grey30",
  label_text       = "grey25",
  empty_panel_text = "grey40",

  # Continuous scales (viridis option names -- see ?viridisLite::viridis).
  density_viridis   = "mako",    # hexbin/density figures: recon scatter, UMAP density
  intensity_viridis = "magma",   # marker-intensity / covariate continuous UMAP colouring
  count_viridis     = "D",       # fig_absolute_counts_qc heatmap (viridis default option)

  # DIVERGING scale, for quantities with a meaningful zero that run both ways:
  # the phenotype heatmap's within-marker z-score, and the batch/threshold-drift
  # diagnostics. Deliberately NOT taken from study_palette -- a colour that means
  # "subject 1" in one figure must not mean "above average" in another. Blue/red
  # about white is the convention for signed heatmaps and keeps the zero point
  # readable as an absence of colour rather than as a hue the eye has to decode.
  heatmap_low      = "#0072F0",
  heatmap_mid      = "#FFFFFF",
  heatmap_high     = "#E8112D"
)

#' Merge a --config `colors:` block over the defaults
#' WHY modifyList() and not replacing wholesale: a config that overrides only
#' `gate_highlight` should not silently blank out every other colour to NULL.
#' @param cfg_colors The cfg colors.
#' @keywords internal
apply_color_config <- function(cfg_colors) {
  if (!length(cfg_colors)) return(default_colors())
  utils::modifyList(default_colors(), cfg_colors)
}

#' Discrete palette that stays distinguishable for many populations
#'
#' WHY the specific ordering: the colours are pre-sorted by greedy
#' farthest-point traversal in CIE Lab space, so the FIRST n entries are a
#' near-maximally separated subset for every n -- a palette that is only
#' well-separated at its full length degrades exactly when a batch happens to
#' yield few populations.
#'
#' WHY THIS LIST REPLACED AN EARLIER 18-COLOUR ONE: that palette was checked
#' for separation only at full opacity, but every UMAP figure draws points at
#' alpha 0.35-0.85 (auto_point_aes) -- 0.35 at the >150k-cell end, which a
#' default-settings embedding of a full cohort routinely reaches -- so points
#' are seen through each other in a dense cluster, and alpha-blending toward
#' the white background compresses colour differences hard. Checked THAT way
#' (each colour rendered at alpha=0.35 over white, then compared in Lab
#' space), the old 18-colour list had five reds/pinks/maroons within 12-18 Lab
#' units of each other (#C81E28/#7F0000/#B03060/#FF6EC7/#D42C88) and two
#' ochre/browns at 12.9 (#8C6218/#A0522D) -- exactly the "these all look the
#' same when they overlap" failure, and brown in particular is the shade a
#' semi-transparent render fades fastest, which is what made it collide with
#' OTHER_GREY too.
#'
#' This list holds exactly ONE colour per hue family with no muddy/low-chroma
#' colour at all. That discipline is deliberate, not automatic: repeatedly
#' adding a "farthest available" next colour to a same-length set kept
#' rebuilding the same failure one family at a time (an added crimson reads
#' fine against everything else individually while still sitting close to a
#' magenta already in the list).
#'
#' WARM-LED, NO DARK SHADES: a pure warm set (red/orange/gold/pink family
#' only, matching how populations read most vividly against a white plot
#' background) tops out at 5-6 mutually distinct colours before pairs start
#' blending on overlap -- warm hues occupy too narrow a slice of colour space
#' for more than that, verified the same way (checked at alpha=0.35, the
#' densest UMAPs actually render at). The first 6 entries below are that warm
#' set; entries 7-9 are BRIGHT cool accents (no navy/black/dark anything, all
#' l >= 55) added only once the warm hues ran out of room, so a study with
#' <=6 populations gets an all-warm palette and larger ones extend into
#' clearly-still-bright, not dark, territory rather than repeating a warm hue
#' that's already taken. Every pair clears 9.2 Lab units at alpha=0.35 (most
#' clear 13+) and every colour stays >= 16.4 units from OTHER_GREY. Beyond 9
#' populations this hands off to hcl.colors(), which spaces hues evenly
#' rather than piling several into one family.
#' @param n Number of items.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @keywords internal
pop_palette <- function(n, colors = fcs_colors()) {
  base <- colors$population_palette
  if (n <= length(base)) return(base[seq_len(n)])
  c(base, grDevices::hcl.colors(n - length(base), "Dark 3"))
}

#' Assign colours to population labels, reserving grey for the catch-all
#'
#' WHY: "Other CD45+" is not a population -- it is everything the gating strategy
#' failed to classify. Colouring it like a lineage invites reading it as one, and
#' when it is large (which is itself the signal to check the thresholds) it
#' dominates the figure in a saturated hue. Grey recedes; that is the point.
#' Fixed meaning-carrying colours, consulted before any palette is issued
#'
#' WHY: a colour that means one group in one figure and a different group in the
#' next is worse than no colour at all -- the reader carries the association
#' across the page whether or not it holds. Sex and study are the two variables
#' drawn repeatedly across these outputs, so both get a fixed assignment
#' (male/healthy green, female/study-1 red, study-2 purple) that every discrete
#' scale checks first. Everything else still falls through to pop_palette().
#'
#' Returns NULL when the levels are not a variable we pin, which is the signal
#' to use the ordinary palette.
#' The study that anchors green is package state, not an argument: it is set once
#' by run_cyraven() from --reference-group (see set_reference_group()) and read
#' here as a DEFAULT ARGUMENT, evaluated at call time. That is what lets every
#' discrete scale colour the cohorts consistently without threading `reference`
#' through a dozen wrapper functions -- the same mechanism fcs_colors() uses.
#' @param levels Character vector of factor levels.
#' @param reference The group every other group is compared against. Default `reference_group()`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @keywords internal
semantic_colours <- function(levels, reference = reference_group(), colors = fcs_colors()) {
  lv <- as.character(levels); if (!length(lv)) return(NULL)
  nm <- tolower(trimws(lv))
  sp <- colors$sex_palette %||% c(male = "#00A651", female = "#E8112D")
  if (all(nm %in% c("male", "female", "m", "f", "w", "man", "woman"))) {
    male <- unname(sp[["male"]]); female <- unname(sp[["female"]])
    map <- c(male = male, m = male, man = male,
             female = female, f = female, w = female, woman = female)
    return(setNames(unname(map[nm]), lv))
  }
  # Study/cohort: the reference group anchors green, the rest take red then
  # purple in sorted order so the assignment is stable across figures and runs.
  if (!is.null(reference) && length(reference) == 1L && reference %in% lv) {
    stp <- colors$study_palette %||% c("#00A651", "#E8112D", "#8E44E8")
    ord <- c(reference, sort(setdiff(lv, reference)))
    return(setNames(rep_len(stp, length(ord)), ord)[lv])
  }
  NULL
}

population_colours <- function(levels, other_pattern = "^Other|unclassified",
                               reference = reference_group(),
                               colors = fcs_colors()) {
  levels <- as.character(levels)
  # Pinned variables (sex, study) win over the palette -- see semantic_colours().
  sem <- semantic_colours(levels, reference = reference, colors = colors)
  if (!is.null(sem)) return(sem)
  is_other <- grepl(other_pattern, levels, ignore.case = TRUE)
  cols <- character(length(levels))
  cols[is_other] <- colors$other_grey
  n_real <- sum(!is_other)
  if (n_real) cols[!is_other] <- pop_palette(n_real, colors = colors)
  setNames(cols, levels)
}

#' Point size / alpha that stay legible across two orders of magnitude of N
#' WHY: a size that reads well at 200k cells is invisible at 3k, and vice versa.
#' Fixed values are the most common cause of unreadable UMAP panels.
#' @param n Number of items.
#' @keywords internal
auto_point_aes <- function(n) {
  if (n <= 5000)        list(size = 0.85, alpha = 0.85)
  else if (n <= 20000)  list(size = 0.55, alpha = 0.70)
  else if (n <= 60000)  list(size = 0.35, alpha = 0.55)
  else if (n <= 150000) list(size = 0.22, alpha = 0.45)
  else                  list(size = 0.14, alpha = 0.35)
}

#' UMAP coloured by a discrete variable
#' @param df The df.
#' @param colour_by The colour by.
#' @param title The title.
#' @param subtitle The subtitle.
#' @param point_size The point size.
#' @param alpha Significance threshold.
#' @param legend_rows The legend rows.
#' @param palette The palette.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @keywords internal
umap_discrete <- function(df, colour_by, title = NULL, subtitle = NULL,
                          point_size = NULL, alpha = NULL, legend_rows = NULL,
                          palette = NULL, colors = fcs_colors()) {
  d <- df[!is.na(df[[colour_by]]), ]
  .a <- auto_point_aes(nrow(d))
  if (is.null(point_size)) point_size <- .a$size
  if (is.null(alpha))      alpha      <- .a$alpha
  lv <- sort(unique(as.character(d[[colour_by]])))
  # Put the catch-all last so it is drawn UNDER the real populations and sits at
  # the end of the legend: a large unclassified fraction must not paint over the
  # lineages it failed to classify.
  is_other <- grepl("^Other|unclassified", lv, ignore.case = TRUE)
  lv <- c(lv[!is_other], lv[is_other])
  d[[colour_by]] <- factor(as.character(d[[colour_by]]), levels = lv)
  cols <- if (is.null(palette)) population_colours(lv, colors = colors) else palette
  # Draw order: catch-all first (bottom), then real populations over it.
  d <- d[order(match(as.character(d[[colour_by]]), rev(lv))), , drop = FALSE]
  # Discrete colouring needs MORE opacity than a continuous scale: hue identity
  # is what the reader is matching against the legend, and a semi-transparent
  # point over white shifts its apparent hue toward the background. Legend keys
  # are drawn fully opaque and large enough to compare against the plot.
  fig <- ggplot(d, aes(umap_1, umap_2, colour = .data[[colour_by]])) +
    geom_point(size = point_size * 1.15, alpha = min(1, alpha + 0.25), stroke = 0) +
    scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1),
                                 nrow = legend_rows)) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    coord_equal() + theme_cyto(colors = colors)
  fig
}

#' UMAP coloured by a continuous variable (marker intensity or covariate)
#' @param df The df.
#' @param colour_by The colour by.
#' @param title The title.
#' @param subtitle The subtitle.
#' @param point_size The point size.
#' @param alpha Significance threshold.
#' @param limits Length-2 numeric vector giving the scale limits.
#' @param legend_name The legend name. Default `"asinh"`.
#' @param compact_bar The compact bar. Default `FALSE`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @keywords internal
umap_continuous <- function(df, colour_by, title = NULL, subtitle = NULL,
                            point_size = NULL, alpha = NULL, limits = NULL,
                            legend_name = "asinh", compact_bar = FALSE,
                            colors = fcs_colors()) {
  d <- df[is.finite(df[[colour_by]]), ]
  .a <- auto_point_aes(nrow(d))
  if (is.null(point_size)) point_size <- .a$size
  if (is.null(alpha))      alpha      <- .a$alpha
  if (is.null(limits)) limits <- quantile(d[[colour_by]], c(0.01, 0.99), na.rm = TRUE)
  # compact_bar: a short horizontal bar labelled only with its two end values,
  # for use inside a small multi-panel grid where a full legend would not fit.
  # Built here rather than layered on afterwards, so only ONE colour scale is
  # ever added to the plot.
  sc <- if (compact_bar) {
    scale_colour_viridis_c(
      option = colors$intensity_viridis, limits = limits, oob = scales::squish, name = NULL,
      breaks = as.numeric(limits), labels = sprintf("%.1f", as.numeric(limits)),
      guide = guide_colourbar(barwidth = unit(2.6, "lines"),
                              barheight = unit(0.35, "lines"),
                              ticks = FALSE, direction = "horizontal",
                              label.position = "bottom"))
  } else {
    scale_colour_viridis_c(option = colors$intensity_viridis, limits = limits,
                           oob = scales::squish, name = legend_name)
  }
  fig <- ggplot(d, aes(umap_1, umap_2, colour = .data[[colour_by]])) +
    geom_point(size = point_size, alpha = alpha, stroke = 0) +
    sc +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    coord_equal() + theme_cyto(colors = colors)
  fig
}

#' Population + sample + covariate colouring panel set
#'
#' WHAT: writes the main UMAP figure (population identity, sample, and any
#'       patient covariates present) for one panel group.
#' WHY:  one shared embedding viewed many ways is what makes populations
#'       comparable across patients; separate embeddings per sample would give
#'       coordinates that cannot be compared.
#' @param cells data.frame with umap_1, umap_2, population_label, sample_id and
#'   optionally covariate columns.
#' @param covariates character vector of covariate column names to render if
#'   present and non-constant.
#' @param outfile Path to write the figure to.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param feature_cols The feature cols. Default `character(0)`.
#' @param width Figure width in inches. Default `15`.
#' @param height Figure height in inches. Default `11`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @return the composed patchwork object (also written to `outfile`).
#' @param covariates NULL (default) to DISCOVER plottable covariates from the
#'   joined table, or an explicit character vector to fix the set.
#'
#'   WHY discovery is the default: a fixed list can only plot the covariates
#'   someone anticipated. The single most important variable in a case-control
#'   study is the group column -- cohort, genotype, treatment arm -- and its name
#'   is study-specific, so a hardcoded list silently omits exactly the panel the
#'   experiment was designed to produce. Anything that arrives from the sample map
#'   or patient table and varies across cells is plottable; the exclusion list
#'   below is structural (coordinates, gate bookkeeping, marker intensities),
#'   because those are plotted by other figures or are not covariates at all.
#' @export
fig_umap_overview <- function(cells, outfile, panel_label = "",
                             covariates = NULL, feature_cols = character(0),
                             width = 15, height = 11, dpi = 300, colors = fcs_colors()) {
  if (is.null(covariates)) {
    structural <- c("umap_1", "umap_2", "population_label", "sample_id", "panel",
                    "event_index", "file", "is_control", "patient_id_matched",
                    feature_cols)
    cand <- setdiff(names(cells), structural)
    # Keep only columns that vary and are not near-unique free text (a per-cell
    # identifier or a comment column produces a legend with one key per sample).
    covariates <- cand[vapply(cand, function(cv) {
      v <- cells[[cv]]
      u <- unique(v[!is.na(v)])
      length(u) >= 2L && (is.numeric(v) || length(u) <= 12L)
    }, logical(1))]
    # patient_id last: useful but usually the busiest legend, and it is exempt
    # from the level cap above (many patients is expected, not a sign of free text).
    if ("patient_id" %in% names(cells) &&
        length(unique(cells$patient_id[!is.na(cells$patient_id)])) >= 2L)
      covariates <- c(setdiff(covariates, "patient_id"), "patient_id")
    covariates <- unique(covariates)
    log_msg("  covariate panels: ",
            if (length(covariates)) paste(covariates, collapse = ", ") else "none")
  }
  plots <- list()
  if ("population_label" %in% names(cells)) {
    npop <- length(unique(cells$population_label))
    plots[["pop"]] <- umap_discrete(
      cells, "population_label",
      title = "Gated population",
      subtitle = paste0(format(nrow(cells), big.mark = ","), " CD45+ cells, ",
                        npop, " labels"),
      legend_rows = ceiling(npop / 3), colors = colors)
  }
  # A one-level discrete panel carries no information (every point one colour) and
  # wastes half the figure. Batches with a single stained sample per panel are
  # common, so skip rather than render it.
  if ("sample_id" %in% names(cells) && length(unique(cells$sample_id)) > 1L) {
    plots[["sample"]] <- umap_discrete(cells, "sample_id", title = "Sample",
                                       subtitle = "equal cells per sample", colors = colors)
  }
  for (cv in covariates) {
    if (!cv %in% names(cells)) next
    v <- cells[[cv]]
    if (all(is.na(v)) || length(unique(v[!is.na(v)])) < 2L) next
    plots[[cv]] <- if (is.numeric(v)) {
      umap_continuous(cells, cv, title = pretty_label(cv), legend_name = pretty_label(cv),
                      colors = colors)
    } else {
      umap_discrete(cells, cv, title = pretty_label(cv), colors = colors)
    }
  }
  if (!length(plots)) { warning("[fig] nothing to plot"); return(invisible(NULL)) }
  fig <- patchwork::wrap_plots(plots, ncol = min(2L, length(plots))) +
    patchwork::plot_annotation(
      title = paste0("Shared CD45+ UMAP", if (nzchar(panel_label)) paste0(", ", panel_label) else ""),
      theme = theme(plot.title = element_text(face = "bold", size = 14)))
  safe_ggsave(outfile, plot = fig, width = width, height = height, dpi = dpi,
             limitsize = FALSE)
  message("[fig] wrote ", outfile)
  invisible(fig)
}

#' Combined + per-group UMAP small multiples
#'
#' WHAT: for each colouring fig_umap_overview() already draws (gated
#' population, sample, and discovered covariates -- excluding the group column
#' itself), one FACETED ROW: "All samples" first (every group pooled -- the
#' same view fig_umap_overview() shows), then one facet per resolved group,
#' each restricted to that group's cells on the same shared axes.
#'
#' WHY A SEPARATE FIGURE, NOT A CHANGE TO fig_umap_overview(): that figure
#' answers "what does the pooled embedding look like"; this one answers "does
#' any group's embedding look different from the others" -- the same split
#' section 8.1 draws between population_frequencies.png and group_comparison.png,
#' applied to the embedding itself instead of to abundance.
#'
#' WHY facet_wrap() PER ROW, NOT ONE PANEL PER (colouring x group) HAND-BUILT
#' WITH patchwork: an earlier version built every small panel as its own
#' ggplot object with its own manually-sized title and legend, then tiled
#' them with patchwork::wrap_plots(). It looked fine on short synthetic group
#' names but collapsed into an unreadable mess on this study's real group
#' names (a long cohort label, say) -- long strip text has nowhere
#' fixed-size manual titles can absorb it. facet_wrap() is the same
#' mechanism umap_density_by_group.png already uses successfully: ONE ggplot
#' object per row, ONE colour scale, ONE legend, and ggplot's own (far more
#' battle-tested) strip-label layout instead of a hand-rolled one.
#'
#' WHY the "All samples" facet is a literal duplicate of every group's rows,
#' not a separate panel bolted on: a single shared scale_colour_*() only
#' exists once you give ggplot one data frame and one aesthetic mapping to
#' compute it from. Duplicating the pooled rows under a synthetic group level
#' means population_colours()/scale_colour_viridis_c() sees every level that
#' exists ANYWHERE, so a population missing from one group keeps its colour
#' rather than shifting every alphabetically-later population in that facet.
#' facet_wrap()'s default `scales = "fixed"` gives every facet the same x/y
#' range for the same reason -- coordinates only mean the same thing across
#' panels if the axes are literally the same.
#'
#' @param cells data.frame with umap_1, umap_2, population_label, sample_id,
#'   and the resolved group column.
#' @param group_col column in `cells` defining the groups (e.g. "cohort").
#'   Returns NULL (writes nothing) when it resolves to fewer than 2 groups --
#'   there is nothing to compare side by side.
#' @param outfile Path to write the figure to.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param covariates Column names in the patient table to screen as confounders.
#' @param feature_cols The feature cols. Default `character(0)`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_umap_overview_by_group <- function(cells, outfile, group_col,
                                       panel_label = "", covariates = NULL,
                                       feature_cols = character(0), dpi = 200,
                                       colors = fcs_colors()) {
  if (!group_col %in% names(cells)) return(invisible(NULL))
  d0 <- cells[!is.na(cells[[group_col]]) & nzchar(as.character(cells[[group_col]])), ,
             drop = FALSE]
  groups <- sort(unique(as.character(d0[[group_col]])))
  if (length(groups) < 2L) return(invisible(NULL))
  n_facet <- length(groups) + 1L

  d_all <- d0; d_all[[group_col]] <- "All samples"
  dd <- rbind(d_all, d0)
  dd[[group_col]] <- factor(dd[[group_col]], levels = c("All samples", groups))
  rm(d_all)

  if (is.null(covariates)) {
    structural <- c("umap_1", "umap_2", "population_label", "sample_id", "panel",
                    "event_index", "file", "is_control", "patient_id_matched",
                    group_col, feature_cols)
    cand <- setdiff(names(d0), structural)
    covariates <- cand[vapply(cand, function(cv) {
      v <- d0[[cv]]
      u <- unique(v[!is.na(v)])
      length(u) >= 2L && (is.numeric(v) || length(u) <= 12L)
    }, logical(1))]
    if ("patient_id" %in% names(d0) &&
        length(unique(d0$patient_id[!is.na(d0$patient_id)])) >= 2L)
      covariates <- c(setdiff(covariates, "patient_id"), "patient_id")
    covariates <- unique(covariates)
  }

  # Long group names wrap onto 2-3 lines
  # inside the facet strip instead of overflowing into neighbouring panels.
  # Appended to whatever umap_discrete()/umap_continuous() already builds --
  # facet is orthogonal to coord/scale/geom, so `+`-ing it on top is safe.
  facet_row <- list(
    facet_wrap(vars(.data[[group_col]]), nrow = 1,
              labeller = label_wrap_gen(width = 16)),
    theme(plot.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 7.5)))

  row_plots <- list()
  add_row <- function(key, g) row_plots[[key]] <<- g

  # umap_discrete()/umap_continuous() are reused rather than reimplemented so
  # this row gets everything they already do correctly for a single panel --
  # NA filtering, point size/alpha for the actual N, and (population_label
  # specifically) drawing "Other CD45+" underneath the real populations and
  # listing it last in the legend, not sorted in wherever "O" falls
  # alphabetically. Passing them the FULL `dd` (every facet's rows at once,
  # "All samples" included) is what makes the colour scale -- computed once,
  # inside them, from dd's complete level set -- identical in every facet.
  if ("population_label" %in% names(dd))
    add_row("Gated population",
      umap_discrete(dd, "population_label", title = "Gated population", colors = colors) + facet_row)

  if ("sample_id" %in% names(dd) && length(unique(d0$sample_id)) > 1L)
    # No legend on this row: sample identity is local to whichever group a
    # sample belongs to (colours are assigned once across ALL samples so
    # they don't collide, but a 20+ entry legend for "which sample is which"
    # is not the question this row answers). It answers "does one sample
    # dominate a cluster inside its own group" -- visible without one.
    add_row("Sample",
      umap_discrete(dd, "sample_id", title = "Sample", colors = colors) + facet_row +
        theme(legend.position = "none"))

  for (cv in covariates) {
    if (!cv %in% names(dd)) next
    v <- dd[[cv]]
    if (all(is.na(v)) || length(unique(v[!is.na(v)])) < 2L) next
    add_row(pretty_label(cv),
      (if (is.numeric(v))
         umap_continuous(dd, cv, title = pretty_label(cv), legend_name = NULL, colors = colors)
       else umap_discrete(dd, cv, title = pretty_label(cv), colors = colors)) + facet_row)
  }
  if (!length(row_plots)) return(invisible(NULL))

  n_by_group <- table(dd[[group_col]])
  fig <- patchwork::wrap_plots(row_plots, ncol = 1) +
    patchwork::plot_annotation(
      title = paste0("Shared CD45+ UMAP: combined vs per-group",
                     if (nzchar(panel_label)) paste0(", ", panel_label) else ""),
      subtitle = paste(
        "First facet in every row pools all groups; the rest are the SAME",
        "shared embedding restricted to one group's cells, a cluster's",
        "position is directly comparable across every facet."),
      caption = paste0("n = ", paste(sprintf("%s: %s", names(n_by_group),
                       format(as.integer(n_by_group), big.mark = ",")),
                       collapse = ", "), " cells"),
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(size = 8),
                    plot.caption = element_text(size = 7, hjust = 0, colour = colors$caption_text)))
  safe_ggsave(outfile, plot = fig, width = max(7, 2.9 * n_facet + 2.2),
             height = 3.3 * length(row_plots) + 1.4, dpi = dpi, limitsize = FALSE)
  message("[fig] wrote ", outfile, " (", length(row_plots), " colouring(s) x ",
          n_facet, " facets)")
  invisible(fig)
}

#' Turn a snake_case column name into a readable axis/panel label
#' @param x A vector of values.
#' @keywords internal
pretty_label <- function(x) {
  x <- gsub("_", " ", x)
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}

#' Marker-expression grid over the shared embedding
#'
#' WHAT: one small UMAP per marker, coloured by that marker's asinh intensity.
#' WHY:  this is the check that the embedding is biologically structured -- each
#'       lineage marker should light up a coherent region. If markers are
#'       smeared uniformly, the embedding is driven by noise (which is exactly
#'       what the ungated template produced).
#' @param cells data.frame with umap_1/umap_2 plus one column per marker.
#' @param markers character vector of marker column names to render.
#' @param outfile Path to write the figure to.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param ncol Number of panel columns; NULL computes one that keeps the canvas roughly square.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param panel_size The panel size. Default `2.5`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_marker_grid <- function(cells, markers, outfile, panel_label = "",
                            ncol = NULL, dpi = 300, panel_size = 2.5,
                            colors = fcs_colors()) {
  markers <- markers[markers %in% names(cells)]
  if (!length(markers)) { warning("[fig] no markers to grid"); return(invisible(NULL)) }
  ncol <- if (is.null(ncol)) ceiling(sqrt(length(markers))) else ncol
  nrow <- ceiling(length(markers) / ncol)
  # grid panels are small, so shrink one notch below the standalone size
  ps <- auto_point_aes(nrow(cells))$size * 0.7
  # Each panel carries its OWN colourbar, labelled with that marker's actual
  # low/high asinh values.
  #
  # WHY per-panel rather than one shared bar: the panels do NOT share a colour
  # scale -- each is independently clipped to its own 1st-99th percentile, because
  # markers differ by orders of magnitude in brightness and a common scale would
  # flatten every dim marker to uniform black. A single figure-level colourbar
  # would therefore be actively wrong: it would imply comparability of colour
  # across panels that does not exist. Numeric end labels make the per-panel
  # range explicit, so "bright yellow" is read as bright *for this marker*.
  plots <- lapply(markers, function(m) {
    umap_continuous(cells, m, title = m, point_size = ps, compact_bar = TRUE,
                    colors = colors) +
      theme(legend.position = "bottom",
            legend.margin = margin(t = -4, b = 0),
            legend.text = element_text(size = 6),
            axis.title = element_blank(), axis.text = element_blank(),
            plot.title = element_text(size = 9))
  })
  fig <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = paste0("Marker expression on the shared embedding",
                     if (nzchar(panel_label)) paste0(", ", panel_label) else ""),
      subtitle = paste("asinh-transformed intensity; each panel has its OWN",
                       "colour scale, clipped to its 1st-99th percentile.",
                       "\nBar labels give that marker's low/high value,",
                       "colours are NOT comparable between panels."),
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(size = 9)))
  safe_ggsave(outfile, plot = fig, width = ncol * panel_size + 0.5,
             height = nrow * (panel_size + 0.28) + 1.2, dpi = dpi, limitsize = FALSE)
  message("[fig] wrote ", outfile, " (", length(markers), " markers)")
  invisible(fig)
}

#' Per-sample (or per-group) density comparison over the shared embedding
#'
#' WHY: side-by-side 2D density of the SAME embedding shows where each sample
#'      or group gains or loses cells -- the comparison the shared embedding
#'      exists for. Called with `facet_by = "sample_id"` (the default) for
#'      per-sample QC, and again with the resolved group column when 2+ groups
#'      exist, so cohorts can be compared the same way sample_id already is.
#' @param facet_by any column present in `cells` -- sample_id, cohort, or any
#'   other covariate the patient table carried through.
#' @param cells Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`, `sample_id`, `population_label` and one column per marker.
#' @param outfile Path to write the figure to.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_density_by_sample <- function(cells, outfile, panel_label = "",
                                 facet_by = "sample_id", dpi = 300,
                                 colors = fcs_colors()) {
  if (!facet_by %in% names(cells)) return(invisible(NULL))
  # A cell with no value for facet_by (e.g. an unmatched patient_id) would
  # otherwise draw its own "NA" panel, which isn't a group and would confuse
  # the comparison this figure exists to support.
  d <- cells[!is.na(cells[[facet_by]]) & nzchar(as.character(cells[[facet_by]])), ,
            drop = FALSE]
  n <- length(unique(d[[facet_by]]))
  if (!n) return(invisible(NULL))
  # One row reads best for a HANDFUL of facets -- a few groups, a modest
  # sample count -- which is the common case this was written for. Beyond
  # max_row it would run off the page long before it ran out of readable
  # width, so wrap into a roughly square grid instead: width and height both
  # grow with sqrt(n) rather than one of them growing linearly without bound
  # (a 300-sample cohort in one row would be a ~1000-inch-wide image).
  max_row <- 12L
  ncol <- if (n <= max_row) n else ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  fig <- ggplot(d, aes(umap_1, umap_2)) +
    stat_bin_hex(bins = 70, aes(fill = after_stat(count))) +
    scale_fill_viridis_c(option = colors$density_viridis, trans = "log10", name = "cells",
                         labels = scales::label_number(accuracy = 1)) +
    facet_wrap(vars(.data[[facet_by]]), ncol = ncol) +
    labs(title = paste0("Cell density per ", pretty_label(facet_by),
                        if (nzchar(panel_label)) paste0(", ", panel_label) else ""),
         subtitle = "shared embedding, log-scaled counts",
         x = "UMAP 1", y = "UMAP 2") +
    coord_equal() + theme_cyto(colors = colors)
  safe_ggsave(outfile, plot = fig, width = max(5, 3.6 * ncol),
             height = max(4.2, 3.3 * nrow + 0.8), dpi = dpi, limitsize = FALSE)
  message("[fig] wrote ", outfile)
  invisible(fig)
}

#' Initial reconnaissance / QC diagnostic figure -- run on EVERY batch
#'
#' WHAT: a three-row diagnostic across all input files:
#'   row 1  log10 FSC-A vs log10 SSC-A occupancy with the derived leukocyte gate
#'   row 2  CD45 (asinh) vs log10 SSC-A with the derived CD45 cutoff
#'   row 3  overlaid density of the lineage markers within the CD45+ parent
#' WHY:  this is the figure that catches the three failure modes that silently
#'       ruin a run: (a) gating debris instead of leukocytes -- visible as a
#'       dense low-FSC blob sitting at CD45 background, (b) a mis-placed CD45
#'       cutoff, and (c) unstained or failed-staining files, visible as marker
#'       densities collapsed into a single background peak. It must be inspected
#'       before any downstream result is trusted.
#'
#' @param recon list of per-file diagnostic records, each a list with:
#'   sample_id, panel, fsc_log10, ssc_log10 (numeric vectors, may be subsampled),
#'   cd45 (asinh vector or NULL), gate (list with fsc_lo/fsc_hi/ssc_lo/ssc_hi),
#'   cd45_threshold (numeric or NA), marker_densities (named list of asinh
#'   vectors within the CD45+ parent), verdict (character, from staining QC).
#' @param outfile Path to write the figure to.
#' @param max_points The max points. Default `40000L`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param dens_markers The dens markers.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @return the composed patchwork object (also written to `outfile`).
#' @export
fig_recon_diagnostics <- function(recon, outfile, max_points = 40000L,
                                  dpi = 300, dens_markers = NULL, colors = fcs_colors()) {
  if (!length(recon)) { warning("[fig] no recon records"); return(invisible(NULL)) }
  # Fixed seed so the thinning below picks the same cells on a re-run, and the
  # stream put back afterwards so drawing a figure cannot change which cells the
  # embedding samples. Same guard, and the same reason, as fig_gating_qc().
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(1L)
  thin <- function(x, n = max_points)
    if (length(x) > n) x[sort(sample.int(length(x), n))] else x

  row1 <- row2 <- row3 <- list()
  for (r in recon) {
    lab <- paste0(r$sample_id, if (!is.null(r$panel)) paste0("  [", r$panel, "]") else "")
    # verdict strings run to ~70 chars; unwrapped they overrun the panel width
    # and collide with the neighbouring column's title.
    sub <- if (!is.null(r$verdict))
      paste(strwrap(r$verdict, width = 38), collapse = "\n") else NULL
    idx <- if (length(r$fsc_log10) > max_points)
      sort(sample.int(length(r$fsc_log10), max_points)) else seq_along(r$fsc_log10)

    # --- row 1: scatter occupancy + derived leukocyte gate
    d1 <- data.frame(fsc = r$fsc_log10[idx], ssc = r$ssc_log10[idx])
    p1 <- ggplot(d1, aes(fsc, ssc)) +
      stat_bin_hex(bins = 80, aes(fill = after_stat(count))) +
      scale_fill_viridis_c(option = colors$density_viridis, trans = "log10", guide = "none") +
      labs(title = lab, subtitle = sub,
           x = expression(log[10]~"FSC-A"), y = expression(log[10]~"SSC-A")) +
      theme_cyto(10, colors = colors)
    g <- r$gate
    if (!is.null(g))
      p1 <- p1 + annotate("rect", xmin = g$fsc_lo, xmax = g$fsc_hi,
                          ymin = g$ssc_lo, ymax = g$ssc_hi,
                          fill = NA, colour = colors$gate_highlight, linewidth = 0.7)
    row1[[length(row1) + 1L]] <- p1

    # --- row 2: CD45 vs SSC with the derived cutoff
    if (!is.null(r$cd45) && length(r$cd45)) {
      d2 <- data.frame(cd45 = r$cd45[idx], ssc = r$ssc_log10[idx])
      p2 <- ggplot(d2, aes(cd45, ssc)) +
        stat_bin_hex(bins = 80, aes(fill = after_stat(count))) +
        scale_fill_viridis_c(option = colors$density_viridis, trans = "log10", guide = "none") +
        labs(x = "CD45 (asinh)", y = expression(log[10]~"SSC-A")) +
        theme_cyto(10, colors = colors)
      if (is.finite(r$cd45_threshold %||% NA))
        p2 <- p2 + geom_vline(xintercept = r$cd45_threshold, colour = colors$gate_highlight,
                              linetype = "22", linewidth = 0.7)
      row2[[length(row2) + 1L]] <- p2
    } else {
      row2[[length(row2) + 1L]] <- empty_panel("CD45 not in panel", colors = colors)
    }

    # --- row 3: lineage marker densities within CD45+
    md <- r$marker_densities
    if (!is.null(dens_markers)) md <- md[intersect(dens_markers, names(md))]
    # The guard is on the ASSEMBLED frame, not on md. md can hold an entry per
    # lineage marker and still yield nothing to draw -- every entry drops out
    # when the parent gate selected no cells, or when a marker is present but
    # all-NA -- and rbind over a list of NULLs returns NULL, not a zero-row
    # frame. Testing length(md) alone then hands ggplot a NULL data argument and
    # the run dies several frames later with "object 'value' not found", which
    # names an aesthetic and points nowhere near the empty gate that caused it.
    d3 <- if (length(md)) do.call(rbind, lapply(names(md), function(m) {
        v <- md[[m]]; v <- v[is.finite(v)]
        if (!length(v)) return(NULL)
        data.frame(marker = m, value = thin(v, 20000L))
      })) else NULL
    if (!is.null(d3) && nrow(d3)) {
      p3 <- ggplot(d3, aes(value, colour = marker)) +
        geom_density(linewidth = 0.5, adjust = 1.2) +
        scale_colour_manual(values = pop_palette(length(unique(d3$marker)), colors = colors),
                            name = NULL) +
        guides(colour = guide_legend(ncol = 2, override.aes = list(linewidth = 1.2))) +
        labs(x = "asinh intensity", y = "density",
             # For a control the "CD45+" parent is a quantile-fallback slice, so
             # any apparent bimodality here is the cut, not real staining. Say so
             # on the panel rather than letting the reader infer staining.
             caption = if (identical(r$cd45_source, "quantile_fallback"))
               "parent = top-N% fallback slice, not a real CD45+ gate" else NULL) +
        theme_cyto(10, colors = colors) +
        theme(legend.text = element_text(size = 6.5),
              legend.key.size = unit(0.6, "lines"),
              plot.caption = element_text(size = 6, colour = colors$caption_text,
                                          hjust = 0))
      row3[[length(row3) + 1L]] <- p3
    } else {
      row3[[length(row3) + 1L]] <- empty_panel("no CD45+ cells to profile", colors = colors)
    }
  }

  # One sample per COLUMN, 3 rows, was fine at cohort sizes this was
  # developed against but scales the canvas width linearly and without bound
  # (a few hundred samples in one row is a page many yards wide). Past
  # max_col samples, wrap into BANDS instead: each band is its own 3-row x
  # max_col grid, stacked vertically. patchwork's byrow fill is sequential
  # and has no notion of "band", so a short final band is padded with blank
  # panels up to max_col -- otherwise its few real panels would bleed into
  # what should be the next row's cells.
  n <- length(recon)
  max_col <- 20L
  n_bands <- ceiling(n / max_col)
  if (n_bands <= 1L) {
    # Single band, no padding needed: the grid is exactly n columns wide, not
    # max_col -- padding to max_col here while declaring only n columns is
    # exactly the mismatch that used to crash wrap_plots() below whenever a
    # run had fewer than max_col samples.
    ncol_grid <- n
    panels_ordered <- c(row1, row2, row3)
  } else {
    ncol_grid <- max_col
    pad_to <- n_bands * max_col
    pad_row <- function(lst) c(lst, if (pad_to > length(lst))
      replicate(pad_to - length(lst), patchwork::plot_spacer(), simplify = FALSE))
    row1 <- pad_row(row1); row2 <- pad_row(row2); row3 <- pad_row(row3)
    panels_ordered <- do.call(c, lapply(seq_len(n_bands), function(b) {
      idx <- ((b - 1L) * max_col + 1L):(b * max_col)
      c(row1[idx], row2[idx], row3[idx])
    }))
  }
  fig <- patchwork::wrap_plots(panels_ordered, ncol = ncol_grid, nrow = 3 * n_bands,
                               byrow = TRUE) +
    patchwork::plot_annotation(
      title = "Initial QC diagnostics, inspect before trusting any downstream result",
      subtitle = paste("Row 1: scatter occupancy with the derived leukocyte gate (outlined).",
                       "Row 2: CD45 vs SSC-A with the derived cutoff (dashed).",
                       "Row 3: lineage-marker densities within CD45+.",
                       if (n_bands > 1)
                         paste0("\n", n, " samples wrapped into ", n_bands,
                               " bands of up to ", max_col, ".") else NULL),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 9, colour = colors$subtitle_text)))
  safe_ggsave(outfile, plot = fig, width = max(9, 4.2 * ncol_grid),
             height = 11.5 * n_bands, dpi = dpi, limitsize = FALSE)
  message("[fig] wrote ", outfile, " (", n, " sample(s)",
          if (n_bands > 1) paste0(" in ", n_bands, " band(s)") else "", ")")
  invisible(fig)
}

#' Placeholder panel used when a diagnostic cannot be drawn
#' @param msg The msg.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @keywords internal
empty_panel <- function(msg, colors = fcs_colors()) {
  ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 3, colour = colors$empty_panel_text) +
    theme_void()
}

#' Per-gate QC figure: every applied gate shown per sample with its threshold
#'
#' WHY THIS IS A DELIVERABLE: a gated analysis is only as trustworthy as its
#' gates, and a threshold that is numerically plausible can still sit in the
#' wrong place. Showing every cut on the distribution it was derived from is the
#' only way a reviewer can check the gating without re-running anything.
#' DO NOT "SIMPLIFY" THIS BACK TO geom_density() ON RAW CELLS. It was written
#' that way and it OOM-killed the container mid-run.
#'
#' The obvious implementation hands every CD45+ cell to geom_density() and lets
#' ggplot compute the curve. That builds one row per CELL per MARKER: with 25
#' samples x 13 markers x ~813k CD45+ cells (--max-events-per-file 300000) it is
#' 10,572,757 rows through do.call(rbind, ...), on top of the ~1.9 GB of
#' reads$exprs + pops$tmat that main() is still holding at this point. Measured
#' under --memory=6g: SIGKILL, exit 137. The failure is SILENT -- the OOM killer
#' leaves no R error, no traceback, nothing in the log. The run simply stops
#' after STEP 4 with recon_diagnostics.png written and gating_qc.png missing.
#' It survived at --max-events-per-file 120000 purely by luck of scale.
#'
#' The cure is NOT to subsample the cells. Measured on this study's own data,
#' thinning to 5,000 cells per panel distorts the drawn curve by up to 21% of
#' its peak height and makes a 2%-of-CD45+ population -- a perfectly ordinary
#' size for Vd1/Vd2, NKT or dendritic cells here -- invisible in 19 of 20 draws.
#' A gating-QC figure that silently drops the rare populations whose gates most
#' need checking is worse than no figure.
#'
#' So the curve is computed here, from EVERY cell, and only the curve is passed
#' to ggplot: stats::density() reduces each (sample, marker) to 512 points, which
#' is what geom_density would have drawn anyway (same default bw.nrd0 bandwidth,
#' same 512-point grid). 10.6M rows become 325 x 512 = 166k, a ~64x reduction,
#' with NO loss of fidelity and no subsampling. It is also RNG-free, so unlike a
#' sampling-based fix it cannot shift the .Random.seed stream that STEP 6's UMAP
#' cell selection (sample(), line ~3559) draws from -- the embedding is
#' bit-identical to before this change.
#' @param recon The recon.
#' @param outfile Path to write the figure to.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_gating_qc <- function(recon, outfile, dpi = 200, colors = fcs_colors()) {
  curves <- list(); labs <- list()
  for (r in recon) {
    thr <- r$thresholds %||% list()
    for (mk in names(thr)) {
      x <- r$marker_densities[[mk]]
      if (is.null(x) || !length(x)) next
      x <- as.numeric(x); x <- x[is.finite(x)]
      # density() cannot pick a bandwidth from a single point.
      if (length(x) < 2L) next
      # from/to are pinned to the data range on purpose, and n = 512 matches
      # stat_density's default. ggplot evaluates its KDE over the panel's data
      # range, whereas stats::density() defaults to range +/- 3*bandwidth; left
      # at the default the curve would be drawn on a wider grid than
      # geom_density used and every panel's tails would shift. With them pinned
      # the rendered layer is IDENTICAL to the old geom_density output to
      # floating-point equality (verified against ggplot_build on this study's
      # own data), so this is purely a memory fix and changes no figure.
      dk <- stats::density(x, n = 512, from = min(x), to = max(x))
      curves[[length(curves) + 1L]] <- data.frame(
        sample_id = r$sample_id, marker = mk,
        value = dk$x, density = dk$y, stringsAsFactors = FALSE)
      labs[[length(labs) + 1L]] <- data.frame(
        sample_id = r$sample_id, marker = mk,
        threshold = thr[[mk]]$threshold %||% NA_real_,
        source = thr[[mk]]$source %||% "unknown",
        needs_review = isTRUE(thr[[mk]]$needs_review),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(curves)) { log_msg("[fig] no gate QC data, skipped"); return(invisible(NULL)) }
  d <- do.call(rbind, curves)
  lab <- do.call(rbind, labs)
  lab$tag <- ifelse(lab$needs_review, paste0(lab$source, " (REVIEW)"), lab$source)
  fig <- ggplot(d, aes(value, density)) +
    geom_area(aes(fill = sample_id), colour = NA, alpha = 0.45) +
    geom_vline(data = lab, aes(xintercept = threshold, colour = needs_review),
               linewidth = 0.5, linetype = "dashed", show.legend = FALSE) +
    geom_text(data = lab, aes(x = threshold, y = Inf, label = tag),
              hjust = -0.05, vjust = 1.6, size = 2.4, colour = colors$label_text) +
    scale_colour_manual(values = c(`FALSE` = colors$threshold_ok, `TRUE` = colors$threshold_review)) +
    scale_fill_manual(values = pop_palette(length(unique(d$sample_id)), colors = colors), name = NULL) +
    facet_grid(sample_id ~ marker, scales = "free") +
    labs(title = "Gating QC, every applied threshold on its own distribution",
         subtitle = paste("dashed line = applied cut, labelled with its source;",
                          "highlighted colour = needs review (quantile fallback)"),
         x = "asinh-transformed intensity", y = "density") +
    theme_cyto(9, colors = colors) + theme(legend.position = "none",
                          strip.text.y = element_text(angle = 0, size = 7),
                          strip.text.x = element_text(size = 7),
                          axis.text = element_text(size = 6))
  nmk <- length(unique(d$marker)); nsm <- length(unique(d$sample_id))
  # facet_grid(sample_id ~ marker) is a true 2D cross-tabulation -- unlike
  # facet_wrap it cannot be wrapped into a bounded grid without breaking the
  # marker-column alignment across samples. safe_ggsave's dpi clamp keeps a
  # large sample x marker count from crashing the raster device; it gets
  # dense rather than unbounded, which is acceptable for a figure meant to be
  # read per sample rather than as a single glance.
  safe_ggsave(outfile, plot = fig, width = max(8, 1.9 * nmk), height = max(4, 1.7 * nsm),
             dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' Split each cluster into subclusters DEFINED ON THE REFERENCE GROUP
#'
#' WHY THE REFERENCE DEFINES THE SUBCLUSTERS, AND NOT ALL CELLS POOLED: the
#' question is "what inside this cluster differs from healthy". Clustering the
#' pooled cells would let the patient cells help draw the very boundaries they
#' are then tested against -- a disease-specific subset would carve out its own
#' subcluster and come back looking like a normal compartment that simply
#' contains patient cells. Fitting k-means on the HEALTHY cells alone makes each
#' subcluster a piece of normal biology, and assigning every study's cells to the
#' nearest healthy centroid asks the honest question: given the compartments a
#' healthy immune system has, where do patient cells fall, and do they look the
#' same once they get there?
#'
#' That yields the two ways a study can deviate, and the figure and table
#' downstream report both:
#'   OCCUPANCY  a subcluster holds a different share of the cluster's cells
#'   INTENSITY  the cells that are there express a marker differently
#'
#' A cluster whose reference has too few cells to subcluster safely is left whole
#' (subcluster 1), which is the honest answer rather than splitting noise.
#'
#' @param reference the group whose cells define the subclusters (e.g. controls)
#' @param k maximum subclusters per cluster; reduced when the reference is small
#' @param cells Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`, `sample_id`, `population_label` and one column per marker.
#' @param markers Character vector of marker names to use.
#' @param group_col Name of the column holding the biological grouping (cohort). Default `"cohort"`.
#' @param min_ref Minimum reference-group cells required to fit. Default `150L`.
#' @param min_cell Minimum cells a population must have before it is split. Default `20L`.
#' @param seed Random seed. The RNG stream is restored afterwards, so this cannot perturb later sampling. Default `42L`.
#' @return integer vector, one subcluster id per row of `cells`
#' @export
subcluster_by_reference <- function(cells, markers, group_col = "cohort",
                                    reference = NULL, k = 3L, min_ref = 150L,
                                    min_cell = 20L, seed = 42L) {
  n <- nrow(cells); sub <- rep(1L, n)
  if (!group_col %in% names(cells) || !"population_label" %in% names(cells)) return(sub)
  markers <- intersect(markers, names(cells))
  markers <- markers[vapply(cells[markers], is.numeric, logical(1))]
  if (length(markers) < 2L) return(sub)
  grp <- as.character(cells[[group_col]])
  if (is.null(reference) || !any(grp == reference, na.rm = TRUE)) return(sub)

  # RNG hygiene, for the reason documented above fig_gating_qc: kmeans() draws
  # from the global stream, and STEP 6 picks UMAP cells with sample() for every
  # LATER panel. Consuming draws here would silently change which cells get
  # embedded in a multi-panel run. Save the stream and put it back.
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  M <- as.matrix(cells[, markers, drop = FALSE])
  for (pp in unique(as.character(cells$population_label))) {
    idx <- which(!is.na(cells$population_label) & cells$population_label == pp)
    if (length(idx) < min_cell) next
    r <- idx[!is.na(grp[idx]) & grp[idx] == reference]
    Xr <- M[r, , drop = FALSE]
    Xr <- Xr[stats::complete.cases(Xr), , drop = FALSE]
    if (nrow(Xr) < min_ref) next
    # `k` may be a single number (the default, every population split the same
    # way) OR a named vector population -> k, as produced by
    # choose_subcluster_k() under --auto-subcluster-k. A population absent from a
    # named vector falls back to the scalar default, so a partially-populated
    # vector cannot silently leave a population unsplit.
    k_pp <- if (!is.null(names(k)) && pp %in% names(k)) k[[pp]] else k[[1]]
    kk <- max(2L, min(as.integer(k_pp), floor(nrow(Xr) / 75)))
    km <- try(stats::kmeans(Xr, centers = kk, nstart = 5L, iter.max = 50L), silent = TRUE)
    if (inherits(km, "try-error")) next
    # Order the centroids so subcluster 1 is always the same compartment between
    # runs; k-means labels are arbitrary and would otherwise permute per run,
    # making "cluster 4 subcluster 2" mean something different each time.
    cen <- km$centers[order(km$centers[, 1]), , drop = FALSE]
    X <- M[idx, , drop = FALSE]
    dm <- vapply(seq_len(nrow(cen)), function(j)
      rowSums((X - matrix(cen[j, ], nrow(X), ncol(X), byrow = TRUE))^2), numeric(nrow(X)))
    a <- max.col(-dm, ties.method = "first")
    a[!stats::complete.cases(X)] <- 1L
    # Letter the subclusters by their POSITION ON THE EMBEDDING (left to right)
    # rather than by a marker coordinate. k-means centroid order is arbitrary
    # with respect to the UMAP, so without this "4a, 4b, 4c" bear no relation to
    # where those compartments actually sit -- and the reader is comparing peak
    # rows against a map. Ordering on median UMAP-1 makes the letters trace the
    # cluster left to right, so adjacent letters are adjacent in space and a
    # densely packed pair like 4a/4c ends up lettered in the order it is seen.
    if (all(c("umap_1", "umap_2") %in% names(cells))) {
      lv_a <- sort(unique(a))
      pos  <- vapply(lv_a, function(s)
        stats::median(cells$umap_1[idx][a == s], na.rm = TRUE), numeric(1))
      pos[!is.finite(pos)] <- Inf
      a <- unname(setNames(rank(pos, ties.method = "first"), lv_a)[as.character(a)])
    }
    sub[idx] <- a
  }
  sub
}

#' Per (cluster, subcluster, marker, study): how far is this study from reference
#'
#' The figure shows the shift; this makes it sortable. Reports BOTH deviation
#' modes named above -- the occupancy of the subcluster and the intensity of the
#' marker within it -- so "which markers and which cells differ" is a sort on a
#' column rather than an eyeball over 288 panels.
#'
#' Effect size is the median difference in asinh units plus Cliff's delta, which
#' is bounded -1 to 1 and rank-based: at these per-subcluster n a difference in
#' means is dominated by outliers, and these distributions are not normal.
#' @param cells Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`, `sample_id`, `population_label` and one column per marker.
#' @param markers Character vector of marker names to use.
#' @param subcluster Integer vector of subcluster assignments, one per cell.
#' @param group_col Name of the column holding the biological grouping (cohort). Default `"cohort"`.
#' @param reference The group every other group is compared against.
#' @param min_cells Minimum cells a sample must contribute before it is used. Default `20L`.
#' @export
stats_subcluster_shifts <- function(cells, markers, subcluster, group_col = "cohort",
                                    reference = NULL, min_cells = 20L) {
  grp <- as.character(cells[[group_col]])
  ok  <- !is.na(grp) & nzchar(trimws(grp)) & !is.na(cells$population_label)
  if (is.null(reference) || !any(grp == reference, na.rm = TRUE)) return(NULL)
  others <- setdiff(sort(unique(grp[ok])), reference)
  if (!length(others)) return(NULL)
  markers <- intersect(markers, names(cells))
  sub <- as.integer(subcluster)
  cliff <- function(a, b) {
    if (!length(a) || !length(b)) return(NA_real_)
    if (length(a) > 2000L) a <- a[sample.int(length(a), 2000L)]
    if (length(b) > 2000L) b <- b[sample.int(length(b), 2000L)]
    r <- rank(c(a, b)); (2 * (sum(r[seq_along(a)]) - length(a) * (length(a) + 1) / 2) /
                          (length(a) * length(b))) - 1
  }
  out <- list()
  for (pp in sort(unique(as.character(cells$population_label[ok])))) {
    inp <- ok & cells$population_label == pp
    for (sb in sort(unique(sub[inp]))) {
      sel <- inp & sub == sb
      for (g in c(reference, others)) {
        if (g == reference) next
        a <- sel & grp == reference; b <- sel & grp == g
        if (sum(a) < min_cells || sum(b) < min_cells) next
        # occupancy: share of THIS study's cells in the cluster that land here
        occ_r <- sum(a) / max(1L, sum(inp & grp == reference)) * 100
        occ_g <- sum(b) / max(1L, sum(inp & grp == g)) * 100
        for (mk in markers) {
          xa <- as.numeric(cells[[mk]][a]); xa <- xa[is.finite(xa)]
          xb <- as.numeric(cells[[mk]][b]); xb <- xb[is.finite(xb)]
          if (length(xa) < min_cells || length(xb) < min_cells) next
          out[[length(out) + 1L]] <- data.frame(
            population = pp, subcluster = letters[pmin(as.integer(sb), 26L)],
            marker = mk, study = g,
            n_reference = length(xa), n_study = length(xb),
            occupancy_pct_reference = round(occ_r, 3), occupancy_pct_study = round(occ_g, 3),
            occupancy_pct_delta = round(occ_g - occ_r, 3),
            median_reference = round(stats::median(xa), 4),
            median_study = round(stats::median(xb), 4),
            median_delta = round(stats::median(xb) - stats::median(xa), 4),
            cliffs_delta = round(cliff(xb, xa), 4), stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  res[order(-abs(res$cliffs_delta)), , drop = FALSE]
}

#' FlowJo-style Multigraph Overlay: every cluster x every marker, one peak per group
#'
#' WHAT THIS IS: a port of FlowJo's Layout Editor "Make Multigraph Overlay ->
#' Histograms" (docs.flowjo.com/flowjo/graphical-reports/graph-options-and-
#' annotation/le-mgo/), which draws a histogram of every parameter with each
#' overlaid subset kept in its own colour. Here the grid is CLUSTER (row) x
#' MARKER (column) and the overlaid subsets are the STUDY GROUPS, so every panel
#' asks one question: inside this population, does this marker sit at the same
#' intensity in every cohort?
#'
#' WHY IT IS WORTH A FIGURE OF ITS OWN: umap_overview_by_group.png shows WHERE a
#' population sits and how big it is; this shows WHAT it expresses. A cohort can
#' carry a normal-sized cluster in a normal position whose marker intensity has
#' shifted, and no abundance figure in this pipeline can surface that -- they all
#' count cells rather than reading them.
#'
#' WHY THE PEAKS ARE MODE-NORMALISED (each curve scaled so its own maximum is
#' 100, which is FlowJo's "Count (%)" convention): the groups differ several-fold
#' in size, and within a rare population they differ far more. Plotted on a
#' shared density axis the largest cohort's curve simply towers over the others
#' and the figure silently becomes a cell-count comparison -- which
#' population_frequencies.png and group_comparison.png already do properly, with
#' statistics. Normalising to the mode puts every group on the same footing so
#' the eye compares SHAPE and POSITION, the only thing this figure is for.
#' Abundance is deliberately somebody else's job.
#'
#' MEMORY: the curves are precomputed with stats::density() rather than handing
#' ggplot the raw per-cell values, for exactly the reason documented above
#' fig_gating_qc -- a (population x marker x group) grid of raw cells rebuilds the
#' multi-million-row frame that OOM-killed this pipeline. Do not "simplify" this
#' into geom_density() either.
#'
#' @param cells embedding cell table (population_label, group_col, marker cols)
#' @param markers marker columns to draw, normally the UMAP feature set
#' @param min_cells a group needs at least this many cells in a population before
#'   its curve is drawn; below it a KDE is noise shaped like a result.
#' @param outfile Path to write the figure to.
#' @param group_col Name of the column holding the biological grouping (cohort). Default `"cohort"`.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `200`.
#' @param subcluster Integer vector of subcluster assignments, one per cell.
#' @param reference The group every other group is compared against.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_multigraph_overlay <- function(cells, outfile, markers, group_col = "cohort",
                                   min_cells = 20L, panel_label = "", dpi = 200,
                                   subcluster = NULL, reference = NULL, colors = fcs_colors()) {
  if (!group_col %in% names(cells) || !"population_label" %in% names(cells))
    return(invisible(NULL))
  markers <- intersect(markers, names(cells))
  markers <- markers[vapply(cells[markers], is.numeric, logical(1))]
  if (!length(markers)) {
    log_msg("[fig] multigraph overlay: no numeric markers, skipped"); return(invisible(NULL)) }

  grp  <- as.character(cells[[group_col]])
  ok   <- !is.na(grp) & nzchar(trimws(grp))
  groups <- sort(unique(grp[ok]))
  if (length(groups) < 2L) {
    log_msg("[fig] multigraph overlay: needs 2+ groups, skipped"); return(invisible(NULL)) }

  sub_of  <- if (is.null(subcluster)) rep(1L, nrow(cells)) else as.integer(subcluster)
  pop_all <- sort(unique(as.character(cells$population_label)))
  curves <- list(); n_skip <- 0L
  for (pp in pop_all) {
    in_pop0 <- ok & !is.na(cells$population_label) & cells$population_label == pp
    if (!any(in_pop0)) next
    for (sb in sort(unique(sub_of[in_pop0]))) {
    in_pop <- in_pop0 & sub_of == sb
    if (sum(in_pop) < min_cells) next
    for (mk in markers) {
      # ONE shared grid per panel, spanning every group's cells: curves drawn on
      # per-group ranges would each be stretched to the panel width and a real
      # shift between cohorts would vanish into the axis.
      xs <- as.numeric(cells[[mk]][in_pop]); xs <- xs[is.finite(xs)]
      if (length(xs) < 2L) next
      rng <- range(xs)
      if (!all(is.finite(rng)) || diff(rng) <= 0) next
      for (g in groups) {
        x <- as.numeric(cells[[mk]][in_pop & grp == g]); x <- x[is.finite(x)]
        if (length(x) < min_cells) { n_skip <- n_skip + 1L; next }
        dk <- stats::density(x, n = 256, from = rng[1], to = rng[2])
        pk <- max(dk$y)
        if (!is.finite(pk) || pk <= 0) next
        curves[[length(curves) + 1L]] <- data.frame(
          population = pp, sub = sb, marker = mk, group = g,
          value = dk$x, pct = dk$y / pk * 100, stringsAsFactors = FALSE)
      }
    }
    }
  }
  if (!length(curves)) {
    log_msg("[fig] multigraph overlay: no population had ", min_cells,
            "+ cells in 2 groups, skipped"); return(invisible(NULL)) }
  d <- do.call(rbind, curves)

  # Cluster NUMBERS over the sorted population names. Deliberately the SAME
  # scheme export_flowjo_fcs.R's code_of() uses for its Population channel, so
  # cluster 4 in this figure is code 4 in population_codes.csv and code 4 inside
  # the exported FCS files -- one numbering across figures, tables and the GUI.
  lv   <- sort(unique(as.character(cells$population_label)))
  num  <- setNames(seq_along(lv), lv)
  cells$.num <- unname(num[as.character(cells$population_label)])
  d$num <- unname(num[d$population])
  keylab  <- setNames(sprintf("%d \u00b7 %s", num, names(num)), names(num))
  pal_pop <- population_colours(lv, colors = colors)

  # --- UMAP row: pooled, then one panel per STUDY GROUP ----------------------
  # All panels are the SAME shared embedding restricted to different cells, so a
  # cluster keeps its position everywhere and the numbers can be read across.
  cells$.sub   <- sub_of
  cells$.panel <- "All samples"
  gsub_cells <- cells[ok, , drop = FALSE]; gsub_cells$.panel <- grp[ok]
  dd <- rbind(cells, gsub_cells)
  dd$.panel <- factor(dd$.panel, levels = c("All samples", groups))

  # One label per (cluster, SUBCLUSTER) -- "4a", "4b" -- so the compartments the
  # peaks below are keyed to are locatable on the map, not just their parent.
  centroids <- do.call(rbind, lapply(
    split(dd, list(dd$.panel, dd$population_label, dd$.sub), drop = TRUE), function(s) {
      if (nrow(s) < min_cells) return(NULL)
      data.frame(.panel = s$.panel[1], num = s$.num[1], sub = s$.sub[1],
                 umap_1 = stats::median(s$umap_1, na.rm = TRUE),
                 umap_2 = stats::median(s$umap_2, na.rm = TRUE), stringsAsFactors = FALSE) }))

  if (!is.null(centroids) && nrow(centroids)) {
    centroids$lab <- sprintf("%d%s", centroids$num, letters[pmin(centroids$sub, 26L)])
    # Subclusters of one cluster share a neighbourhood by construction, so their
    # median positions land almost on top of each other and the labels overprint
    # into an unreadable smear. Rather than hiding or jittering them at random,
    # collapse each colliding set onto its own centre and lay its members out
    # ADJACENT on a row: the group still marks one spot, and every member stays
    # legible and in a predictable order.
    # Iterative pairwise repulsion, not a single grouping pass: separating one
    # colliding pair routinely pushes a member into a THIRD label, and a
    # one-shot pass never revisits it -- which is what left "4a7a" and
    # "12c12b13c" overprinted. Each pass nudges every still-too-close pair
    # apart along the line between them and repeats until nothing moves, so
    # labels settle ADJACENT to their true position rather than being thrown
    # away from it.
    #
    # The threshold is elliptical because the text is far wider than it is
    # tall: two labels can sit 0.6 units apart vertically and read cleanly,
    # but need ~1.6 horizontally. A circular radius large enough to fix the
    # horizontal case would scatter vertically-stacked labels for no reason.
    #
    # Deliberately RNG-free (fixed trig offset for exact ties) -- this runs
    # inside a figure and must not touch the stream STEP 6 samples from.
    # The 2px gap the caller asks for is in PIXELS, but the separation below
    # works in UMAP data units, so convert once: the UMAP row is drawn
    # (W - legend) inches wide across (n groups + 1) panels at `dpi`, spanning
    # the x-range of the embedding. Anything that changes the figure width or
    # dpi therefore keeps the same visual gap instead of drifting.
    .W_est   <- max(15, 1.95 * length(markers) + 2.4)
    .pan_in  <- max(1, (.W_est - 2.6) / (length(groups) + 1L))
    .px_unit <- (.pan_in * dpi) / max(1e-6, diff(range(cells$umap_1, na.rm = TRUE)))
    .pad     <- 2 / .px_unit                      # 2 px, in UMAP units
    spread_labels <- function(cc, min_dx = 1.65 + .pad, min_dy = 0.62 + .pad,
                              iters = 400L) {
      n <- nrow(cc); if (n < 2L) return(cc)
      cc <- cc[order(cc$umap_1, cc$umap_2), , drop = FALSE]
      x <- cc$umap_1; y <- cc$umap_2
      for (it in seq_len(iters)) {
        moved <- FALSE
        for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
          dx <- (x[j] - x[i]) / min_dx; dy <- (y[j] - y[i]) / min_dy
          d  <- sqrt(dx * dx + dy * dy)
          if (d >= 1) next
          if (d < 1e-8) { dx <- cos(i + j); dy <- sin(i + j); d <- 1e-8 }
          ux <- dx / max(d, 1e-8); uy <- dy / max(d, 1e-8)
          push <- (1 - d) / 2 * 1.05
          x[i] <- x[i] - ux * push * min_dx; y[i] <- y[i] - uy * push * min_dy
          x[j] <- x[j] + ux * push * min_dx; y[j] <- y[j] + uy * push * min_dy
          moved <- TRUE
        }
        if (!moved) break
      }
      cc$umap_1 <- x; cc$umap_2 <- y; cc
    }
    centroids <- do.call(rbind, lapply(split(centroids, centroids$.panel, drop = TRUE),
                                       spread_labels))
  }

  ap <- auto_point_aes(nrow(cells))
  p_umap <- ggplot(dd, aes(umap_1, umap_2, colour = population_label)) +
    geom_point(size = ap$size, alpha = ap$alpha, stroke = 0) +
    scale_colour_manual(values = pal_pop, labels = keylab, name = "Cluster",
      guide = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
    facet_wrap(~ .panel, nrow = 1) +
    labs(title = "Shared CD45+ embedding, pooled, then one panel per study",
         subtitle = paste("labels are cluster+subcluster (4a, 4b, ...) and key the peak panels",
                          "below; they sit at the same position in every panel, nudged apart",
                          "only where they would overprint"),
         x = "UMAP 1", y = "UMAP 2") +
    theme_cyto(9, colors = colors) +
    theme(legend.position = "right", strip.text = element_text(size = 8))
  # Plain text, no box. The boxes were opaque enough to hide the very points they
  # sat on, and with one label per subcluster there are three times as many of
  # them. A short halo -- the same string drawn in the background colour at four
  # small offsets underneath -- keeps the text readable over dense points without
  # occluding anything.
  if (!is.null(centroids) && nrow(centroids)) {
    halo <- 0.055
    for (o in list(c(1,1), c(1,-1), c(-1,1), c(-1,-1)))
      p_umap <- p_umap + geom_text(
        data = centroids, aes(umap_1 + o[1]*halo, umap_2 + o[2]*halo, label = lab),
        inherit.aes = FALSE, size = 2.6, fontface = "bold",
        colour = colors$reference_fill, show.legend = FALSE)
    p_umap <- p_umap + geom_text(
      data = centroids, aes(umap_1, umap_2, label = lab), inherit.aes = FALSE,
      size = 2.6, fontface = "bold", colour = colors$label_text, show.legend = FALSE)
  }

  # --- peak row: one panel per (cluster, marker), titled "<cluster> - <marker>"
  # Panel titles are "cluster-subcluster-marker" (e.g. 4-2-CD3). Levels are laid
  # out cluster-major then subcluster then marker, so with ncol = one-per-marker
  # the grid still reads as a table: a row is one (cluster, subcluster), a column
  # is one marker, and consecutive rows are the subclusters of the same cluster.
  # Subclusters are lettered so a panel title cannot be misread: "4-2-CD3" has
  # two numbers meaning different things, "4-b-CD3" reads as cluster 4,
  # subcluster b, marker CD3 at a glance.
  sub_ltr <- function(i) letters[pmin(as.integer(i), 26L)]
  mk_lv <- intersect(markers, unique(d$marker))
  key   <- unique(d[, c("num", "sub")])
  key   <- key[order(key$num, key$sub), , drop = FALSE]
  d$lab <- factor(sprintf("%d-%s-%s", d$num, sub_ltr(d$sub), d$marker),
                  levels = unlist(lapply(seq_len(nrow(key)), function(i)
                    sprintf("%d-%s-%s", key$num[i], sub_ltr(key$sub[i]), mk_lv))))
  d$lab <- droplevels(d$lab)
  # Fixed, meaningful study colours rather than palette positions: the REFERENCE
  # group is always green (the baseline every other curve is read against) and
  # the patient studies take red then purple, in sorted order, so a colour means
  # the same study in every one of the ~300 panels and across reruns. Purple is
  # used here deliberately even though the population palette excludes violet --
  # three thick overlaid lines are a different discrimination problem from a
  # dozen point clouds, and green/red/purple are maximally separated for it.
  sp  <- colors$study_palette %||% c("#00A651", "#E8112D", "#8E44E8")
  ord <- c(intersect(reference, groups), setdiff(groups, reference))
  pal <- setNames(rep_len(sp, length(ord)), ord)

  # npp counts cluster-SUBCLUSTER rows, not clusters: with subclustering on, the
  # grid has one row per compartment (37 here, not 12), and sizing off the
  # cluster count alone collapses the peak panels to nothing.
  nmk <- length(mk_lv); npp <- nrow(key)
  p_pk <- ggplot(d, aes(value, pct, colour = group)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ lab, ncol = nmk, scales = "free_x") +
    scale_colour_manual(values = pal, name = "Study") +
    labs(title = "Multigraph overlay, one peak per study, for every cluster x marker",
         subtitle = paste("each curve is one study's distribution of that marker inside that",
                          "cluster, scaled to its own mode (FlowJo \"Count (%)\");",
                          "read peak POSITION and SHAPE, not height, abundance is in",
                          "group_comparison.png"),
         x = "asinh-transformed intensity", y = "% of mode") +
    theme_cyto(9, colors = colors) +
    theme(legend.position = "right",
          strip.text = element_text(size = 6.5),
          axis.text = element_text(size = 5.5))

  # ---- compose with grid, NOT patchwork ------------------------------------
  # DO NOT put these two back under patchwork's `/`. Once subclustering is on
  # the peaks plot carries one facet per (cluster, subcluster, marker) -- 288
  # here -- and patchwork composes that to a BLANK region: the UMAP row draws,
  # the peaks pane comes out empty at every dpi from 100 to 200, with no error
  # and no warning. The same peaks plot saved on its own renders correctly
  # (verified), so it is the composition step, not the plot. Laying the two
  # gtables into an explicit grid layout sidesteps patchwork's panel-alignment
  # logic entirely and is deterministic at any panel count.
  h_umap <- 4.6; h_pk <- max(4, 1.30 * npp)
  W <- max(15, 1.95 * nmk + 2.4); H <- h_umap + h_pk + 1.15
  # Same raster ceiling safe_ggsave enforces, applied here since we drive the
  # device directly.
  dpi_use <- max(36L, min(as.integer(dpi), floor(getOption("cyRAVEN.max_raster_px", 30000L) / max(W, H))))
  ttl <- paste0("Clusters and their marker profiles",
                if (nzchar(panel_label)) paste0(", ", panel_label) else "")
  cap <- paste0("panels are cluster-subcluster-marker; subclusters are lettered and were fitted on the ",
                "reference group only, then every study's cells assigned to the nearest reference centroid. ",
                "Cluster numbers match population_codes.csv and the Population channel of the FlowJo export. ",
                "Curves omitted where a study had fewer than ", min_cells, " cells in the compartment.")

  grDevices::png(outfile, width = W, height = H, units = "in", res = dpi_use, type = "cairo")
  ok_dev <- TRUE
  tryCatch({
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      4, 1, heights = grid::unit(c(0.42, h_umap, h_pk, 0.46),
                                 c("in", "null", "null", "in")))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(ttl, x = grid::unit(0.006, "npc"), hjust = 0,
                    gp = grid::gpar(fontface = "bold", cex = 1.35))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    grid::grid.draw(ggplotGrob(p_umap)); grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 3))
    grid::grid.draw(ggplotGrob(p_pk));   grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 4))
    grid::grid.text(paste(strwrap(cap, width = 200), collapse = "\n"),
                    x = grid::unit(0.006, "npc"), hjust = 0,
                    gp = grid::gpar(cex = 0.62, col = colors$caption_text))
    grid::popViewport(2)
  }, error = function(e) { ok_dev <<- FALSE
    log_msg("[fig] multigraph overlay: compose failed: ", conditionMessage(e)) })
  grDevices::dev.off()
  if (!ok_dev) return(invisible(NULL))

  log_msg("[fig] wrote ", outfile, " (", npp, " cluster-subcluster compartment(s) x ",
          nmk, " marker(s) x ", length(groups), " study/studies; UMAP panels: pooled + ",
          length(groups), "; ", dpi_use, " dpi",
          if (n_skip) paste0("; ", n_skip, " curve(s) below the ",
                             min_cells, "-cell floor") else "", ")")
  invisible(list(umap = p_umap, peaks = p_pk))
}

# =============================================================================
