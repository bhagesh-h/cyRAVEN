# SECTION 9b -- BETWEEN-GROUP ABUNDANCE COMPARISON
# =============================================================================

#' Pick the best available abundance measure, and say what it means
#'
#' WHY THIS IS A FUNCTION AND NOT A HARDCODED COLUMN: three quantities in the
#' frequency table look interchangeable and are not.
#'
#'   `count`           EVENT count. Depends on how long the operator ran the tube
#'                     and on --max-events-per-file. NEVER comparable between
#'                     samples; never a valid y-axis for a group comparison.
#'   `pct_of_cd45_pos` Comparable, but COMPOSITIONAL: populations are constrained
#'                     to sum to 100%, so one lineage expanding mathematically
#'                     forces the others down. A significant fall in "% NK" is
#'                     equally consistent with NK loss and with granulocyte gain.
#'   `cells_per_ul`    Absolute. Independent per population, so a change in one
#'                     says nothing about the others. Requires wbc_per_ul.
#'
#' Prefer absolute when available and fall back to frequency otherwise, carrying
#' the axis label and the compositional caveat with the choice so the figure
#' cannot silently misrepresent which of the two it is showing.
#' @param freq Population frequency table, one row per sample x population.
#' @param prefer Preferred choice; auto picks the best available. Default `c("auto", "absolute", "frequency")`.
#' @export
abundance_measure <- function(freq, prefer = c("auto", "absolute", "frequency")) {
  prefer <- match.arg(prefer)
  freq_m <- list(col = "pct_of_cd45_pos", label = "% of CD45+", absolute = FALSE,
                 caveat = paste("frequencies are compositional: a fall in one",
                                "population may reflect expansion of another, not",
                                "its own loss."))
  abs_ok <- "cells_per_ul" %in% names(freq) && any(is.finite(freq$cells_per_ul))
  if (prefer == "frequency") return(freq_m)
  if (abs_ok)
    return(list(col = "cells_per_ul", label = "cells / \u00b5L", absolute = TRUE,
                caveat = NULL))
  if (prefer == "absolute")
    freq_m$caveat <- paste(freq_m$caveat, "Absolute concentrations were requested",
                           "but wbc_per_ul is absent from the patient table, so",
                           "frequencies are shown instead.")
  freq_m
}

#' Test each population's abundance between groups
#'
#' WHAT: for every population, compares each non-reference group against the
#' reference group, and (with 3+ groups) also reports an omnibus test.
#'
#' WHY THESE TESTS: group sizes in this design are single-digit to low-double-digit
#' and cell-abundance distributions are right-skewed with outliers, so normality is
#' not a safe assumption and a t-test's nominal p-value would be unreliable.
#' Wilcoxon rank-sum (Mann-Whitney) assumes only independence and ordinal
#' comparability. Kruskal-Wallis is its multi-group extension. Both are exact-rank
#' based, so they behave sensibly at n = 6.
#'
#' WHY BOTH RAW AND ADJUSTED p: testing 12 populations x 2 groups is 24 tests, at
#' which ~1 in 5 chance of a spurious p < 0.05 somewhere is expected under the
#' null. Reporting only raw p invites reading noise as signal; reporting only
#' adjusted p hides a real trend in an underpowered study. Both are written, with
#' Benjamini-Hochberg FDR (not Bonferroni) because these tests are correlated and
#' the goal is a ranked shortlist, not a single confirmatory decision.
#'
#' @param freq population_frequencies table
#' @param group_of named character vector sample_id -> group label
#' @param reference the group all others are compared against (e.g. controls)
#' @param min_n minimum samples per group for a test to be attempted
#' @param value_col override the abundance measure: test this column instead
#'   of pct_of_cd45_pos/cells_per_ul. Used for tables shaped like `freq`
#'   (sample_id, population, ...) but carrying a different quantity, e.g.
#'   functional-marker pct_positive or a derived population ratio.
#' @param measure Which abundance measure to use. Default `"auto"`.
#' @return data.frame, one row per population x comparison
#' @export
stats_group_comparison <- function(freq, group_of, reference = NULL, min_n = 3L,
                                   measure = "auto", value_col = NULL) {
  meas <- if (!is.null(value_col)) list(col = value_col) else abundance_measure(freq, measure)
  d <- freq[!is.na(freq[[meas$col]]), , drop = FALSE]
  d <- qc_pass_rows(d)
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  groups <- sort(unique(d$group))
  if (length(groups) < 2L) return(NULL)
  if (is.null(reference) || !reference %in% groups) reference <- groups[1]
  others <- setdiff(groups, reference)

  rows <- list()
  for (pop in sort(unique(d$population))) {
    dp <- d[d$population == pop, , drop = FALSE]
    vals <- split(dp[[meas$col]], dp$group)
    n_by <- vapply(groups, function(g) length(vals[[g]] %||% numeric(0)), integer(1))
    # Omnibus first, so a reader can see whether ANY group differs before reading
    # the pairwise columns.
    p_omni <- NA_real_
    if (length(groups) > 2L && all(n_by >= min_n)) {
      kw <- try(stats::kruskal.test(dp[[meas$col]], factor(dp$group)), silent = TRUE)
      if (!inherits(kw, "try-error")) p_omni <- kw$p.value
    }
    for (g in others) {
      a <- vals[[reference]]; b <- vals[[g]]
      if (is.null(a) || is.null(b) || length(a) < min_n || length(b) < min_n) next
      wt <- try(stats::wilcox.test(a, b, exact = FALSE), silent = TRUE)
      p <- if (inherits(wt, "try-error")) NA_real_ else wt$p.value
      # Report the effect on the scale the reader sees, plus a scale-free version.
      # Cliff's delta (rank-biserial) is bounded \code{[-1, 1]} and is interpretable at
      # small n where a ratio of means is dominated by single outliers.
      cliff <- if (length(a) && length(b))
        mean(outer(b, a, ">")) - mean(outer(b, a, "<")) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        population = pop, measure = meas$col,
        reference_group = reference, comparison_group = g,
        n_reference = length(a), n_comparison = length(b),
        median_reference = median(a), median_comparison = median(b),
        mean_reference = mean(a), mean_comparison = mean(b),
        sd_reference = stats::sd(a), sd_comparison = stats::sd(b),
        fold_change = if (median(a) > 0) median(b) / median(a) else NA_real_,
        cliffs_delta = round(cliff, 3),
        test = "Wilcoxon rank-sum", p_value = p,
        p_omnibus_kruskal = p_omni,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  # FDR across ALL tests in the family, which is every population x group
  # comparison in this figure -- the family a reader scans when hunting for a hit.
  out$p_adj_BH <- stats::p.adjust(out$p_value, method = "BH")
  out$significant_raw <- !is.na(out$p_value) & out$p_value < 0.05
  out$significant_BH  <- !is.na(out$p_adj_BH) & out$p_adj_BH < 0.05
  out[order(out$p_value), ]
}

#' Format a p-value the way these figures conventionally annotate them
#' @param p The p.
#' @param style The style. Default `c("stars_and_trend", "stars", "numeric")`.
#' @keywords internal
p_annotation <- function(p, style = c("stars_and_trend", "stars", "numeric")) {
  style <- match.arg(style)
  if (!is.finite(p)) return(NA_character_)
  if (style == "numeric") return(sprintf("p: %.4f", p))
  stars <- if (p < 0.0001) "****" else if (p < 0.001) "***" else
           if (p < 0.01) "**" else if (p < 0.05) "*" else NA_character_
  if (!is.na(stars)) return(stars)
  # Trends are printed as an exact number rather than "ns": at these group sizes
  # p = 0.06 is a result worth a reader's eye, and hiding it behind "ns" loses
  # information that an exact value preserves without overclaiming.
  if (style == "stars_and_trend" && p < 0.10) return(sprintf("p: %.4f", p))
  NA_character_
}

#' Grouped abundance figure: one panel per population, bar + SD + every sample
#'
#' WHAT: mean bar, SD whiskers, one point per sample, significance brackets from
#' the reference group to each other group. Panels are lettered.
#'
#' WHY POINTS AND NOT JUST BARS: at n = 6-11 a bar and a whisker hide whether a
#' "difference" rests on one outlier or on a consistent shift. Every sample is
#' drawn, so the reader can see the shape of the evidence -- which is also why the
#' bar is left unfilled for the reference group: the fill carries group identity,
#' not emphasis.
#'
#' @param freq population_frequencies table
#' @param group_of named vector sample_id -> group
#' @param stats output of stats_group_comparison(); NULL to skip brackets
#' @param p_source "raw" or "BH" -- which p-value the brackets display
#' @param value_col,value_label,value_caveat override the abundance measure
#'   with an arbitrary numeric column of `freq` (must still be one row per
#'   sample_id x population). Used to draw the same bar/median/whisker/points
#'   layout, with the same Wilcoxon-bracket machinery, for tables that share
#'   `freq`'s shape but not its meaning -- functional-marker positivity, or a
#'   derived population ratio.
#' @param title_noun the figure title's subject, e.g. "Population abundance"
#'   or "Functional marker positivity"; " by group" is appended automatically
#'   when 2+ groups are plotted.
#' @param outfile Path to write the figure to.
#' @param reference The group every other group is compared against.
#' @param measure Which abundance measure to use. Default `"auto"`.
#' @param ncol Number of panel columns; NULL computes one that keeps the canvas roughly square. Default `3L`.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_group_comparison <- function(freq, outfile, group_of, stats = NULL,
                                 reference = NULL, p_source = c("raw", "BH"),
                                 measure = "auto", ncol = 3L,
                                 panel_label = "", dpi = 300,
                                 value_col = NULL, value_label = NULL,
                                 value_caveat = NULL,
                                 title_noun = "Population abundance",
                                 colors = fcs_colors()) {
  p_source <- match.arg(p_source)
  meas <- if (!is.null(value_col))
            list(col = value_col, label = value_label %||% value_col,
                 absolute = FALSE, caveat = value_caveat)
          else abundance_measure(freq, measure)
  n_ctrl <- if ("is_control" %in% names(freq))
    length(unique(freq$sample_id[freq$is_control])) else 0L
  d <- freq[!is.na(freq[[meas$col]]), , drop = FALSE]
  d <- qc_pass_rows(d)
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) { log_msg("[fig] no grouped samples, group comparison skipped"); return(invisible(NULL)) }
  groups <- sort(unique(d$group))
  # A single group is a legitimate layout, not a failure: it is exactly the
  # "all samples pooled" view used for the frequency overview, where the value of
  # the panel is the mean/SD/median/spread of one group rather than a comparison.
  if (length(groups) < 2L) stats <- NULL
  if (is.null(reference) || !reference %in% groups) reference <- groups[1]
  # Reference group first on the x-axis, so every bracket reads left-to-right.
  glev <- c(reference, setdiff(groups, reference))
  d$group <- factor(d$group, levels = glev)

  # Reference group unfilled (white), others take hues. Matches the convention in
  # this literature and keeps the reference visually neutral.
  fills <- if (length(glev) == 1L) setNames(colors$reference_fill, glev)
           # Pass the reference this figure actually resolved, not the global
           # default. population_colours() only reaches the separated
           # study palette when the reference is one of the levels, and the
           # global default is whatever some earlier call left set -- so the
           # same cohort could come out green/red/purple in one run and
           # red/orange/yellow in another, purely on call order.
           else population_colours(glev, reference = reference, colors = colors)
  # One group needs no key -- it would label a single bar with what the title
  # already says -- and a full-width bar reads as a filled panel rather than a bar.
  show_key <- length(glev) > 1L
  bar_w <- if (length(glev) == 1L) 0.30 else 0.62

  # Order panels by overall abundance so the eye moves from common to rare
  # populations rather than alphabetically.
  ord <- stats::aggregate(stats::reformulate("population", meas$col), d, median)
  pops <- ord$population[order(-ord[[meas$col]])]

  # ncol = 3 reads well at the population counts this was developed against.
  # fig_functional_markers() can hand this population x marker panels -- a
  # full marker panel across several functional blocks easily passes 100 --
  # and a fixed 3-wide column would turn that into a page many yards tall.
  # Only auto-widen when the caller left ncol at its default: an explicit
  # ncol is a deliberate choice and stays honoured.
  if (missing(ncol) && length(pops) > 30L)
    ncol <- ceiling(sqrt(length(pops) * 1.6))

  panels <- list()
  for (i in seq_along(pops)) {
    pop <- pops[i]
    # The panel key is the y-axis title, and a y-axis title is ROTATED, so the
    # space it has is the panel's HEIGHT, not its width. A population name like
    # "CD4 T cells" fits; the key fig_functional_markers() supplies does not --
    # "Homing receptors on T cells: CD8 T cells, CXCR3" is 47 characters and
    # overflowed into the neighbouring panel, landing on its axis numbers and
    # its tag letter. Wrapping trades length for lines, which the panel has
    # room for.
    ylab_wrapped <- wrap_axis_title(pop)
    ylab_lines <- length(strsplit(ylab_wrapped, "\n", fixed = TRUE)[[1]])
    dp <- d[d$population == pop, , drop = FALSE]
    smry <- do.call(rbind, lapply(glev, function(g) {
      v <- dp[[meas$col]][dp$group == g]
      if (!length(v)) return(NULL)
      data.frame(group = g, mean = mean(v), sd = stats::sd(v),
                 median = stats::median(v), n = length(v),
                 stringsAsFactors = FALSE)
    }))
    if (is.null(smry)) next
    smry$group <- factor(smry$group, levels = glev)
    smry$sd[is.na(smry$sd)] <- 0
    ymax <- max(c(dp[[meas$col]], smry$mean + smry$sd), na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1

    # Brackets: reference vs each other group, only where an annotation exists.
    ann <- NULL
    if (!is.null(stats)) {
      st <- stats[stats$population == pop, , drop = FALSE]
      if (nrow(st)) {
        pv <- if (p_source == "BH") st$p_adj_BH else st$p_value
        lab <- vapply(pv, p_annotation, character(1))
        keep <- which(!is.na(lab))
        if (length(keep)) {
          ann <- data.frame(
            x1 = 1L,
            x2 = match(st$comparison_group[keep], glev),
            label = lab[keep], stringsAsFactors = FALSE)
          # Stack brackets so multiple comparisons cannot collide (section 9.1).
          ann <- ann[order(ann$x2), , drop = FALSE]
          ann$y <- ymax * (1.10 + 0.13 * (seq_len(nrow(ann)) - 1L))
        }
      }
    }
    headroom <- if (!is.null(ann)) max(ann$y) * 1.10 else ymax * 1.06

    g <- ggplot(smry, aes(group, mean)) +
      geom_col(aes(fill = group), colour = colors$bar_outline, width = bar_w, linewidth = 0.35) +
      geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
                    width = 0.20, linewidth = 0.35, colour = colors$bar_outline) +
      # Median as a rule that OVERHANGS the bar on both sides. The bar height is
      # the mean, so mean and median are both readable off one panel and their
      # divergence -- the signature of a skewed population, or of one outlier
      # dragging the mean -- is visible directly instead of having to be inferred.
      # The overhang, not a colour, is what distinguishes it: a hue would risk
      # colliding with a group fill and would vanish under colour-vision
      # deficiency, whereas the geometry reads at any n of groups (section 4.5).
      geom_segment(aes(x = as.integer(group) - bar_w * 0.65,
                       xend = as.integer(group) + bar_w * 0.65,
                       y = median, yend = median),
                   linewidth = 0.7, colour = colors$bracket) +
      # Raw samples on top of the summary; jittered only horizontally so no point's
      # vertical position -- its actual value -- is distorted.
      geom_point(data = dp, aes(group, .data[[meas$col]]),
                 position = position_jitter(width = 0.13, height = 0),
                 size = 1.15, colour = colors$bracket, inherit.aes = FALSE) +
      scale_fill_manual(values = fills, name = NULL, drop = FALSE) +
      scale_y_continuous(labels = sci_labels, expand = expansion(mult = c(0, 0.02)),
                         limits = c(0, headroom)) +
      labs(x = NULL, y = ylab_wrapped) +
      theme_cyto(colors = colors) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position = "none",
            axis.title.y = element_text(size = 7.5, lineheight = 0.95,
                                        margin = margin(r = 3)),
            # A rotated title grows LEFTWARD into whatever sits beside it, so
            # each extra wrapped line needs its own strip of left margin or the
            # text lands on the neighbouring panel's axis numbers and tag.
            plot.margin = margin(2, 5, 2, 2 + 7 * (ylab_lines - 1L)))
    if (!is.null(ann)) {
      for (r in seq_len(nrow(ann))) {
        yy <- ann$y[r]; tick <- ymax * 0.030
        g <- g +
          annotate("segment", x = ann$x1[r], xend = ann$x2[r], y = yy, yend = yy,
                   linewidth = 0.3, colour = colors$bracket) +
          annotate("segment", x = ann$x1[r], xend = ann$x1[r],
                   y = yy - tick, yend = yy, linewidth = 0.3, colour = colors$bracket) +
          annotate("segment", x = ann$x2[r], xend = ann$x2[r],
                   y = yy - tick, yend = yy, linewidth = 0.3, colour = colors$bracket) +
          annotate("text", x = (ann$x1[r] + ann$x2[r]) / 2, y = yy + tick * 0.6,
                   label = ann$label[r], size = 2.35, vjust = 0, colour = colors$bracket)
      }
    }
    panels[[length(panels) + 1L]] <- g
  }
  if (!length(panels)) return(invisible(NULL))

  # One shared legend for the whole figure: group identity is figure-wide, so
  # repeating it in 12 panels would spend the label budget over and over (section 2.2).
  # patchwork collects it from the panels, so the guide always matches the fills
  # actually drawn rather than a separately-constructed copy that could drift.
  if (show_key)
    # TOP-LEFT, NOT RIGHT. patchwork collects the guide once for the whole
    # composition and centres it vertically on the right. On a tall grid --
    # functional_markers.png runs to 26 panels and six thousand pixels -- that
    # puts the key thousands of pixels below the top of the image, so a reader
    # has to scroll away from the figure to find out what the colours mean and
    # scroll back. Top-left puts it where reading starts.
    panels[[1]] <- panels[[1]] +
      theme(legend.position = "top", legend.justification = "left",
            legend.direction = "horizontal",
            legend.key.size = unit(9, "pt"),
            legend.margin = margin(b = 2),
            legend.text = element_text(size = 7))

  nr <- ceiling(length(panels) / ncol)
  fig <- patchwork::wrap_plots(panels, ncol = ncol, guides = "collect") +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = paste0(if (length(glev) > 1L) paste0(title_noun, " by group")
                     else title_noun,
                     if (nzchar(panel_label)) paste0(", ", panel_label) else ""),
      # The subtitle states the measure and what every mark means, because a bar
      # chart with points is ambiguous otherwise: bar height could be mean, median
      # or sum, and whiskers could be SD, SEM or CI (section 1.3, section 3.2).
      subtitle = paste0(
        meas$label, "; bar = mean, black rule = median, whiskers = SD, ",
        "points = individual samples",
        if (!is.null(stats)) paste0(
          "; brackets: Wilcoxon rank-sum vs ", reference,
          if (p_source == "BH") " (BH-adjusted)" else " (unadjusted)") else "",
        if (n_ctrl > 0) sprintf("; %d unstained control%s excluded", n_ctrl,
                                if (n_ctrl > 1) "s" else "") else "",
        if (!is.null(meas$caveat)) paste0("\n", meas$caveat) else ""),
      # Newline-separated, not one long line: a single-line caption is silently
      # CLIPPED at the device edge, losing the significance key entirely.
      caption = paste0(
        sprintf("n = %s samples.", paste(sprintf("%s: %d", glev,
                vapply(glev, function(g) length(unique(d$sample_id[d$group == g])),
                       integer(1))), collapse = ", ")),
        if (length(glev) > 1L)
          paste0("  Bars left to right: ", paste(glev, collapse = ", "), ".") else "",
        if (!is.null(stats))
          "\n* p<0.05   ** p<0.01   *** p<0.001   **** p<0.0001;   0.05 <= p < 0.10 printed as an exact value."
        else ""),
      theme = ggplot2::theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, lineheight = 1.15),
        plot.caption = element_text(size = 7, hjust = 0, colour = colors$caption_text))) &
    ggplot2::theme(plot.tag = element_text(size = 10, face = "bold"))
  safe_ggsave(outfile, plot = fig, width = 2.75 * ncol + 0.6, height = 2.15 * nr + 1.1,
             dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile, " (", length(panels), " populations, ",
          length(glev), " groups, y = ", meas$col, ")")
  invisible(fig)
}

#' Wrap a panel's axis title so a long key does not overrun its neighbour
#'
#' WHY THIS IS NEEDED. The panel key becomes the y-axis title, and a y-axis
#' title is rotated, so the room it has is the panel's HEIGHT rather than its
#' width. A population name fits comfortably. The key
#' [fig_functional_markers()] supplies does not: it is
#' "block: population, marker", which reaches 47 characters and more, and an
#' unwrapped title of that length runs out of the panel and lands on the axis
#' numbers and tag letter of the panel beside it.
#'
#' WHY IT REBALANCES RATHER THAN WRAPPING AT A FIXED WIDTH. `strwrap()` at a
#' fixed width leaves the last line short, so a 47-character title at width 30
#' gives one full line and one stub. Choosing the line count first and then
#' dividing gives lines of roughly equal length, which reads better rotated and
#' keeps the block from looking ragged.
#'
#' @param x Character scalar, the title to wrap.
#' @param width Target characters per line before wrapping starts.
#' @param max_lines Ceiling on line count. Beyond this the lines get longer
#'   rather than more numerous, because a title taller than the panel is wide is
#'   worse than a slightly long one.
#' @return `x` with newlines inserted, unchanged when it already fits.
#' @keywords internal
wrap_axis_title <- function(x, width = 30L, max_lines = 3L) {
  x <- as.character(x)
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(x)
  if (nchar(x) <= width) return(x)
  n <- min(max_lines, ceiling(nchar(x) / width))
  w <- max(width, ceiling(nchar(x) / n))
  paste(strwrap(x, width = w), collapse = "\n")
}

#' Axis labels in the 1x10^7 form these figures conventionally use
#'
#' WHY NOT scales::label_scientific(): it renders "1e+07", which is code notation.
#' WHY NOT plain numerals: absolute cell concentrations span 10^3 to 10^7 across
#' populations, so full numerals make the tick column wider than the panel.
#' Returns an expression vector so the exponent renders as a true superscript;
#' a plain "10^7" string would print the caret literally.
#' @param x A vector of values.
#' @keywords internal
sci_labels <- function(x) {
  txt <- vapply(x, function(v) {
    if (!is.finite(v)) return("")
    if (v == 0) return("0")
    if (abs(v) < 1000 && abs(v) >= 0.01)
      return(format(v, trim = TRUE, drop0trailing = TRUE))
    e <- floor(log10(abs(v))); m <- v / 10^e
    if (abs(m - 1) < 1e-8) sprintf("10^%d", e)
    else sprintf("%s%%*%%10^%d", format(round(m, 1), trim = TRUE, drop0trailing = TRUE), e)
  }, character(1))
  parse(text = ifelse(nzchar(txt), txt, "\"\""))
}

#' Population frequency figure -- one lettered panel per population
#'
#' Delegates to fig_group_comparison() so the frequency overview and the
#' between-group comparison are THE SAME LAYOUT. The previous version was one
#' dodged bar per sample per population on a shared axis, which showed each
#' sample's value but no summary at all: mean, median and spread could not be
#' read off it, and the rare populations were invisible next to the abundant ones
#' because every population shared one x-axis.
#'
#' Giving each population its own panel gives each its own y-scale, so a
#' population at 0.3% is as legible as one at 40%.
#'
#' @param freq population_frequencies table; rows flagged `is_control` are
#'   DROPPED, because a control's percentages are fractions of a quantile-
#'   fallback slice rather than of a real CD45+ population -- plotting them
#'   beside stained samples invites a false comparison.
#'   split by group with tests; when absent all samples form one group and the
#'   panels show the pooled mean/median/SD/points.
#' @param outfile Path to write the figure to.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_population_frequencies <- function(freq, outfile, panel_label = "", dpi = 300,
                                       colors = fcs_colors()) {
  # Deliberately POOLED and always in per-cent, even when a group column and a
  # blood count exist: this figure answers "what does the cohort look like", and
  # group_comparison.png answers "do the groups differ". Splitting this one by
  # group too would emit the same figure twice under two filenames.
  ids <- unique(freq$sample_id)
  group_of <- setNames(rep("All samples", length(ids)), ids)
  fig_group_comparison(freq, outfile, group_of = group_of, stats = NULL,
                       measure = "frequency", panel_label = panel_label, dpi = dpi,
                       colors = colors)
}

#' Collapse functional_markers' (block, population, marker) into one panel
#' key, shared by fig_functional_markers() and the stats call in main() so a
#' figure's bracket always matches the right row in functional_markers_stats.csv
#' @param fx The fx.
#' @keywords internal
fx_panel_key <- function(fx) {
  fx$population <- paste0(fx$block, ": ", fx$population, ", ", fx$marker)
  fx
}

#' Functional-marker positivity figure
#'
#' Delegates to fig_group_comparison(), exactly as fig_population_frequencies()
#' does, so a monocyte HLA-DR shift gets the same bar/median/whisker/points
#' treatment and the same Wilcoxon brackets as a population-abundance shift --
#' before this figure existed, functional_markers.csv was written but never
#' plotted, so a functional difference between cohorts was visible only as
#' rows in a CSV.
#'
#' @param fx functional_markers table (sample_id, block, population, marker,
#'   pct_positive, is_control). Population x marker becomes the panel key,
#'   because a positivity threshold only means something inside the
#'   population it was derived for.
#' @param group_of optional sample_id -> group map; NULL pools every sample
#'   into "All samples" (mirrors fig_population_frequencies()).
#' @param outfile Path to write the figure to.
#' @param stats Statistics table from the matching stats_ function, used to annotate the figure.
#' @param reference The group every other group is compared against.
#' @param p_source Which p-value the figure annotates: raw or BH-adjusted. Default `c("raw", "BH")`.
#' @param ncol Number of panel columns; NULL computes one that keeps the canvas roughly square. Default `3L`.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_functional_markers <- function(fx, outfile, group_of = NULL, stats = NULL,
                                   reference = NULL, p_source = c("raw", "BH"),
                                   ncol = 3L, panel_label = "", dpi = 300,
                                   colors = fcs_colors()) {
  p_source <- match.arg(p_source)
  if (is.null(fx) || !nrow(fx)) return(invisible(NULL))
  d <- fx_panel_key(fx)
  if (is.null(group_of)) {
    ids <- unique(d$sample_id)
    group_of <- setNames(rep("All samples", length(ids)), ids)
    stats <- NULL
  }
  fig_group_comparison(
    d, outfile, group_of = group_of, stats = stats, reference = reference,
    p_source = p_source, ncol = ncol, panel_label = panel_label, dpi = dpi,
    value_col = "pct_positive", value_label = "% positive",
    value_caveat = paste(
      "positivity = fraction of cells above the per-sample gating threshold",
      "for that marker (thresholds_used.csv); population x marker",
      "combinations with fewer than 10 cells are omitted."),
    title_noun = "Functional marker positivity", colors = colors)
}

#' Derived abundance ratios (e.g. CD4:CD8), from the config's `ratios:` block
#'
#' WHY CONFIG-DRIVEN, NOT HARDCODED: population names live in `populations:`
#' (section 5.4) and are study-specific -- a ratio hardcoded to "CD4 T cells"/"CD8 T
#' cells" would silently mean nothing (or error) under a config using
#' different labels. Each ratio supplies its own numerator/denominator
#' population name, the same way each functional_blocks entry supplies its
#' own marker list, so this stays generic across cohorts and panels.
#'
#' @param freq population_frequencies table (sample_id, population,
#'   pct_of_cd45_pos, is_control)
#' @param ratios named list from cfg$ratios; each entry has `numerator` and
#'   `denominator` (population names from `populations:`) and an optional
#'   `label`
#' @return data.frame(sample_id, ratio, population (=label), numerator,
#'   denominator, numerator_pct, denominator_pct, value, is_control), or NULL
#' @export
compute_population_ratios <- function(freq, ratios) {
  if (!length(ratios) || is.null(freq)) return(NULL)
  d <- qc_pass_rows(freq)
  # sample_id -> panel, so a sample's ratio row can be filtered/split by panel
  # the same way freq/fx already can (a sample belongs to exactly one panel,
  # so numerator and denominator are always scored within the same one).
  panel_of <- if ("panel" %in% names(d)) setNames(d$panel, d$sample_id) else NULL
  rows <- lapply(names(ratios), function(rn) {
    r <- ratios[[rn]]
    num <- d[d$population == r$numerator,   c("sample_id", "pct_of_cd45_pos")]
    den <- d[d$population == r$denominator, c("sample_id", "pct_of_cd45_pos")]
    if (!nrow(num) || !nrow(den)) {
      log_msg("NOTE ratio '", rn, "' skipped: population '",
              if (!nrow(num)) r$numerator else r$denominator,
              "' was not scored in any sample")
      return(NULL)
    }
    names(num)[2] <- "numerator_pct"; names(den)[2] <- "denominator_pct"
    m <- merge(num, den, by = "sample_id")
    if (!nrow(m)) return(NULL)
    # Denominator = 0 is a real gate outcome (no CD8 T cells found), not a
    # data error, but a ratio against it is undefined rather than infinite.
    m$value <- ifelse(m$denominator_pct > 0,
                      m$numerator_pct / m$denominator_pct, NA_real_)
    data.frame(sample_id = m$sample_id,
               panel = if (is.null(panel_of)) NA_character_ else panel_of[m$sample_id],
               ratio = rn, population = r$label %||% rn,
               numerator = r$numerator, denominator = r$denominator,
               numerator_pct = m$numerator_pct,
               denominator_pct = m$denominator_pct,
               value = m$value, is_control = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) return(NULL)
  out
}

#' Population-ratio figure -- same layout as fig_group_comparison(), one panel
#' per ratio defined in the config's `ratios:` block
#' @param rt The rt.
#' @param outfile Path to write the figure to.
#' @param group_of Named character vector mapping sample_id to group label.
#' @param stats Statistics table from the matching stats_ function, used to annotate the figure.
#' @param reference The group every other group is compared against.
#' @param p_source Which p-value the figure annotates: raw or BH-adjusted. Default `c("raw", "BH")`.
#' @param ncol Number of panel columns; NULL computes one that keeps the canvas roughly square. Default `3L`.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `300`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_population_ratios <- function(rt, outfile, group_of = NULL, stats = NULL,
                                  reference = NULL, p_source = c("raw", "BH"),
                                  ncol = 3L, panel_label = "", dpi = 300,
                                  colors = fcs_colors()) {
  p_source <- match.arg(p_source)
  if (is.null(rt) || !nrow(rt)) return(invisible(NULL))
  if (is.null(group_of)) {
    ids <- unique(rt$sample_id)
    group_of <- setNames(rep("All samples", length(ids)), ids)
    stats <- NULL
  }
  fig_group_comparison(
    rt, outfile, group_of = group_of, stats = stats, reference = reference,
    p_source = p_source, ncol = ncol, panel_label = panel_label, dpi = dpi,
    value_col = "value", value_label = "ratio",
    value_caveat = paste(
      "ratio = numerator's pct_of_cd45_pos / denominator's pct_of_cd45_pos,",
      "per sample (population_ratios.csv); samples missing either",
      "population are dropped."),
    title_noun = "Population ratio", colors = colors)
}

# =============================================================================
