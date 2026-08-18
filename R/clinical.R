# =============================================================================
# SECTION 14 -- CLINICAL VARIABLES AGAINST POPULATIONS AND MARKERS
# =============================================================================
#
# WHY THIS FILE EXISTS. A sample sheet may carry clinical variables that are
# neither the study group nor a nuisance covariate: a severity score, a
# laboratory value, an outcome. Before this, the sheet recognised a fixed list of
# subject attributes -- patient_id, sex, age_years, height_cm, weight_kg,
# infection_focus, cohort, wbc_per_ul -- and anything outside it was carried as a
# "study column", usable only to GROUP samples. A SOFA score cannot sensibly
# group anything, and a 28-day survival flag grouped as a cohort loses the fact
# that it is an outcome rather than a design variable. So those columns were read
# from the sheet, written into the manifest, and never analysed.
#
# WHAT THIS ADDS. One question, asked per population and per marker: does this
# quantity move with that clinical variable?
#
#   numeric variable  (SOFA, CRP, lactate, BMI)  -> Spearman's rho
#   two-level variable (survival yes/no)          -> Wilcoxon rank-sum + Cliff's delta
#   many-level variable (infection focus)         -> Kruskal-Wallis
#
# WHY RANK METHODS THROUGHOUT. Same reason the rest of the package uses them.
# Clinical scores are ordinal by construction, laboratory values are skewed and
# carry outliers that are real rather than erroneous, and cohorts are small. A
# Pearson correlation on nine points with one extreme creatinine is a statement
# about that one patient.
#
# WHY ONE VALUE PER SAMPLE. Every test here runs on the per-sample frequency or
# the per-sample median intensity, never on pooled cells. Correlating a clinical
# score against tens of thousands of events treats one deeply acquired patient as
# tens of thousands of independent observations, and the p-value that comes out
# is a statement about acquisition depth.
#
# WHAT IS DELIBERATELY NOT CLAIMED. This is association, not survival analysis.
# A 28-day flag is tested as a two-group comparison, which is what it is; there
# is no time-to-event model here, because the sheet carries no follow-up time and
# a Cox model on a binary column with no time is a category error. If follow-up
# times are recorded later, that is the point to add one.

#' Is a vector usable as a numeric clinical variable?
#' @param v vector.
#' @keywords internal
clin_is_numeric <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(FALSE)
  if (is.numeric(v)) return(TRUE)
  suppressWarnings(all(!is.na(as.numeric(as.character(v)))))
}

#' Cliff's delta, the rank effect size for two groups
#'
#' Reported beside every two-group p-value because a p-value on nine samples
#' says mostly how many samples there were.
#' @param a,b numeric vectors.
#' @keywords internal
clin_cliff <- function(a, b) {
  if (!length(a) || !length(b)) return(NA_real_)
  mean(outer(b, a, ">")) - mean(outer(b, a, "<"))
}

#' Epsilon-squared, the rank effect size for Kruskal-Wallis
#'
#' WHY IT WAS ADDED. A three-level variable such as infection focus produced a
#' p-value and nothing else, so it appeared in the heatmap as a grey tile: the
#' reader could see that a test had been run and not whether it had found
#' anything. Epsilon-squared is `H / (n - 1)`, the standard rank effect size for
#' this test (Tomczak & Tomczak 2014), and is the proportion of rank variance the
#' grouping accounts for. It is bounded 0 to 1 and UNSIGNED, because a variable
#' with three levels has no single direction -- which is why it is carried in its
#' own column and never coloured on the signed scale.
#' @param h Kruskal-Wallis statistic.
#' @param n total number of observations.
#' @keywords internal
clin_epsilon2 <- function(h, n) {
  if (!is.finite(h) || !is.finite(n) || n < 2L) return(NA_real_)
  min(1, max(0, h / (n - 1)))
}

#' Percentile bootstrap interval for a rank effect size
#'
#' WHY AN INTERVAL AND NOT ONLY A POINT ESTIMATE. On nine samples a rho of 0.61
#' and a rho of 0.05 can be the same underlying quantity, and the point estimate
#' alone cannot say so. The interval is the part of the answer that carries how
#' little the cohort constrains the effect, and reporting the effect with its
#' interval beside the p-value is the standard recommendation for exactly this
#' situation.
#'
#' WHAT IT IS NOT. A percentile bootstrap at this sample size is itself
#' approximate: below about nine observations the effective resample is smaller
#' than the nominal one and coverage falls short of 95%. The interval is written
#' anyway, because a visibly wide interval is a better description of the
#' evidence than a bare point estimate, and it is suppressed below `min_n` rather
#' than being quoted at a width nobody should trust.
#'
#' The RNG stream is saved and restored, so the interval is reproducible for a
#' given table no matter what sampled before it, and nothing downstream of here
#' sees a different stream than it would have.
#' @param stat function of an index vector returning the effect.
#' @param n number of observations to resample.
#' @param B resamples.
#' @param min_n fewest observations for an interval to be reported at all.
#' @param seed fixed so the same table always yields the same interval.
#' @keywords internal
clin_boot_ci <- function(stat, n, B = 2000L, min_n = 6L, seed = 42L) {
  if (n < min_n) return(c(NA_real_, NA_real_))
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  v <- vapply(seq_len(B), function(i) {
    r <- try(stat(sample.int(n, n, replace = TRUE)), silent = TRUE)
    if (inherits(r, "try-error") || !is.finite(r)) NA_real_ else r
  }, numeric(1))
  v <- v[is.finite(v)]
  # A resample that degenerates -- every draw the same patient, so no variance to
  # correlate -- is dropped rather than counted as zero. Too few survivors means
  # the statistic is not estimable here at all, and no interval is claimed.
  if (length(v) < B * 0.5) return(c(NA_real_, NA_real_))
  unname(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))
}

#' Wrap a caption to the canvas it will be drawn on
#'
#' WHY IT TAKES THE WIDTH. A caption is drawn at one text size with no wrapping of
#' its own, so a line longer than the canvas is silently CLIPPED -- the figure
#' still writes and the end of the sentence is gone. A fixed wrap width therefore
#' only works for a fixed canvas, and several of these figures size themselves from
#' the data. Roughly seventeen characters fit per inch at the caption's size, with
#' the constant kept slightly conservative because the exact figure depends on the
#' device's font metrics.
#' @param txt one or more strings, pasted together.
#' @param width_in the canvas width in inches, the same value passed to ggsave.
#' @keywords internal
cap_wrap <- function(txt, width_in) {
  n <- max(40L, floor(width_in * 17))
  paste(strwrap(paste(txt, collapse = " "), width = n), collapse = "\n")
}

#' Well-separated colours for a small set of categories
#'
#' WHY NOT [population_colours()]. That spreads the population palette by index,
#' and the population palette holds two greens: at three levels it returns red,
#' green and spring-green, whose last two pair is the closest in the whole
#' palette. `study_palette` is the set chosen for exactly this job -- few
#' categories, each about a sixth of the hue wheel from the last -- so batches
#' and clinical strata come out as distinguishable as cohorts already do.
#' @param levels category labels.
#' @param colors palette.
#' @keywords internal
clin_cat_palette <- function(levels, colors = fcs_colors()) {
  levels <- as.character(levels)
  base <- colors$study_palette
  # TWO LEVELS GET BLUE AND AMBER, NOT GREEN AND RED. The study palette leads with
  # green then red, and levels arrive alphabetically, so a survival flag came out
  # "no" green and "yes" red -- a value judgement the wrong way round, asserted by
  # the palette rather than by anything in the data. Blue and amber carry no
  # good/bad reading, separate cleanly under red-green colour vision deficiency,
  # and leave green and red to mean reference-and-study where they already do.
  if (length(levels) == 2L)
    return(stats::setNames(c("#0072F0", "#FFA400"), levels))
  cols <- if (length(levels) <= length(base)) base[seq_along(levels)]
          else c(base, grDevices::hcl.colors(length(levels) - length(base),
                                             "Dark 3"))
  stats::setNames(cols, levels)
}

#' Associate per-sample values with clinical variables
#'
#' @param d long data frame with `sample_id`, a key column and a value column.
#' @param key_col column naming the thing measured: "population" or "marker".
#' @param value_col the per-sample quantity.
#' @param clin named list: variable -> named vector of sample_id -> value.
#' @param min_n fewest samples for a test to run at all.
#' @return data frame, one row per key x variable, BH-adjusted within variable.
#' @keywords internal
clin_associate <- function(d, key_col, value_col, clin, min_n = 4L) {
  if (is.null(d) || !nrow(d) || !length(clin)) return(NULL)
  if (!all(c("sample_id", key_col, value_col) %in% names(d))) return(NULL)
  d <- d[is.finite(d[[value_col]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  rows <- list()
  for (cv in names(clin)) {
    v_all <- clin[[cv]]
    numeric_cv <- clin_is_numeric(v_all)
    for (k in sort(unique(d[[key_col]]))) {
      dk <- d[d[[key_col]] == k, , drop = FALSE]
      dk <- dk[!duplicated(dk$sample_id), , drop = FALSE]
      dk$cv <- v_all[dk$sample_id]
      dk <- dk[!is.na(dk$cv), , drop = FALSE]
      if (nrow(dk) < min_n) next
      y <- dk[[value_col]]

      if (numeric_cv) {
        x <- suppressWarnings(as.numeric(as.character(dk$cv)))
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < min_n || length(unique(x[ok])) < 3L) next
        ct <- try(stats::cor.test(x[ok], y[ok], method = "spearman",
                                  exact = FALSE), silent = TRUE)
        if (inherits(ct, "try-error")) next
        # Resample the PAIRS, not the two vectors independently: the association
        # is a property of which value went with which patient, and shuffling
        # that away would bootstrap the null instead of the estimate.
        xo <- x[ok]; yo <- y[ok]
        ci <- clin_boot_ci(function(i) suppressWarnings(
                stats::cor(xo[i], yo[i], method = "spearman")), length(xo))
        rows[[length(rows) + 1L]] <- data.frame(
          variable = cv, key = k, kind = "numeric", test = "Spearman",
          n = sum(ok), estimate = unname(round(ct$estimate, 4)),
          effect = "rho", signed = TRUE,
          ci_low = round(ci[1], 4), ci_high = round(ci[2], 4),
          p_value = ct$p.value,
          levels_or_range = sprintf("%.4g to %.4g", min(x[ok]), max(x[ok])),
          stringsAsFactors = FALSE)
      } else {
        f <- factor(trimws(as.character(dk$cv)))
        if (nlevels(f) < 2L || min(table(f)) < 2L) next
        if (nlevels(f) == 2L) {
          a <- y[f == levels(f)[1]]; b <- y[f == levels(f)[2]]
          wt <- try(stats::wilcox.test(b, a, exact = FALSE), silent = TRUE)
          if (inherits(wt, "try-error")) next
          # Stratified resample: within each arm separately, so every resample
          # keeps both groups non-empty and the group sizes it was asked about.
          ci <- clin_boot_ci(function(i) {
            ia <- sample.int(length(a), length(a), replace = TRUE)
            ib <- sample.int(length(b), length(b), replace = TRUE)
            clin_cliff(a[ia], b[ib])
          }, length(y))
          rows[[length(rows) + 1L]] <- data.frame(
            variable = cv, key = k, kind = "two-level",
            test = "Wilcoxon rank-sum", n = length(y),
            estimate = round(clin_cliff(a, b), 4), effect = "Cliff's delta",
            signed = TRUE,
            ci_low = round(ci[1], 4), ci_high = round(ci[2], 4),
            p_value = wt$p.value,
            levels_or_range = paste0(levels(f)[1], " (n=", sum(f == levels(f)[1]),
                                     ") vs ", levels(f)[2], " (n=",
                                     sum(f == levels(f)[2]), ")"),
            stringsAsFactors = FALSE)
        } else {
          kw <- try(stats::kruskal.test(y, f), silent = TRUE)
          if (inherits(kw, "try-error")) next
          # UNSIGNED, so no interval is drawn on the signed scale and `signed` is
          # FALSE: three levels have a magnitude but no direction.
          rows[[length(rows) + 1L]] <- data.frame(
            variable = cv, key = k, kind = "multi-level",
            test = "Kruskal-Wallis", n = length(y),
            estimate = round(clin_epsilon2(unname(kw$statistic), length(y)), 4),
            effect = "epsilon squared", signed = FALSE,
            ci_low = NA_real_, ci_high = NA_real_, p_value = kw$p.value,
            levels_or_range = paste(levels(f), collapse = " / "),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  names(out)[names(out) == "key"] <- key_col
  # BH WITHIN a variable, across the populations tested against it. Each clinical
  # variable is its own question asked of every population, which is the family;
  # pooling across variables would penalise a well-powered variable for the
  # company it keeps.
  out$p_adj_BH <- stats::ave(out$p_value, out$variable,
                             FUN = function(p) stats::p.adjust(p, method = "BH"))
  out$significant_BH <- out$p_adj_BH < 0.05
  # An honest flag rather than a silent caveat: rank tests on this many samples
  # detect only very large effects, so a null result is uninformative.
  out$underpowered <- out$n < 10L
  out[order(out$variable, out$p_value), , drop = FALSE]
}

#' Clinical variables against population frequencies and marker intensities
#'
#' @param freq population_frequencies table.
#' @param mfi population_marker_mfi table, optional.
#' @param clin named list of clinical variables, sample_id -> value.
#' @param value_col the frequency column to use.
#' @param min_n fewest samples for a test.
#' @return list(populations, markers)
#' @export
stats_clinical_association <- function(freq, mfi = NULL, clin = list(),
                                       value_col = NULL, min_n = 4L) {
  if (!length(clin)) return(NULL)
  pop <- NULL
  if (!is.null(freq) && nrow(freq)) {
    if (is.null(value_col)) value_col <- abundance_measure(freq)$col
    if (value_col %in% names(freq))
      pop <- clin_associate(qc_pass_rows(freq), "population", value_col,
                            clin, min_n = min_n)
  }
  mk <- NULL
  if (!is.null(mfi) && nrow(mfi) && "median_asinh" %in% names(mfi)) {
    d <- qc_pass_rows(mfi)
    # One row per sample x marker x population would test the same marker many
    # times over. The question here is the marker across the whole sample, so the
    # populations are collapsed to one median per sample and marker.
    if (all(c("sample_id", "marker") %in% names(d))) {
      agg <- stats::aggregate(list(median_asinh = d$median_asinh),
                              by = list(sample_id = d$sample_id, marker = d$marker),
                              FUN = stats::median, na.rm = TRUE)
      mk <- clin_associate(agg, "marker", "median_asinh", clin, min_n = min_n)
    }
  }
  if (is.null(pop) && is.null(mk)) return(NULL)
  list(populations = pop, markers = mk)
}

#' Heatmap of clinical association across populations or markers
#'
#' WHAT IS ON THE TILE. The signed effect where one exists -- Spearman's rho, or
#' Cliff's delta -- so the colour carries direction as well as strength. A
#' Kruskal-Wallis row has no signed effect and is drawn grey with its p-value,
#' because a three-level variable has no single direction to show.
#'
#' Significance is marked on the tile rather than encoded in the colour, so a
#' large effect that did not survive correction still reads as large. On a cohort
#' this size that distinction is most of the message.
#' @param assoc data frame from [clin_associate()].
#' @param key_col "population" or "marker".
#' @param outfile path.
#' @param title figure title.
#' @param dpi resolution.
#' @param colors palette.
#' @export
fig_clinical_heatmap <- function(assoc, key_col, outfile,
                                 title = "Clinical association",
                                 dpi = 200, colors = fcs_colors()) {
  if (is.null(assoc) || !nrow(assoc) || !key_col %in% names(assoc))
    return(invisible(NULL))
  d <- assoc
  d$.key <- factor(d[[key_col]], levels = sort(unique(d[[key_col]])))
  d$.var <- factor(d$variable, levels = sort(unique(d$variable)))
  # SIGNED IS A COLUMN, NOT "the estimate happens to be finite". A multi-level
  # variable now carries epsilon-squared, which is finite and always positive, so
  # testing finiteness would have coloured every infection-focus tile red as
  # though the association had a direction.
  sgn <- if ("signed" %in% names(d)) !is.na(d$signed) & d$signed
         else is.finite(d$estimate)
  d$.eff <- ifelse(sgn & is.finite(d$estimate), d$estimate, NA_real_)
  # Built from code points rather than typed: a source file this package ships
  # has to stay ASCII for R CMD check. 03B5 is epsilon, 00B2 superscript two.
  .eps2 <- paste0(intToUtf8(c(0x03B5L, 0x00B2L)), " %.2f")
  d$.lab <- ifelse(sgn & is.finite(d$estimate), sprintf("%.2f", d$estimate),
                   ifelse(is.finite(d$estimate),
                          sprintf(.eps2, d$estimate),
                          sprintf("p=%.3g", d$p_value)))
  d$.mark <- ifelse(d$significant_BH, "*", "")

  fig <- ggplot(d, aes(.var, .key, fill = .eff)) +
    geom_tile(colour = colors$tile_border, linewidth = 0.3) +
    geom_text(aes(label = paste0(.lab, .mark)), size = 2.6,
              colour = colors$label_text) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, na.value = colors$na_fill,
                         name = "rho / delta", limits = c(-1, 1)) +
    labs(title = title,
         subtitle = paste("signed effect per tile; * marks BH-adjusted p < 0.05",
                          "within a variable. Grey = no signed effect",
                          "(Kruskal-Wallis), labelled with its p-value"),
         x = NULL, y = NULL,
         caption = paste("Rank tests on per-sample values, so the replicates are",
                         "subjects and not cells. Read the effect before the",
                         "asterisk: on a small cohort a large effect that misses",
                         "correction is the more useful signal.")) +
    theme_cyto(9, colors = colors) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "top", legend.justification = "left",
          legend.direction = "horizontal",
          plot.caption = element_text(size = 7, hjust = 0,
                                      colour = colors$caption_text))
  safe_ggsave(outfile, plot = fig,
              width = max(7, 1.1 * nlevels(d$.var) + 3.4),
              height = max(4, 0.32 * nlevels(d$.key) + 2.2),
              dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' Per-population detail against one clinical variable
#'
#' The heatmap says which associations are worth looking at; this shows the data
#' behind them, because a rank correlation of 0.8 on nine points can be one
#' outlier and the scatter is the only place that is visible.
#' @param freq population_frequencies table.
#' @param clin_values named vector sample_id -> value.
#' @param var_name the variable's name.
#' @param outfile path.
#' @param value_col frequency column.
#' @param ncol panels per row.
#' @param dpi resolution.
#' @param colors palette.
#' @export
fig_clinical_detail <- function(freq, clin_values, var_name, outfile,
                                value_col = NULL, ncol = 4L, dpi = 200,
                                colors = fcs_colors()) {
  if (is.null(freq) || !nrow(freq) || !length(clin_values)) return(invisible(NULL))
  if (is.null(value_col)) value_col <- abundance_measure(freq)$col
  if (!value_col %in% names(freq)) return(invisible(NULL))
  d <- qc_pass_rows(freq)
  d <- d[!duplicated(d[, c("sample_id", "population")]), , drop = FALSE]
  d$cv <- clin_values[d$sample_id]
  d <- d[!is.na(d$cv) & is.finite(d[[value_col]]), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  d$.y <- d[[value_col]]
  npop <- length(unique(d$population))

  if (clin_is_numeric(d$cv)) {
    d$.x <- suppressWarnings(as.numeric(as.character(d$cv)))
    d <- d[is.finite(d$.x), , drop = FALSE]
    if (!nrow(d)) return(invisible(NULL))
    fig <- ggplot(d, aes(.x, .y)) +
      geom_point(size = 1.5, alpha = 0.85, colour = colors$accent %||% "#0a7d4a") +
      # A GUIDE FOR THE EYE, AND THE SUBTITLE SAYS SO. The test was Spearman on
      # ranks and this line is least squares on the raw values, so its slope is
      # not the tested quantity: it is here to make the direction of a cloud of
      # nine points readable, and a reader who takes its steepness for the effect
      # size has been told otherwise on the figure itself.
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5,
                  colour = colors$threshold_review %||% "firebrick") +
      facet_wrap(~ population, scales = "free_y", ncol = ncol) +
      labs(title = paste0("Population abundance against ", var_name),
           subtitle = paste("one point per sample; line is a least-squares guide,",
                            "the test was Spearman on ranks"),
           x = var_name, y = value_col)
  } else {
    d$.x <- factor(trimws(as.character(d$cv)))
    fig <- ggplot(d, aes(.x, .y)) +
      geom_boxplot(outlier.shape = NA, fill = colors$na_fill %||% "grey85",
                   colour = colors$bar_outline %||% "grey20", linewidth = 0.3) +
      geom_jitter(width = 0.15, height = 0, size = 1.4, alpha = 0.85,
                  colour = colors$accent %||% "#0a7d4a") +
      facet_wrap(~ population, scales = "free_y", ncol = ncol) +
      labs(title = paste0("Population abundance by ", var_name),
           subtitle = "one point per sample; box is median and quartiles",
           x = var_name, y = value_col)
  }
  fig <- fig + theme_cyto(8, colors = colors) +
    theme(panel.border = element_rect(colour = colors$axis_ticks, fill = NA,
                                      linewidth = 0.4),
          panel.spacing = grid::unit(0.6, "lines"),
          strip.background = element_rect(fill = colors$grid_major, colour = NA),
          strip.text = element_text(size = 7))
  nr <- ceiling(npop / ncol)
  safe_ggsave(outfile, plot = fig, width = 2.6 * min(ncol, npop) + 1.2,
              height = 2.3 * nr + 1.1, dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

# -----------------------------------------------------------------------------
# HOW THESE FIGURES WERE CHOSEN
#
# The heatmap and the per-variable detail panels answer "which pairs moved
# together" and "what do the points look like". Three questions a reader of a
# small clinical cohort asks next had no figure at all:
#
#   How well is each effect pinned down?     -> fig_clinical_forest()
#     The ordered effect with its bootstrap interval. On nine patients the
#     interval usually spans zero, and seeing that is the finding: a rho of 0.61
#     whose interval runs -0.1 to 0.9 is a lead, not a result. Reporting the
#     effect with an interval beside the p-value rather than the p-value alone is
#     the standard recommendation for this situation.
#
#   Do the variables agree with each other?  -> fig_clinical_correlogram()
#     p-values are adjusted WITHIN each variable, which treats the variables as
#     separate questions. When SOFA and survival are themselves correlated they
#     are not separate questions, and the same finding is being counted twice.
#     Nothing else in the run would show that.
#
#   Where does one patient sit in all of it? -> fig_clinical_landscape()
#     One column per sample, ordered by the clinical variable, with every
#     clinical variable drawn as a strip above the populations. This is the
#     standard immunoprofiling layout -- a z-scored matrix with clinical
#     annotation bars -- and it is the only figure here that shows severity,
#     outcome, timepoint, infection focus and the cell populations at once, which
#     is how a gradient across the cohort becomes visible rather than inferred
#     from a table of correlations.
#
#   And where a cohort was sampled repeatedly:  fig_clinical_trajectory()
#     Per-patient lines across timepoints, split by outcome. In the sepsis
#     literature this is the figure that separates survivors from non-survivors:
#     the difference is in the direction of travel, not in any single timepoint,
#     and a boxplot per timepoint averages exactly that away.
# -----------------------------------------------------------------------------

#' Ordered effect sizes with bootstrap intervals for one clinical variable
#'
#' WHAT IS ON IT. One row per population (or marker), ordered by effect: the
#' point is Spearman's rho or Cliff's delta, the bar is its 95% percentile
#' bootstrap interval, and the dashed line is no effect. Rows whose interval
#' clears the line are the ones a follow-up cohort should be powered on.
#'
#' A three-level variable is drawn on epsilon-squared instead, which is unsigned,
#' so its axis starts at zero and no interval is drawn: there is no direction to
#' put an interval around.
#'
#' WHY IT IS ORDERED BY EFFECT AND NOT BY p. At these sample sizes the p-value
#' ordering is close to the effect ordering but not identical, and where they
#' disagree the effect is the one worth reading: p carries how many samples there
#' were as much as how large the difference is.
#' @param assoc data frame from [stats_clinical_association()].
#' @param key_col "population" or "marker".
#' @param var_name which clinical variable to draw.
#' @param outfile path.
#' @param dpi resolution.
#' @param colors palette.
#' @export
fig_clinical_forest <- function(assoc, key_col, var_name, outfile, dpi = 200,
                                colors = fcs_colors()) {
  if (is.null(assoc) || !nrow(assoc) || !key_col %in% names(assoc))
    return(invisible(NULL))
  d <- assoc[assoc$variable == var_name & is.finite(assoc$estimate), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  sgn <- if ("signed" %in% names(d)) all(!is.na(d$signed) & d$signed) else TRUE
  d$.key <- factor(d[[key_col]], levels = d[[key_col]][order(d$estimate)])
  d$.lo <- if ("ci_low"  %in% names(d)) d$ci_low  else NA_real_
  d$.hi <- if ("ci_high" %in% names(d)) d$ci_high else NA_real_
  # Three states where there is an interval to speak of, two where there is not.
  # "Survived correction" and "interval clears zero" are different claims and both
  # are weaker than they look at n < 12, so the strongest thing the figure asserts
  # is which rows are worth a second look. An unsigned effect has no interval and
  # no zero to clear, so labelling its rows "not distinguished from zero" would be
  # asserting the result of a test that was never run.
  d$.sig <- if (sgn)
    ifelse(!is.na(d$significant_BH) & d$significant_BH, "BH-adjusted p < 0.05",
           ifelse(is.finite(d$.lo) & is.finite(d$.hi) & (d$.lo > 0 | d$.hi < 0),
                  "interval excludes zero", "not distinguished from zero"))
  else
    ifelse(!is.na(d$significant_BH) & d$significant_BH, "BH-adjusted p < 0.05",
           "not significant after correction")
  key_cols <- stats::setNames(
    c(colors$gate_highlight %||% "#D62728", colors$study_palette[3] %||% "#8E44E8",
      colors$axis_ticks %||% "grey30", colors$axis_ticks %||% "grey30"),
    c("BH-adjusted p < 0.05", "interval excludes zero",
      "not distinguished from zero", "not significant after correction"))
  eff_lab <- if (sgn) unique(d$effect)[1] else "epsilon squared (unsigned)"

  fig <- ggplot(d, aes(estimate, .key, colour = .sig))
  if (sgn)
    fig <- fig +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4,
                 colour = colors$axis_ticks) +
      # geom_linerange rather than geom_errorbarh, which ggplot2 soft-deprecated
      # in 3.5.0; the caps it would add carry no information here anyway.
      geom_linerange(aes(xmin = .lo, xmax = .hi), linewidth = 0.55,
                     na.rm = TRUE)
  fig <- fig +
    geom_point(size = 2.1) +
    # The sample count sits INSIDE the panel on the right. At Inf with the default
    # expansion it landed in the plot margin, outside the panel border, which reads
    # as a stray annotation rather than as a column of the figure.
    geom_text(aes(label = paste0("n=", n)), x = Inf, hjust = 1.15, size = 2.3,
              colour = colors$label_text, show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.16))) +
    scale_colour_manual(values = key_cols, name = NULL, drop = TRUE) +
    labs(title = paste0(tools::toTitleCase(key_col), " effects against ", var_name),
         subtitle = if (sgn)
           paste("point is the effect, bar is its 95% percentile bootstrap",
                 "interval, dashed line is no effect")
         else
           paste("unsigned effect: a variable with three or more levels has a",
                 "magnitude but no direction, so no interval is drawn"),
         x = eff_lab, y = NULL,
         caption = cap_wrap(if (sgn) c(
           "Rank tests on per-sample values; the replicates are subjects.",
           "The interval is a percentile bootstrap and is itself approximate",
           "below about nine samples, where the effective resample is smaller",
           "than the nominal one -- read a wide interval as 'this cohort does",
           "not constrain the effect', which is the honest reading, rather than",
           "as a precise range.")
           else c(
           "Rank tests on per-sample values; the replicates are subjects.",
           "Epsilon-squared is the proportion of rank variance the grouping",
           "accounts for: H / (n - 1), bounded 0 to 1. There is no interval",
           "because there is no direction to put one around, and the ordering",
           "here says which populations differ most between the levels, not",
           "which level is higher. The box plot for this variable is where",
           "that is visible."), 8.4)) +
    theme_cyto(9, colors = colors) +
    theme(legend.position = "top", legend.justification = "left",
          legend.direction = "horizontal",
          panel.grid.major.y = element_line(colour = colors$grid_major,
                                            linewidth = 0.3),
          plot.caption = element_text(size = 7, hjust = 0,
                                      colour = colors$caption_text))
  safe_ggsave(outfile, plot = fig, width = 8.4,
              height = max(3.4, 0.26 * nlevels(d$.key) + 2.4), dpi = dpi,
              limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' Spearman correlation between the clinical variables themselves
#'
#' WHY THIS IS NOT AN AFTERTHOUGHT. [stats_clinical_association()] adjusts
#' p-values within each variable, on the grounds that each variable is its own
#' question. That is only true when the variables carry different information. A
#' cohort where the sickest patients are also the ones who died has one gradient
#' and two columns describing it, and an association found against both is one
#' finding reported twice. This figure is where that is visible.
#'
#' WHAT IS INCLUDED. Numeric variables, and two-level variables coded 0/1 -- for
#' which Spearman is the rank-biserial correlation, a legitimate quantity.
#' Variables with three or more unordered levels are excluded and named in the
#' caption: there is no ordering to correlate, and coding them 1/2/3 would invent
#' one.
#' @param clin named list of clinical variables, sample_id -> value.
#' @param outfile path.
#' @param dpi resolution.
#' @param colors palette.
#' @export
fig_clinical_correlogram <- function(clin, outfile, dpi = 200,
                                     colors = fcs_colors()) {
  if (length(clin) < 2L) return(invisible(NULL))
  ids <- sort(unique(unlist(lapply(clin, names))))
  if (length(ids) < 4L) return(invisible(NULL))
  num <- list(); dropped <- character(0)
  for (cv in names(clin)) {
    v <- clin[[cv]][ids]
    if (clin_is_numeric(v)) {
      num[[cv]] <- suppressWarnings(as.numeric(as.character(v)))
    } else {
      f <- factor(trimws(as.character(v)))
      if (nlevels(f) == 2L) {
        # 0/1 for a two-level variable, so the coefficient is the rank-biserial
        # correlation. The level that sorts first is 0, which is stated in the
        # caption because the SIGN is meaningless without it.
        num[[paste0(cv, " [", levels(f)[2], "=1]")]] <- as.numeric(f) - 1
      } else dropped <- c(dropped, cv)
    }
  }
  if (length(num) < 2L) return(invisible(NULL))
  m <- as.matrix(as.data.frame(num, check.names = FALSE))
  nv <- ncol(m)
  rows <- list()
  for (i in seq_len(nv)) for (j in seq_len(nv)) {
    ok <- is.finite(m[, i]) & is.finite(m[, j])
    r <- NA_real_; p <- NA_real_
    if (sum(ok) >= 4L && length(unique(m[ok, i])) > 1L &&
        length(unique(m[ok, j])) > 1L) {
      ct <- try(stats::cor.test(m[ok, i], m[ok, j], method = "spearman",
                                exact = FALSE), silent = TRUE)
      if (!inherits(ct, "try-error")) { r <- unname(ct$estimate); p <- ct$p.value }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      .a = colnames(m)[i], .b = colnames(m)[j], .r = r, .p = p, .n = sum(ok),
      stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  lv <- colnames(m)
  d$.a <- factor(d$.a, levels = lv); d$.b <- factor(d$.b, levels = rev(lv))
  d$.absr <- ifelse(is.finite(d$.r), abs(d$.r), 0)
  diagonal <- as.character(d$.a) == as.character(d$.b)
  # The coefficient is printed only where it reached nominal significance, and the
  # circle is drawn either way: an unlabelled small circle says "measured, not
  # distinguishable from zero", which is different from a blank cell saying
  # "not measured". The diagonal is skipped -- a variable correlates with itself at
  # 1.00 by construction, and printing it there collided with the sample count,
  # which is the one thing the diagonal is carrying.
  d$.lab <- ifelse(!diagonal & is.finite(d$.p) & d$.p < 0.05,
                   sprintf("%.2f", d$.r), "")
  # The canvas is sized from the number of variables, so the caption has to be
  # wrapped to that width and not to a constant -- at three variables a 108-column
  # caption ran off the right edge and lost the end of every line.
  cell <- 0.7
  wd <- cell * nv + 3.1

  fig <- ggplot(d, aes(.a, .b)) +
    geom_point(aes(size = .absr, fill = .r), shape = 21, colour = colors$tile_border,
               stroke = 0.3, na.rm = TRUE) +
    geom_text(aes(label = .lab), size = 2.4, colour = colors$label_text) +
    geom_text(data = d[diagonal, , drop = FALSE],
              aes(label = paste0("n=", .n)), size = 2.2,
              colour = colors$label_text) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), na.value = colors$na_fill,
                         name = "Spearman rho") +
    # A FLOOR ON THE SIZE, not scale_size_area, which maps zero to zero: a rho of
    # 0.03 became an invisible dot, indistinguishable from a cell where the test
    # could not run. Every measured cell now shows something.
    scale_size(range = c(1.8, 11), limits = c(0, 1), guide = "none") +
    labs(title = "The clinical variables against each other",
         subtitle = paste("circle area is |rho|, fill is its sign; the number is",
                          "printed where p < 0.05 unadjusted"),
         x = NULL, y = NULL,
         caption = cap_wrap(paste0(
           "Read this BEFORE the association heatmap. p-values there are ",
           "adjusted within each variable, which assumes the variables are ",
           "separate questions; two variables strongly correlated here are one ",
           "question asked twice, and a population associated with both is one ",
           "finding, not two. A two-level variable is coded 0/1 with the level ",
           "named in its axis label as 1, so its sign is read against that.",
           if (length(dropped))
             paste0(" Excluded, having three or more unordered levels: ",
                    paste(dropped, collapse = ", "),
                    " -- there is no ordering to correlate.")
           else ""), wd)) +
    theme_cyto(9, colors = colors) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.major = element_line(colour = colors$grid_major,
                                          linewidth = 0.25),
          legend.position = "top", legend.justification = "left",
          legend.direction = "horizontal",
          plot.caption = element_text(size = 7, hjust = 0,
                                      colour = colors$caption_text))
  # A FIXED INCH BUDGET PER CELL rather than a square canvas. A correlation matrix
  # is read as a grid, so the cells have to be about as tall as they are wide, and
  # sizing the canvas by the number of variables is what achieves that. Sizing it
  # square instead left three rows spread over six inches, with the circles lost
  # in the middle of their cells.
  safe_ggsave(outfile, plot = fig, width = wd, height = cell * nv + 3.6,
              dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' The cohort as one picture: populations, samples and every clinical variable
#'
#' WHAT IT IS. One column per sample, ordered along the clinical variable named in
#' `order_by`. Below the strips, one row per population, filled by that sample's
#' value expressed as a z-score WITHIN the population, so a row is readable
#' whether the population is 40% of CD45+ or 0.4%. Above them, one strip per
#' clinical variable.
#'
#' WHY Z-SCORE WITHIN A ROW AND NOT ACROSS THE MATRIX. Populations differ by two
#' orders of magnitude in abundance. On a shared scale the frequent lineages
#' saturate the palette and every rare one is the same colour, so the figure shows
#' which populations are large -- which is already known -- instead of which
#' samples are unusual, which is the question.
#'
#' WHAT IT CANNOT SHOW. Nothing here is a test. A gradient that reads clearly
#' across twelve ordered columns can still be the ordering the eye was given, and
#' the tests in `clinical_association.csv` are what decide.
#' @param freq population_frequencies table.
#' @param clin named list of clinical variables, sample_id -> value.
#' @param outfile path.
#' @param order_by which clinical variable orders the columns; the first numeric
#'   one by default.
#' @param value_col frequency column.
#' @param max_vars most clinical strips to draw, so the key cannot outgrow the
#'   figure.
#' @param dpi resolution.
#' @param colors palette.
#' @export
fig_clinical_landscape <- function(freq, clin, outfile, order_by = NULL,
                                   value_col = NULL, max_vars = 6L, dpi = 200,
                                   colors = fcs_colors()) {
  if (is.null(freq) || !nrow(freq) || !length(clin)) return(invisible(NULL))
  if (is.null(value_col)) value_col <- abundance_measure(freq)$col
  if (!value_col %in% names(freq)) return(invisible(NULL))
  d <- qc_pass_rows(freq)
  d <- d[!duplicated(d[, c("sample_id", "population")]), , drop = FALSE]
  d <- d[is.finite(d[[value_col]]), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))

  vars <- utils::head(names(clin), max_vars)
  # Order by the first numeric variable, because ordering is only meaningful for
  # something with an order; with no numeric variable the columns stay in
  # sample-id order and the strips still carry the categories.
  if (is.null(order_by)) {
    nm <- vars[vapply(vars, function(v) clin_is_numeric(clin[[v]]), logical(1))]
    order_by <- if (length(nm)) nm[1] else NA_character_
  }
  ids <- sort(unique(d$sample_id))
  if (length(ids) < 3L) return(invisible(NULL))
  if (!is.na(order_by) && order_by %in% names(clin)) {
    ov <- suppressWarnings(as.numeric(as.character(clin[[order_by]][ids])))
    if (!all(is.na(ov))) ids <- ids[order(ov, na.last = TRUE)]
  }
  d$.samp <- factor(d$sample_id, levels = ids)
  d <- d[!is.na(d$.samp), , drop = FALSE]

  # Z within population. A population measured in one sample only has no spread
  # and is dropped rather than drawn as a row of zeros.
  d$.z <- stats::ave(d[[value_col]], d$population, FUN = function(v) {
    s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(NA_real_, length(v))
    else (v - mean(v, na.rm = TRUE)) / s
  })
  d <- d[is.finite(d$.z), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  # Most abundant population at the top, so the eye starts where the counting
  # uncertainty is smallest.
  ord <- stats::aggregate(stats::reformulate("population", value_col), d, stats::median)
  d$.pop <- factor(d$population, levels = ord$population[order(ord[[value_col]])])

  wd <- max(8, 0.34 * length(ids) + 4.6)
  main <- ggplot(d, aes(.samp, .pop, fill = .z)) +
    geom_tile(colour = colors$tile_border, linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, na.value = colors$na_fill,
                         name = "z within population") +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL,
         caption = cap_wrap(paste0(
           "Descriptive. Columns are samples", if (!is.na(order_by))
             paste0(" ordered by ", order_by) else "",
           "; the fill is a z-score within each row, so colours compare samples ",
           "within a population and never populations with each other. A visible ",
           "gradient is a reason to read clinical_association.csv, not a result ",
           "on its own -- twelve columns in a chosen order will often look like ",
           "a trend."), wd)) +
    theme_cyto(8, colors = colors) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
          panel.grid = element_blank(),
          plot.caption = element_text(size = 7, hjust = 0,
                                      colour = colors$caption_text))

  # ---- one strip per clinical variable, each with its own scale --------------
  # Separate plots rather than extra rows of the matrix: ggplot allows one fill
  # scale per plot, and a severity score and an infection focus cannot share one.
  strip <- function(cv) {
    v <- clin[[cv]][ids]
    s <- data.frame(.samp = factor(ids, levels = ids), .lab = cv,
                    stringsAsFactors = FALSE)
    base <- function(p) p +
      labs(x = NULL, y = NULL) +
      scale_x_discrete(drop = FALSE) +
      theme_cyto(8, colors = colors) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_rect(fill = colors$grid_major, colour = NA),
            strip.text.y.left = element_text(angle = 0, size = 6.5, hjust = 1),
            legend.key.size = grid::unit(0.3, "cm"),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 6.5)) +
      facet_grid(rows = vars(.lab), switch = "y")
    if (clin_is_numeric(v)) {
      s$.val <- suppressWarnings(as.numeric(as.character(v)))
      base(ggplot(s, aes(.samp, 1, fill = .val)) +
             geom_tile(colour = colors$tile_border, linewidth = 0.25) +
             # Trimmed at both ends: magma runs to near-black, and a black tile
             # for the highest severity score reads as a missing value, which is
             # what na.value is for.
             scale_fill_viridis_c(option = "magma", direction = -1,
                                  begin = 0.12, end = 0.94,
                                  na.value = colors$na_fill, name = cv))
    } else {
      s$.val <- factor(trimws(as.character(v)))
      base(ggplot(s, aes(.samp, 1, fill = .val)) +
             geom_tile(colour = colors$tile_border, linewidth = 0.25) +
             scale_fill_manual(values = clin_cat_palette(levels(s$.val),
                                                         colors = colors),
                               na.value = colors$na_fill, name = cv))
    }
  }
  strips <- lapply(vars, strip)
  npop <- nlevels(d$.pop)
  # GUIDES COLLECTED, ON THE RIGHT, AND NOT ON TOP. Every other figure in this
  # package puts its key above the panel, because there it costs nothing. Here
  # there are as many keys as there are clinical variables plus one, and each of
  # the strips is a plot of its own -- patchwork lays a collected guide out once
  # for the whole composition, so collecting them is also what keeps the strips'
  # columns aligned with the matrix below. Left uncollected, a strip with a wide
  # key would have a narrower panel than the matrix and its columns would no
  # longer stand over the samples they describe.
  fig <- patchwork::wrap_plots(
    c(strips, list(main)), ncol = 1, guides = "collect",
    heights = c(rep(0.34, length(strips)), max(2.6, 0.2 * npop))) +
    patchwork::plot_annotation(
      title = "The cohort in one picture: clinical variables above, populations below",
      subtitle = paste0("one column per sample",
                        if (!is.na(order_by)) paste0(", ordered by ", order_by) else "",
                        "; strips are the clinical variables, tiles are ",
                        "population abundance as a z-score within each row"),
      theme = theme_cyto(9, colors = colors) +
        theme(legend.position = "right", legend.key.size = grid::unit(0.32, "cm"),
              legend.text = element_text(size = 6.5),
              legend.title = element_text(size = 7)))
  safe_ggsave(outfile, plot = fig, width = wd,
              height = 0.235 * npop + 0.46 * length(strips) + 3.2,
              dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' Per-patient trajectories across timepoints, split by outcome
#'
#' WHY LINES AND NOT A BOX PER TIMEPOINT. In a repeated-measures cohort the signal
#' is usually the DIRECTION each patient moves, and a box per timepoint averages
#' exactly that away: two patients rising and two falling look like a flat cohort.
#' In the sepsis literature this is the figure that separates survivors from
#' non-survivors, because they differ in trajectory rather than at any one
#' timepoint.
#'
#' WHAT IT IS NOT. Descriptive. The test for a repeated-measures design is the
#' paired comparison in `paired_comparison_stats.csv`; a trajectory figure shows
#' the shape of the data that test was run on.
#' @param freq population_frequencies table.
#' @param time_of named vector sample_id -> timepoint.
#' @param patient_of named vector sample_id -> patient.
#' @param outfile path.
#' @param outcome_of optional named vector sample_id -> outcome, used to colour
#'   the lines and to split the mean.
#' @param outcome_name label for the outcome in the key.
#' @param value_col frequency column.
#' @param ncol panels per row.
#' @param time_levels order of the timepoints; sorted unique values by default.
#' @param dpi resolution.
#' @param colors palette.
#' @export
fig_clinical_trajectory <- function(freq, time_of, patient_of, outfile,
                                    outcome_of = NULL, outcome_name = "outcome",
                                    value_col = NULL, ncol = 4L,
                                    time_levels = NULL, dpi = 200,
                                    colors = fcs_colors()) {
  if (is.null(freq) || !nrow(freq) || !length(time_of) || !length(patient_of))
    return(invisible(NULL))
  if (is.null(value_col)) value_col <- abundance_measure(freq)$col
  if (!value_col %in% names(freq)) return(invisible(NULL))
  d <- qc_pass_rows(freq)
  d <- d[!duplicated(d[, c("sample_id", "population")]), , drop = FALSE]
  d$.t <- as.character(time_of[d$sample_id])
  d$.pid <- as.character(patient_of[d$sample_id])
  d$.y <- d[[value_col]]
  d <- d[!is.na(d$.t) & !is.na(d$.pid) & is.finite(d$.y), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  tl <- time_levels %||% sort(unique(d$.t))
  d$.t <- factor(d$.t, levels = tl)
  d <- d[!is.na(d$.t), , drop = FALSE]
  # A single timepoint is not a trajectory. Saying so and writing nothing is
  # better than writing a figure of unconnected points called a trajectory.
  if (nlevels(droplevels(d$.t)) < 2L) {
    log_msg("  NOTE trajectory figure skipped: only one timepoint resolved")
    return(invisible(NULL))
  }
  has_out <- !is.null(outcome_of) && length(outcome_of)
  if (has_out) {
    d$.out <- as.character(outcome_of[d$sample_id])
    d$.out[is.na(d$.out)] <- "unknown"
    d$.out <- factor(d$.out, levels = sort(unique(d$.out)))
  } else d$.out <- factor("all samples")
  # THE KEY CARRIES THE PATIENT COUNT. The thick line is a median, and a median
  # over one patient is that patient drawn heavier than their own thin line --
  # which asserts a cohort summary that does not exist. Putting n in the key is
  # what stops "Yes" reading as a group when it is one person. Patients, not
  # samples: the same patient contributes a point at every timepoint.
  .npat <- vapply(levels(d$.out), function(lv)
    length(unique(d$.pid[d$.out == lv])), integer(1))
  .out_lab <- sprintf("%s (%d patient%s)", levels(d$.out), .npat,
                      ifelse(.npat == 1L, "", "s"))
  npop <- length(unique(d$population))
  wd <- max(7, 2.4 * min(ncol, npop) + 1.2)

  # The cohort mean per timepoint drawn over the individual lines: the individual
  # lines are the evidence and the mean is the summary, so the summary is the one
  # drawn heavier and on top.
  mn <- stats::aggregate(list(.y = d$.y),
                         by = list(.t = d$.t, .out = d$.out,
                                   population = d$population),
                         FUN = stats::median, na.rm = TRUE)
  fig <- ggplot(d, aes(.t, .y)) +
    geom_line(aes(group = interaction(.pid, .out), colour = .out),
              linewidth = 0.35, alpha = 0.55) +
    geom_point(aes(colour = .out), size = 1.1, alpha = 0.8) +
    geom_line(data = mn, aes(group = .out, colour = .out), linewidth = 1.15) +
    scale_colour_manual(values = clin_cat_palette(levels(d$.out), colors = colors),
                        labels = .out_lab,
                        name = if (has_out) outcome_name else NULL) +
    facet_wrap(~ population, scales = "free_y", ncol = ncol) +
    labs(title = "Population abundance across timepoints, one line per patient",
         subtitle = paste0("thin lines are patients, thick line is the median at ",
                           "each timepoint", if (has_out)
                             paste0(" within each ", outcome_name) else ""),
         x = NULL, y = value_col,
         caption = cap_wrap(c(
           "Descriptive. The test for this design is the paired comparison in",
           "paired_comparison_stats.csv, which uses only patients present at",
           "both timepoints; a line that stops early is a patient who was not",
           "sampled again, and reading the median at the last timepoint as the",
           "cohort's endpoint would be reading the survivors of the sampling",
           "rather than of the illness."), wd)) +
    theme_cyto(8, colors = colors) +
    theme(panel.border = element_rect(colour = colors$axis_ticks, fill = NA,
                                      linewidth = 0.4),
          panel.spacing = grid::unit(0.6, "lines"),
          strip.background = element_rect(fill = colors$grid_major, colour = NA),
          strip.text = element_text(size = 7),
          legend.position = "top", legend.justification = "left",
          legend.direction = "horizontal",
          plot.caption = element_text(size = 7, hjust = 0,
                                      colour = colors$caption_text))
  nr <- ceiling(npop / ncol)
  safe_ggsave(outfile, plot = fig, width = wd,
              height = 2.2 * nr + 1.6, dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}

#' Population abundance per batch
#'
#' WHY THIS EXISTS SEPARATELY from the batch diagnostics. Those ask whether the
#' EMBEDDING separates by batch, which is a question about the cells. This asks
#' whether the reported numbers differ by batch, which is a question about the
#' result: a population that shifts with acquisition batch is a population whose
#' between-group difference may be an acquisition difference.
#' @param freq population_frequencies table.
#' @param batch_of named vector sample_id -> batch.
#' @param outfile path.
#' @param value_col frequency column.
#' @param ncol panels per row.
#' @param dpi resolution.
#' @param colors palette.
#' @export
fig_populations_by_batch <- function(freq, batch_of, outfile, value_col = NULL,
                                     ncol = 4L, dpi = 200, colors = fcs_colors()) {
  if (is.null(freq) || !nrow(freq) || !length(batch_of)) return(invisible(NULL))
  if (is.null(value_col)) value_col <- abundance_measure(freq)$col
  if (!value_col %in% names(freq)) return(invisible(NULL))
  d <- qc_pass_rows(freq)
  d <- d[!duplicated(d[, c("sample_id", "population")]), , drop = FALSE]
  d$batch <- batch_of[d$sample_id]
  d <- d[!is.na(d$batch) & is.finite(d[[value_col]]), , drop = FALSE]
  if (!nrow(d) || length(unique(d$batch)) < 2L) return(invisible(NULL))
  d$.y <- d[[value_col]]
  d$batch <- factor(d$batch, levels = sort(unique(d$batch)))
  npop <- length(unique(d$population))
  nb <- nlevels(d$batch)
  wd <- max(7, 2.2 * min(ncol, npop) + 1.2)

  fig <- ggplot(d, aes(batch, .y, fill = batch)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3,
                 colour = colors$bar_outline %||% "grey20") +
    geom_jitter(width = 0.15, height = 0, size = 1.2, alpha = 0.85,
                colour = colors$label_text %||% "grey25") +
    # A KEY, NOT ONLY AXIS TEXT. facet_wrap draws x-axis labels on the bottom
    # panel of each column and nowhere else, so on a three-row grid two thirds of
    # the panels had coloured boxes with nothing naming them: the batch could only
    # be recovered by counting across from a panel three rows down. The key names
    # every batch once, above the panels, and the axis text stays for the panels
    # that have it.
    scale_fill_manual(values = clin_cat_palette(levels(d$batch), colors = colors),
                      name = "Batch") +
    facet_wrap(~ population, scales = "free_y", ncol = ncol) +
    labs(title = "Population abundance by acquisition batch",
         subtitle = paste0(nb, " batches; one point per sample. A population that ",
                           "steps with batch is one whose group difference may be ",
                           "an acquisition difference"),
         x = NULL, y = value_col,
         caption = cap_wrap(c(
           "Descriptive. batch_group_confounding.csv reports whether batch and",
           "study group can be told apart at all; where they cannot, no",
           "correction here or elsewhere can separate them."), wd)) +
    theme_cyto(8, colors = colors) +
    theme(panel.border = element_rect(colour = colors$axis_ticks, fill = NA,
                                      linewidth = 0.4),
          panel.spacing = grid::unit(0.6, "lines"),
          strip.background = element_rect(fill = colors$grid_major, colour = NA),
          strip.text = element_text(size = 7),
          legend.position = "top", legend.justification = "left",
          legend.direction = "horizontal",
          plot.caption = element_text(size = 7, hjust = 0,
                                      colour = colors$caption_text))
  nr <- ceiling(npop / ncol)
  safe_ggsave(outfile, plot = fig, width = wd,
              height = 2.3 * nr + 1.55, dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile)
  invisible(fig)
}
