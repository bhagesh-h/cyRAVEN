# =============================================================================
# SECTION 9c -- PARAMETRIC TESTS AND POST-HOC COMPARISONS
#
# WHY THESE EXIST ALONGSIDE THE RANK TESTS RATHER THAN INSTEAD OF THEM.
#
# The rank tests stay primary. They assume least, they do not break when one
# donor is extreme, and at the sample sizes an immunophenotyping cohort actually
# has they lose little against their parametric counterparts.
#
# What they cannot do is answer a reviewer who asks for a t-test or an ANOVA,
# which is still what most immunology papers report. Refusing to compute one and
# refusing to report it are different things: the first leaves the user to run it
# somewhere else, by hand, on exported numbers, with no record of whether its
# assumptions held.
#
# So the parametric result is computed and written to its own file, with the two
# assumption checks attached to every row. A reader can then see that a t-test
# was possible, or see exactly which assumption failed. Nothing here changes
# group_comparison_stats.csv.
#
# THE PROPORTION PROBLEM. Population frequencies are proportions. They are bound
# to [0, 100], their variance depends on their mean, and near 0 or 100 they are
# strongly skewed. A t-test on raw percentages violates its own assumptions
# hardest exactly where rare populations live. The arcsine square root transform
# stabilises that variance and is the standard remedy for proportion data, so the
# parametric tests run on transformed values and the table records that they did.
# The effect sizes are reported on both scales, because a difference in arcsine
# units is not something anyone can interpret.
# =============================================================================

#' Arcsine square root transform for percentages
#'
#' Maps [0, 100] onto [0, pi/2]. Variance-stabilising for proportion data.
#'
#' @param pct Numeric vector of percentages.
#' @return Numeric vector on the arcsine scale.
#' @keywords internal
asin_sqrt_pct <- function(pct) {
  p <- pmin(pmax(as.numeric(pct) / 100, 0), 1)
  asin(sqrt(p))
}

#' Cohen's d, pooled standard deviation
#' @param a,b Numeric vectors.
#' @return Numeric, or NA when either group has fewer than two finite values.
#' @keywords internal
cohens_d <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 2L || length(b) < 2L) return(NA_real_)
  s <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) /
              (length(a) + length(b) - 2))
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(b) - mean(a)) / s
}

#' Games-Howell post-hoc comparisons
#'
#' The post-hoc test that pairs with Welch's ANOVA: it uses Welch's standard
#' error and pair-specific degrees of freedom, so it does not assume the equal
#' variances or equal group sizes that Tukey's HSD does.
#'
#' @param values Numeric vector.
#' @param groups Grouping vector the same length.
#' @return data.frame of pairwise comparisons, or NULL when fewer than two
#'   groups have at least two observations.
#' @keywords internal
games_howell <- function(values, groups) {
  ok <- is.finite(values) & !is.na(groups)
  values <- values[ok]; groups <- as.character(groups)[ok]
  sp <- split(values, groups)
  sp <- sp[vapply(sp, function(x) length(x) >= 2L, logical(1))]
  k <- length(sp)
  if (k < 2L) return(NULL)
  nm <- names(sp)
  n <- vapply(sp, length, integer(1))
  m <- vapply(sp, mean, numeric(1))
  v <- vapply(sp, stats::var, numeric(1))

  rows <- list()
  for (i in seq_len(k - 1L)) for (j in seq(i + 1L, k)) {
    se2 <- v[i] / n[i] + v[j] / n[j]
    if (!is.finite(se2) || se2 <= 0) next
    # Welch-Satterthwaite degrees of freedom for this pair alone.
    df <- se2^2 / ((v[i] / n[i])^2 / (n[i] - 1) + (v[j] / n[j])^2 / (n[j] - 1))
    t <- abs(m[i] - m[j]) / sqrt(se2)
    # The studentised range distribution, which is what makes this a post-hoc
    # test rather than a series of unadjusted t-tests.
    p <- stats::ptukey(t * sqrt(2), nmeans = k, df = df, lower.tail = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      group_a = nm[i], group_b = nm[j], n_a = n[i], n_b = n[j],
      mean_a = m[i], mean_b = m[j], difference = m[j] - m[i],
      df = round(df, 2), statistic = round(t, 4),
      p_value = min(max(p, 0), 1), test = "Games-Howell",
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Dunn's post-hoc comparisons with a tie correction
#'
#' The rank-based post-hoc test that follows a Kruskal-Wallis result. Compares
#' mean ranks computed over the whole data, which is what makes it the correct
#' follow-up rather than a set of pairwise Wilcoxon tests.
#'
#' @param values Numeric vector.
#' @param groups Grouping vector the same length.
#' @return data.frame of pairwise comparisons, or NULL when there are fewer than
#'   two usable groups.
#' @keywords internal
dunn_test <- function(values, groups) {
  ok <- is.finite(values) & !is.na(groups)
  values <- values[ok]; groups <- as.character(groups)[ok]
  N <- length(values)
  if (N < 3L) return(NULL)
  r <- rank(values)
  sp <- split(r, groups)
  sp <- sp[vapply(sp, length, integer(1)) > 0L]
  k <- length(sp)
  if (k < 2L) return(NULL)
  nm <- names(sp)
  n <- vapply(sp, length, integer(1))
  mr <- vapply(sp, mean, numeric(1))

  # Tie correction: without it, ranks that repeat inflate the variance term and
  # the test becomes conservative.
  tie <- table(values)
  ties <- sum(tie^3 - tie)
  sigma2 <- (N * (N + 1) / 12) - (ties / (12 * (N - 1)))

  rows <- list()
  for (i in seq_len(k - 1L)) for (j in seq(i + 1L, k)) {
    se <- sqrt(sigma2 * (1 / n[i] + 1 / n[j]))
    if (!is.finite(se) || se <= 0) next
    z <- (mr[j] - mr[i]) / se
    rows[[length(rows) + 1L]] <- data.frame(
      group_a = nm[i], group_b = nm[j], n_a = n[i], n_b = n[j],
      mean_rank_a = mr[i], mean_rank_b = mr[j],
      difference = mr[j] - mr[i], df = NA_real_,
      statistic = round(z, 4),
      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      test = "Dunn", stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Parametric group comparison with its assumptions recorded
#'
#' Runs the test an immunology paper usually reports, on arcsine square root
#' transformed percentages, and states beside each result whether the two
#' assumptions behind it held.
#'
#' Two groups give Welch's t-test, which does not assume equal variances, plus
#' Student's t-test for reference. More than two give Welch's ANOVA plus the
#' classical one-way ANOVA. Which one is defensible is decided by the
#' Brown-Forsythe column, not by the analyst after seeing the p-values.
#'
#' @param freq Long frequency table with `sample_id`, `population` and a value
#'   column.
#' @param group_of Named vector mapping sample_id to group.
#' @param reference Reference group, or NULL for the first alphabetically.
#' @param value_col Value column; defaults to the resolved abundance measure.
#' @param min_n Minimum samples per group.
#' @param transform "asin_sqrt" (default) or "none".
#' @return data.frame, or NULL when nothing is testable.
#' @keywords internal
parametric_group_tests <- function(freq, group_of, reference = NULL,
                                   value_col = NULL, min_n = 3L,
                                   transform = "asin_sqrt") {
  meas <- if (!is.null(value_col)) list(col = value_col) else
    abundance_measure(freq, "auto")
  d <- freq[!is.na(freq[[meas$col]]), , drop = FALSE]
  d <- qc_pass_rows(d)
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  groups <- sort(unique(as.character(d$group)))
  if (length(groups) < 2L) return(NULL)
  if (is.null(reference) || !reference %in% groups) reference <- groups[1]

  use_asin <- identical(transform, "asin_sqrt")
  rows <- list()
  for (pop in sort(unique(d$population))) {
    dp <- d[d$population == pop, , drop = FALSE]
    raw <- dp[[meas$col]]
    val <- if (use_asin) asin_sqrt_pct(raw) else raw
    g <- as.character(dp$group)
    n_by <- vapply(groups, function(x) sum(g == x & is.finite(val)), integer(1))
    if (any(n_by < min_n)) next

    # ---- assumptions, recorded whatever they say -------------------------
    bf <- brown_forsythe_test(val, g)
    equal_var <- is.finite(bf$p_value) && bf$p_value >= 0.05
    # Shapiro-Wilk on the within-group residuals rather than the raw values: it
    # is the residuals a t-test assumes are normal, not the pooled data, which
    # are bimodal whenever the groups genuinely differ.
    resid <- unlist(lapply(split(val, g), function(x) x - mean(x)), use.names = FALSE)
    sw_p <- if (length(resid) >= 3L && length(resid) <= 5000L &&
                stats::sd(resid) > 0)
      suppressWarnings(stats::shapiro.test(resid)$p.value) else NA_real_
    normal <- is.finite(sw_p) && sw_p >= 0.05

    if (length(groups) == 2L) {
      a <- val[g == reference]; b <- val[g != reference]
      ar <- raw[g == reference]; br <- raw[g != reference]
      wt <- try(stats::t.test(b, a, var.equal = FALSE), silent = TRUE)
      st <- try(stats::t.test(b, a, var.equal = TRUE), silent = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        population = pop, measure = meas$col,
        transform = if (use_asin) "arcsine sqrt" else "none",
        n_groups = 2L, reference_group = reference,
        comparison_group = setdiff(groups, reference)[1],
        n_reference = length(a), n_comparison = length(b),
        mean_reference_pct = mean(ar), mean_comparison_pct = mean(br),
        primary_test = "Welch's t-test",
        statistic = if (inherits(wt, "try-error")) NA_real_ else unname(wt$statistic),
        df = if (inherits(wt, "try-error")) NA_real_ else unname(wt$parameter),
        p_value = if (inherits(wt, "try-error")) NA_real_ else wt$p.value,
        p_student_t = if (inherits(st, "try-error")) NA_real_ else st$p.value,
        cohens_d = round(cohens_d(a, b), 3),
        shapiro_p = sw_p, residuals_normal = normal,
        brown_forsythe_p = bf$p_value, equal_variance = equal_var,
        assumptions_met = normal && equal_var,
        recommended = if (!normal) "rank test (see group_comparison_stats.csv)"
                      else if (!equal_var) "Welch's t-test"
                      else "Student's or Welch's t-test",
        stringsAsFactors = FALSE)
    } else {
      wa <- try(stats::oneway.test(val ~ factor(g), var.equal = FALSE), silent = TRUE)
      ca <- try(stats::oneway.test(val ~ factor(g), var.equal = TRUE), silent = TRUE)
      # eta squared from the classical decomposition, reported because an ANOVA
      # p-value alone says nothing about how much of the variance is between
      # groups.
      gm <- mean(val)
      ssb <- sum(vapply(split(val, g), function(x) length(x) * (mean(x) - gm)^2, numeric(1)))
      sst <- sum((val - gm)^2)
      rows[[length(rows) + 1L]] <- data.frame(
        population = pop, measure = meas$col,
        transform = if (use_asin) "arcsine sqrt" else "none",
        n_groups = length(groups), reference_group = reference,
        comparison_group = NA_character_,
        n_reference = unname(n_by[reference]),
        n_comparison = sum(n_by) - unname(n_by[reference]),
        mean_reference_pct = mean(raw[g == reference]),
        mean_comparison_pct = mean(raw[g != reference]),
        primary_test = "Welch's ANOVA",
        statistic = if (inherits(wa, "try-error")) NA_real_ else unname(wa$statistic),
        df = if (inherits(wa, "try-error")) NA_real_ else unname(wa$parameter[1]),
        p_value = if (inherits(wa, "try-error")) NA_real_ else wa$p.value,
        p_student_t = if (inherits(ca, "try-error")) NA_real_ else ca$p.value,
        cohens_d = NA_real_,
        shapiro_p = sw_p, residuals_normal = normal,
        brown_forsythe_p = bf$p_value, equal_variance = equal_var,
        assumptions_met = normal && equal_var,
        recommended = if (!normal) "Kruskal-Wallis (see group_comparison_stats.csv)"
                      else if (!equal_var) "Welch's ANOVA with Games-Howell"
                      else "one-way ANOVA with Tukey HSD",
        stringsAsFactors = FALSE)
      rows[[length(rows)]]$eta_squared <- if (sst > 0) round(ssb / sst, 3) else NA_real_
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, lapply(rows, function(r) {
    if (is.null(r$eta_squared)) r$eta_squared <- NA_real_
    r
  }))
  out$p_adjusted_BH <- stats::p.adjust(out$p_value, method = "BH")
  out$significant_raw <- !is.na(out$p_value) & out$p_value < 0.05
  out$significant_BH <- !is.na(out$p_adjusted_BH) & out$p_adjusted_BH < 0.05
  out[order(out$p_value, na.last = TRUE), , drop = FALSE]
}

#' Post-hoc comparisons for every population with more than two groups
#'
#' Runs all three: Games-Howell for the unequal-variance case, Tukey HSD for the
#' equal-variance case, and Dunn for the rank case. All three are written, with
#' the assumption columns from the parametric table saying which one to read, so
#' the choice is documented rather than made silently.
#'
#' @param freq,group_of,value_col,min_n As in [parametric_group_tests()].
#' @param transform "asin_sqrt" or "none".
#' @return data.frame, or NULL.
#' @keywords internal
posthoc_group_tests <- function(freq, group_of, value_col = NULL, min_n = 3L,
                                transform = "asin_sqrt") {
  meas <- if (!is.null(value_col)) list(col = value_col) else
    abundance_measure(freq, "auto")
  d <- freq[!is.na(freq[[meas$col]]), , drop = FALSE]
  d <- qc_pass_rows(d)
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  if (length(unique(d$group)) < 3L) return(NULL)

  use_asin <- identical(transform, "asin_sqrt")
  out <- list()
  for (pop in sort(unique(d$population))) {
    dp <- d[d$population == pop, , drop = FALSE]
    val <- if (use_asin) asin_sqrt_pct(dp[[meas$col]]) else dp[[meas$col]]
    g <- as.character(dp$group)
    if (any(table(g) < min_n)) next

    gh <- games_howell(val, g)
    dn <- dunn_test(val, g)
    tk <- tryCatch({
      fit <- stats::aov(val ~ factor(g))
      th <- stats::TukeyHSD(fit)[[1]]
      pr <- do.call(rbind, strsplit(rownames(th), "-", fixed = TRUE))
      data.frame(group_a = pr[, 2], group_b = pr[, 1],
                 n_a = NA_integer_, n_b = NA_integer_,
                 mean_a = NA_real_, mean_b = NA_real_,
                 difference = unname(th[, "diff"]), df = NA_real_,
                 statistic = NA_real_, p_value = unname(th[, "p adj"]),
                 test = "Tukey HSD", stringsAsFactors = FALSE)
    }, error = function(e) NULL)

    for (tb in list(gh, tk, dn)) {
      if (is.null(tb) || !nrow(tb)) next
      tb$population <- pop
      tb$measure <- meas$col
      out[[length(out) + 1L]] <- tb
    }
  }
  if (!length(out)) return(NULL)
  cols <- Reduce(intersect, lapply(out, names))
  res <- do.call(rbind, lapply(out, function(x) x[, cols, drop = FALSE]))
  # Adjusted within each test family, because the three are alternatives rather
  # than one set of comparisons.
  res$p_adjusted_BH <- stats::ave(res$p_value, res$test,
                                  FUN = function(p) stats::p.adjust(p, method = "BH"))
  res[order(res$population, res$test, res$p_value), , drop = FALSE]
}
