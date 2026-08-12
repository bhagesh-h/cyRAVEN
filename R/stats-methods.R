# The statistical methods catalogue, and the evidence for the ones chosen.
#
# WHY THIS EXISTS
#
# Immunophenotyping papers report a small, well-known set of tests: Student's or
# Welch's t-test, one-way and two-way ANOVA with Tukey, Mann-Whitney, Kruskal-
# Wallis with Dunn, and in the cytometry-specific literature the diffcyt family
# (edgeR, limma-voom, GLMM). A reader arriving at a cyRAVEN results folder
# reasonably asks why they are seeing Wilcoxon and Kruskal-Wallis and not the
# others.
#
# `statistical_methods.csv` answers that in one table: every commonly reported
# method, whether this run computed it, and where it did not, the reason.
#
# WHAT IS DELIBERATELY NOT DONE
#
# The obvious alternative is to run all of them and report every p-value. That
# is worse than useless. With four to ten donors per group, running a t-test and
# a Mann-Whitney and an ANOVA on the same data and reporting whichever is
# smallest is p-hacking with extra steps, and a table of six p-values per
# population invites exactly that. cyRAVEN computes one test per question and
# says why it is that one.
#
# WHAT IS ADDED INSTEAD
#
# The evidence for the choice. `normality_tests.csv` carries a Shapiro-Wilk test
# per population per group and a Brown-Forsythe test for equal variances across
# groups. Those are the two assumptions a t-test or ANOVA rests on, so the
# reader can check the rank test was the right call rather than taking it on
# trust -- and can see the cases where it made no difference.
#
# Note what the normality table usually shows at these group sizes: Shapiro-Wilk
# on 4 to 10 values has almost no power, so a non-significant result is NOT
# evidence of normality. That is itself the argument for the rank test, and the
# table says so in its own `interpretation` column rather than leaving a reader
# to conclude "p > 0.05, so normal, so a t-test would have been fine".

#' Brown-Forsythe test for equality of variances
#'
#' Levene's test computed on absolute deviations from the group MEDIAN rather
#' than the mean, which is the variant that holds up when the data are skewed --
#' the usual case for population frequencies bounded at zero.
#'
#' @param x Numeric values.
#' @param g Grouping factor.
#' @return list(statistic, p_value, df1, df2), or NULLs when undefined.
#' @keywords internal
brown_forsythe_test <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- factor(g[ok])
  k <- nlevels(g); n <- length(x)
  if (k < 2L || n - k < 1L) return(list(statistic = NA_real_, p_value = NA_real_,
                                        df1 = NA_integer_, df2 = NA_integer_))
  z <- abs(x - stats::ave(x, g, FUN = function(v) stats::median(v, na.rm = TRUE)))
  zbar <- mean(z)
  zg <- tapply(z, g, mean)
  ng <- tapply(z, g, length)
  num <- sum(ng * (zg - zbar)^2) / (k - 1)
  den <- sum((z - stats::ave(z, g, FUN = mean))^2) / (n - k)
  if (!is.finite(den) || den <= 0) return(list(statistic = NA_real_,
                                               p_value = NA_real_,
                                               df1 = k - 1L, df2 = n - k))
  f <- num / den
  list(statistic = f, p_value = stats::pf(f, k - 1, n - k, lower.tail = FALSE),
       df1 = k - 1L, df2 = n - k)
}

#' Normality and equal-variance evidence, per population and group
#'
#' The diagnostics behind cyRAVEN's use of rank tests. One Shapiro-Wilk per
#' population per group, plus one Brown-Forsythe per population across groups.
#'
#' @param freq Long data frame with `sample_id`, `population` and a value column.
#' @param group_of Named character vector, sample id -> group.
#' @param value_col Column holding the per-sample value.
#' @return data.frame, or NULL when there is nothing testable.
#' @keywords internal
normality_report <- function(freq, group_of, value_col = "pct_of_cd45_pos") {
  if (is.null(freq) || !nrow(freq) || !value_col %in% names(freq)) return(NULL)
  d <- freq[is.finite(freq[[value_col]]), , drop = FALSE]
  if ("qc_status" %in% names(d)) d <- d[d$qc_status %in% c("pass", TRUE, "TRUE"), , drop = FALSE]
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  rows <- list()
  for (pop in sort(unique(d$population))) {
    dp <- d[d$population == pop, , drop = FALSE]
    for (g in sort(unique(dp$group))) {
      v <- dp[[value_col]][dp$group == g]
      n <- length(v)
      sw_w <- NA_real_; sw_p <- NA_real_
      # shapiro.test needs 3..5000 values and a non-zero variance.
      if (n >= 3L && n <= 5000L && stats::sd(v) > 0) {
        st <- try(stats::shapiro.test(v), silent = TRUE)
        if (!inherits(st, "try-error")) {
          sw_w <- unname(st$statistic); sw_p <- st$p.value
        }
      }
      interp <- if (!is.finite(sw_p)) {
        if (n < 3L) "too few donors to test" else "no variance to test"
      } else if (sw_p < 0.05) {
        "departs from normal; a t-test or ANOVA is not supported here"
      } else if (n < 15L) {
        paste0("not rejected, but Shapiro-Wilk on ", n,
               " values has little power - this is NOT evidence of normality")
      } else {
        "consistent with normal"
      }
      rows[[length(rows) + 1L]] <- data.frame(
        population = pop, group = g, n_donors = n,
        median = round(stats::median(v), 6),
        iqr = round(stats::IQR(v), 6),
        shapiro_w = round(sw_w, 4), shapiro_p = signif(sw_p, 4),
        interpretation = interp, stringsAsFactors = FALSE)
    }
    bf <- brown_forsythe_test(dp[[value_col]], dp$group)
    rows[[length(rows) + 1L]] <- data.frame(
      population = pop, group = "(across groups)",
      n_donors = nrow(dp), median = NA_real_, iqr = NA_real_,
      shapiro_w = round(bf$statistic, 4), shapiro_p = signif(bf$p_value, 4),
      interpretation = if (!is.finite(bf$p_value)) "variance test undefined"
        else if (bf$p_value < 0.05)
          "Brown-Forsythe: variances differ between groups; Student's t-test would be invalid"
        else "Brown-Forsythe: no evidence of unequal variance",
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' The methods catalogue
#'
#' Every commonly reported test in the immunophenotyping and cytometry
#' literature, what this run did about it, and why. Written on every run so the
#' reasoning travels with the numbers.
#'
#' @param n_groups Number of group levels compared, or NA.
#' @param paired Whether a paired design was detected.
#' @param n_tests Number of tests corrected together, or NA.
#' @return data.frame.
#' @keywords internal
statistical_methods_table <- function(n_groups = NA_integer_, paired = FALSE,
                                      n_tests = NA_integer_) {
  used <- function(what, applied, why) data.frame(
    method = what, role = "used", applied_to = applied, rationale = why,
    stringsAsFactors = FALSE)
  not <- function(what, why) data.frame(
    method = what, role = "not used", applied_to = NA_character_,
    rationale = why, stringsAsFactors = FALSE)

  two <- is.finite(n_groups) && n_groups == 2L
  many <- is.finite(n_groups) && n_groups > 2L

  rows <- list(
    used("Wilcoxon rank-sum (Mann-Whitney U)",
         "each group against the reference, per population",
         paste("Rank based, so it assumes neither normality nor equal variance.",
               "At 4 to 10 donors per group normality cannot be verified, which",
               "is what rules out the parametric alternatives.")),
    if (many) used("Kruskal-Wallis",
         "omnibus across all groups, per population",
         "Asks whether any group differs before the pairwise tests are read.")
    else not("Kruskal-Wallis",
         "Needs three or more groups; this run compares fewer."),
    if (isTRUE(paired)) used("Wilcoxon signed-rank",
         "paired samples, per population",
         "The paired counterpart of the rank-sum test.")
    else not("Wilcoxon signed-rank / paired t-test",
         "No paired design was declared. Requires repeated samples per donor."),
    used("Cliff's delta",
         "reported beside every p-value",
         paste("Effect size. A p-value states how surprising a difference is,",
               "never how large. Delta is the probability a random donor from",
               "one group exceeds one from the other, minus the reverse.")),
    used("Benjamini-Hochberg",
         if (is.finite(n_tests)) paste(n_tests, "tests corrected together")
         else "every test in the family",
         paste("Controls the false discovery rate. Twenty tests at p < 0.05",
               "produce about one false positive by chance alone, so the raw",
               "p-value is not the number to read.")),
    used("Shapiro-Wilk",
         "per population per group, in normality_tests.csv",
         paste("The evidence for choosing a rank test. Reported so the choice",
               "can be checked rather than taken on trust.")),
    used("Brown-Forsythe",
         "per population across groups, in normality_tests.csv",
         paste("Equality of variances, the second assumption behind a t-test",
               "or ANOVA. Median-centred, which survives skew.")),
    used("Wilson score interval",
         "counting uncertainty on every frequency",
         paste("Confidence interval for a proportion that stays sensible near",
               "0% and 100% and at small counts, where the textbook normal",
               "approximation does not.")),
    used("Cramer's V",
         "batch against group, in batch_group_confounding.csv",
         paste("Association between two categorical variables. Above a set",
               "value, batch correction is refused rather than attempted.")),
    used("Permutation null (iLISI)",
         "batch mixing on the embedding",
         "Gives the score expected under no batch structure, by shuffling labels."),
    used("Centred log-ratio",
         "compositional re-test of every abundance finding",
         paste("Frequencies of a whole sum to 100%, so one population rising",
               "forces others down. The CLR removes that constraint.")),

    not("Student's t-test",
        paste("Assumes normality and equal variances. Neither can be verified",
              "at these group sizes - see normality_tests.csv - and the rank",
              "test costs little power while assuming neither.")),
    not("Welch's t-test",
        paste("Drops the equal-variance assumption but keeps normality.",
              "Same objection.")),
    not("One-way ANOVA with Tukey",
        if (many) paste("The parametric counterpart of Kruskal-Wallis with",
                        "Dunn, which this run reports instead, for the reason",
                        "given against Student's t-test.")
        else "Needs three or more groups."),
    not("Two-way / repeated-measures ANOVA",
        paste("For factorial or within-donor designs. cyRAVEN diagnoses",
              "covariates rather than modelling them; see",
              "confounding_diagnostics.csv.")),
    not("diffcyt: edgeR, limma-voom, GLMM",
        paste("Cluster-level cytometry methods that borrow strength across",
              "clusters via an empirical Bayes prior. That prior pays off at",
              "mass-cytometry marker counts and many clusters; on a declared",
              "panel of about a dozen populations there is little to borrow,",
              "and the moderation is harder to explain than it is worth.",
              "Available through --external-labels if wanted.")),
    not("Chi-squared / Fisher's exact",
        paste("For counts of donors in categories, not for per-donor",
              "frequencies. Cramer's V above is the chi-squared based",
              "association actually used.")),
    not("Pearson correlation",
        "Assumes linearity and normality; Spearman is used where rank correlation is needed."),
    not("Bonferroni",
        paste("Controls the family-wise error rate, which is stricter than",
              "needed for a screen across populations. BH is used instead;",
              "multiply the raw p-value by the test count if Bonferroni is",
              "required by a journal."))
  )
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  out$n_groups <- n_groups
  out
}
