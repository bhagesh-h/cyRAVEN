# =============================================================================
# cyRAVEN -- compositional testing and study-design statistics
#   * compositional (CLR) frequency testing
#   * confounding diagnostics, covariate adjustment, paired designs
#
# These run alongside the primary abundance test rather than replacing it. Each
# produces its own table, so the primary result and its compositional or
# covariate-aware counterpart can be read against each other instead of one
# silently standing in for the other.
# =============================================================================


# =============================================================================
# COMPOSITIONAL FREQUENCIES: testing percentages forced to sum to 100
# =============================================================================
#
# THE GAP THIS CLOSES. abundance_measure() in the baseline states the problem
# exactly right -- "frequencies are compositional: a fall in one population may
# reflect expansion of another, not its own loss" -- and then nothing acts on it.
# Every population's percentage is tested as if it could move independently, and
# it cannot: the D populations of a sample lie on a simplex, so a granulocyte
# expansion MECHANICALLY depresses every lymphocyte percentage, and those
# depressions test significant with no lymphocyte having changed.
#
# THE STANDARD REMEDY (Aitchison 1986; what diffcyt sidesteps by modelling counts
# with a library-size offset instead) is to leave the simplex before testing. The
# centred log-ratio expresses each population's log abundance relative to the
# geometric mean of all populations in that same sample:
#
#     clr(x)_i = log(x_i) - (1/D) * sum_j log(x_j)
#
# After it, a change in population i is a change in i RELATIVE TO that sample's
# overall composition, which is the quantity the design can actually identify.
#
# WHAT IT DOES NOT FIX, stated plainly because a transform that looks like a
# solution invites over-reading: CLR does not recover absolute cell numbers. If
# every population doubles, the composition is unchanged and CLR sees nothing.
# Only cells/uL separates "expanded" from "expanded relative to the rest", which
# is exactly why abundance_measure() already prefers cells_per_ul when the
# patient table supplies wbc_per_ul. CLR is the right TEST for frequency data;
# absolute counts remain the better DATA.
#
# WHY BOTH TESTS ARE REPORTED RATHER THAN THE CLR ONE REPLACING THE RAW ONE: the
# comparison between them is itself the finding. A population significant on raw
# percentage but not on CLR is a candidate artefact of the constraint; one
# significant on both is robust to it. Swapping the test silently would hide that
# distinction, so the CLR results are written alongside and the concordance is
# tabulated per population.

#' Zero-safe centred log-ratio of each sample's population composition
#'
#' WHY ZEROS NEED HANDLING: log(0) is -Inf, and a genuinely absent population
#' (0 cells scored) is common for rare subsets in a small sample. Dropping those
#' rows would drop exactly the populations most likely to differ between cohorts.
#'
#' WHY MULTIPLICATIVE REPLACEMENT: the standard treatment for rounded zeros in
#' compositional data (Martin-Fernandez et al. 2003). Each zero takes a small
#' delta and the non-zero parts are scaled down to preserve the total, so the
#' result is still a composition. Delta defaults to 65% of the smallest non-zero
#' part observed, the usual choice when no detection limit is declared.
#'
#' WHY THE GEOMETRIC MEAN IS PER SAMPLE: the constraint is per sample -- one
#' sample's populations sum to 100, not one population's samples. Centring across
#' samples instead would be a different, and wrong, operation.
#'
#' @param freq population_frequencies table
#' @param value_col the compositional column (percentages)
#' @param populations optional restriction to a closed set; defaults to every
#'   population present, the right closure when they partition CD45+
#' @param delta_frac Fraction of the smallest non-zero part used to replace zeros. Default `0.65`.
#' @return `freq` plus a `clr` column, or NULL
#' @export
clr_frequencies <- function(freq, value_col = "pct_of_cd45_pos",
                            populations = NULL, delta_frac = 0.65) {
  if (is.null(freq) || !nrow(freq) || !value_col %in% names(freq)) return(NULL)
  d <- freq[is.finite(freq[[value_col]]), , drop = FALSE]
  if (!is.null(populations)) d <- d[d$population %in% populations, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  pos <- d[[value_col]][d[[value_col]] > 0]
  if (!length(pos)) return(NULL)
  delta <- delta_frac * min(pos)
  n_zero <- sum(d[[value_col]] == 0, na.rm = TRUE)

  parts <- split(seq_len(nrow(d)), d$sample_id)
  d$clr <- NA_real_
  for (idx in parts) {
    x <- d[[value_col]][idx]
    if (length(x) < 2L) next            # a composition of one part says nothing
    nz <- x > 0
    if (!any(nz)) next
    if (any(!nz)) {
      tot <- sum(x)
      x[!nz] <- delta
      x[nz]  <- x[nz] * (1 - sum(!nz) * delta / max(tot, .Machine$double.eps))
      x[x <= 0] <- delta
    }
    lg <- log(x)
    d$clr[idx] <- lg - mean(lg)
  }
  d <- d[is.finite(d$clr), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  attr(d, "clr_delta") <- delta
  attr(d, "clr_n_zeros") <- n_zero
  d
}

#' Concordance between the raw-percentage test and the CLR test
#'
#' WHAT THE READER GETS: for every population x comparison, whether the
#' conclusion survives removing the compositional constraint, with a one-word
#' verdict naming which of the four cases it falls into.
#'
#' WHY A TABLE AND NOT A LOG LINE: "significant on raw but not on CLR" is a
#' per-population fact and populations behave differently. A single summary
#' sentence would let a reader carry a global impression onto the one population
#' where it does not hold.
#' @param raw_stats The raw stats.
#' @param clr_stats The clr stats.
#' @param alpha Significance threshold. Default `0.05`.
#' @param p_col Which p-value column to compare on. Default `c("p_value", "p_adj_BH")`.
#' @export
compositional_concordance <- function(raw_stats, clr_stats, alpha = 0.05,
                                      p_col = c("p_value", "p_adj_BH")) {
  p_col <- match.arg(p_col)
  if (is.null(raw_stats) || is.null(clr_stats)) return(NULL)
  if (!all(c("population", "comparison_group") %in% names(raw_stats))) return(NULL)
  if (!all(c("population", "comparison_group") %in% names(clr_stats))) return(NULL)
  kf <- function(x) paste(x$population, x$comparison_group, sep = "\r")
  a <- data.frame(k = kf(raw_stats), population = raw_stats$population,
                  comparison_group = raw_stats$comparison_group,
                  p_raw = raw_stats[[p_col]],
                  median_reference = raw_stats$median_reference,
                  median_comparison = raw_stats$median_comparison,
                  cliffs_delta_raw = raw_stats$cliffs_delta,
                  stringsAsFactors = FALSE)
  b <- data.frame(k = kf(clr_stats), p_clr = clr_stats[[p_col]],
                  cliffs_delta_clr = clr_stats$cliffs_delta,
                  stringsAsFactors = FALSE)
  m <- merge(a, b, by = "k", all = FALSE)
  if (!nrow(m)) return(NULL)
  m$k <- NULL
  sig_raw <- !is.na(m$p_raw) & m$p_raw < alpha
  sig_clr <- !is.na(m$p_clr) & m$p_clr < alpha
  m$verdict <- ifelse(sig_raw & sig_clr,  "robust_to_composition",
               ifelse(sig_raw & !sig_clr, "raw_only__possible_composition_artefact",
               ifelse(!sig_raw & sig_clr, "clr_only__was_masked_by_composition",
                                          "not_significant_either")))
  m$p_column <- p_col
  m$alpha <- alpha
  m[order(m$verdict, m$p_clr), , drop = FALSE]
}


# =============================================================================
# CONFOUNDING AND DESIGN: age, sex, and repeated measures
# =============================================================================
#
# THE GAP THIS CLOSES. The patient table carries age and sex; the baseline uses
# them for colouring and never for inference. A cohort difference that is really
# an age difference is reported as a cohort difference. The usual remedy is to
# route the covariates through a diffcyt design matrix (~ group + age + sex) with
# a random effect on sample.
#
# WHY THE PRIMARY DELIVERABLE IS A DIAGNOSTIC, NOT AN ADJUSTMENT. At the group
# sizes this pipeline is built for -- single digits per cohort -- covariate
# ADJUSTMENT is close to unusable. A model with group + age + sex spends most of
# its residual degrees of freedom on nuisance terms, and if age is strongly
# associated with cohort (which in a syndrome-vs-adult-control design it very
# often is) the two are not separable at ANY n: the adjusted estimate is
# extrapolation, and it arrives wearing a confidence interval.
#
# So the order is deliberate:
#   1. ALWAYS -- report whether each covariate differs between cohorts, and
#      whether it tracks the outcome. Those two facts TOGETHER are what makes
#      something a confounder, and both are estimable at small n.
#   2. ON REQUEST, and only where n permits -- a rank-based ANCOVA, labelled
#      exploratory in the output itself.
# A confounder you can see is worth more than an adjustment you cannot trust.

#' Is a covariate confounded with cohort, and does it track the outcome?
#'
#' WHY SPEARMAN / KRUSKAL / FISHER: the same reasoning as everywhere else here --
#' small n, no free normality assumption. Numeric covariates get Spearman
#' (monotone, rank-based). Categorical covariates get Kruskal-Wallis against the
#' outcome and a Fisher exact test against cohort, which stays exact at the cell
#' counts a 2x3 sex-by-cohort table produces where chi-squared does not.
#'
#' @param freq population_frequencies table
#' @param patients patient table (patient_id + the covariates)
#' @param smap sample map (sample_id -> patient_id)
#' @param group_of named vector sample_id -> cohort
#' @param covariates column names in `patients` to examine
#' @param value_col outcome column in `freq`; defaults to the abundance measure
#' @param min_n Minimum samples per group before a test is attempted. Default `3L`.
#' @return data.frame; rows with population = NA carry the covariate-vs-cohort
#'   balance test, the rest carry covariate-vs-outcome association
#' @export
stats_confounding <- function(freq, patients, smap, group_of,
                              covariates = c("age", "sex"),
                              value_col = NULL, min_n = 3L) {
  if (is.null(freq) || !nrow(freq) || is.null(patients) || is.null(smap)) return(NULL)
  if (is.null(value_col)) value_col <- abundance_measure(freq)$col
  if (!value_col %in% names(freq)) return(NULL)
  covariates <- intersect(covariates, names(patients))
  if (!length(covariates)) return(NULL)
  if (!all(c("sample_id", "patient_id") %in% names(smap))) return(NULL)
  if (!"patient_id" %in% names(patients)) return(NULL)

  # sample -> patient -> covariate, joined with norm_id() on both sides because
  # the sample map is hand-typed and the patient table is clinically exported --
  # the same normalisation the baseline uses everywhere these two tables meet.
  pid  <- setNames(norm_id(smap$patient_id), smap$sample_id)
  prow <- match(pid, norm_id(patients$patient_id))
  cov_of <- setNames(lapply(covariates, function(cv)
    setNames(patients[[cv]][prow], names(pid))), covariates)

  d <- qc_pass_rows(freq)
  d <- d[is.finite(d[[value_col]]), , drop = FALSE]
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  is_numeric_cov <- function(v) {
    v <- v[!is.na(v)]
    if (!length(v)) return(FALSE)
    if (is.numeric(v)) return(TRUE)
    suppressWarnings(all(!is.na(as.numeric(as.character(v)))))
  }

  rows <- list()
  for (cv in covariates) {
    v_all <- cov_of[[cv]]
    numeric_cov <- is_numeric_cov(v_all)

    # ---- (a) covariate vs cohort: unbalanced across groups? ------------------
    one <- d[!duplicated(d$sample_id), c("sample_id", "group"), drop = FALSE]
    one$cov <- v_all[one$sample_id]
    one <- one[!is.na(one$cov), , drop = FALSE]
    p_bal <- NA_real_; bal_test <- NA_character_; summ <- "n/a"
    if (nrow(one) >= 2L * min_n && length(unique(one$group)) > 1L) {
      if (numeric_cov) {
        one$cov_n <- suppressWarnings(as.numeric(as.character(one$cov)))
        one <- one[is.finite(one$cov_n), , drop = FALSE]
        if (nrow(one) >= 2L * min_n) {
          kw <- try(stats::kruskal.test(one$cov_n, factor(one$group)), silent = TRUE)
          if (!inherits(kw, "try-error")) {
            p_bal <- kw$p.value; bal_test <- "Kruskal-Wallis"
          }
          sp <- split(one$cov_n, one$group)
          summ <- paste(sprintf("%s=%.1f", names(sp),
                                vapply(sp, stats::median, numeric(1))),
                        collapse = " / ")
        }
      } else {
        tb <- table(as.character(one$cov), one$group)
        if (nrow(tb) > 1L && ncol(tb) > 1L) {
          ft <- try(stats::fisher.test(tb, simulate.p.value = TRUE, B = 20000),
                    silent = TRUE)
          if (!inherits(ft, "try-error")) {
            p_bal <- ft$p.value; bal_test <- "Fisher exact (simulated)"
          }
          summ <- paste(colnames(tb), apply(tb, 2, paste, collapse = ":"),
                        sep = "=", collapse = " / ")
        }
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      covariate = cv, population = NA_character_,
      comparison = "covariate vs cohort", test = bal_test,
      rho = NA_real_, p_value = p_bal, n = nrow(one),
      note = paste0("by cohort: ", summ), stringsAsFactors = FALSE)

    # ---- (b) covariate vs outcome, per population ---------------------------
    for (pop in sort(unique(d$population))) {
      dp <- d[d$population == pop, , drop = FALSE]
      dp$cov <- v_all[dp$sample_id]
      dp <- dp[!is.na(dp$cov), , drop = FALSE]
      if (nrow(dp) < 2L * min_n) next
      if (numeric_cov) {
        cn <- suppressWarnings(as.numeric(as.character(dp$cov)))
        keep <- is.finite(cn) & is.finite(dp[[value_col]])
        if (sum(keep) < 2L * min_n) next
        ct <- try(stats::cor.test(cn[keep], dp[[value_col]][keep],
                                  method = "spearman", exact = FALSE), silent = TRUE)
        if (inherits(ct, "try-error")) next
        rows[[length(rows) + 1L]] <- data.frame(
          covariate = cv, population = pop, comparison = "covariate vs outcome",
          test = "Spearman", rho = unname(round(ct$estimate, 3)),
          p_value = ct$p.value, n = sum(keep),
          note = paste0("outcome = ", value_col), stringsAsFactors = FALSE)
      } else {
        f <- factor(as.character(dp$cov))
        if (nlevels(f) < 2L || min(table(f)) < min_n) next
        kw <- try(stats::kruskal.test(dp[[value_col]], f), silent = TRUE)
        if (inherits(kw, "try-error")) next
        rows[[length(rows) + 1L]] <- data.frame(
          covariate = cv, population = pop, comparison = "covariate vs outcome",
          test = "Kruskal-Wallis", rho = NA_real_, p_value = kw$p.value,
          n = nrow(dp), note = paste0("outcome = ", value_col),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj_BH <- stats::p.adjust(out$p_value, method = "BH")
  # The verdict a reader actually wants. A covariate is a CONFOUNDER for a
  # population only when it is BOTH unbalanced across cohorts AND associated with
  # that population's abundance. Either fact alone is harmless, and reporting
  # either alone as a warning would cry wolf on every study with unequal ages.
  unbal_cov <- unique(out$covariate[out$comparison == "covariate vs cohort" &
                                    !is.na(out$p_value) & out$p_value < 0.05])
  out$confounder_risk <- ifelse(
    out$comparison == "covariate vs outcome" & out$covariate %in% unbal_cov &
      !is.na(out$p_value) & out$p_value < 0.05,
    "HIGH - unbalanced across cohorts AND associated with this population",
    ifelse(out$comparison == "covariate vs outcome" & out$covariate %in% unbal_cov,
           "covariate unbalanced across cohorts, but not associated here", ""))
  out[order(out$comparison, out$covariate, out$p_value), , drop = FALSE]
}

#' Rank ANCOVA: cohort effect after adjusting for covariates (EXPLORATORY)
#'
#' WHAT: fits lm(rank(y) ~ covariates + group) per population and reports the
#' partial F test for `group`. Ranking the response first is the Conover-Iman
#' rank-transform approach: it keeps the robustness of a rank test while allowing
#' the design matrix a parametric test cannot otherwise get from these n.
#'
#' WHY IT IS LABELLED EXPLORATORY IN THE OUTPUT ITSELF, not merely in the docs:
#' the rank transform is not exact in the presence of covariates, its Type I
#' error drifts at small n, and -- the binding constraint here -- if a covariate is
#' near-collinear with cohort the adjusted estimate is extrapolation. The
#' `min_resid_df` guard refuses to fit rather than return a number that looks
#' like the others in the results folder but is not comparable to them.
#'
#' READ IT ALONGSIDE stats_confounding(): if that table shows no covariate is
#' both unbalanced and associated, this adjustment should barely move the
#' unadjusted p-value, and agreement between the two is the reassurance. A large
#' disagreement is a finding about the DESIGN, not about the biology.
#' @param freq Population frequency table, one row per sample x population.
#' @param patients Patient metadata table, keyed by patient_id.
#' @param smap Sample map joining sample_id to patient_id.
#' @param group_of Named character vector mapping sample_id to group label.
#' @param covariates Column names in the patient table to screen as confounders. Default `c("age", "sex")`.
#' @param value_col Column to test or plot instead of the default measure.
#' @param min_resid_df Minimum residual degrees of freedom before a model is fitted. Default `5L`.
#' @param min_n Minimum samples per group before a test is attempted. Default `3L`.
#' @export
stats_rank_ancova <- function(freq, patients, smap, group_of,
                              covariates = c("age", "sex"),
                              value_col = NULL, min_resid_df = 5L, min_n = 3L) {
  if (is.null(freq) || !nrow(freq) || is.null(patients) || is.null(smap)) return(NULL)
  if (is.null(value_col)) value_col <- abundance_measure(freq)$col
  if (!value_col %in% names(freq)) return(NULL)
  covariates <- intersect(covariates, names(patients))
  if (!length(covariates)) return(NULL)
  if (!all(c("sample_id", "patient_id") %in% names(smap))) return(NULL)

  pid  <- setNames(norm_id(smap$patient_id), smap$sample_id)
  prow <- match(pid, norm_id(patients$patient_id))

  d <- qc_pass_rows(freq)
  d <- d[is.finite(d[[value_col]]), , drop = FALSE]
  d$group <- unname(group_of[d$sample_id])
  d <- d[!is.na(d$group), , drop = FALSE]
  if (!nrow(d) || length(unique(d$group)) < 2L) return(NULL)
  for (cv in covariates) {
    v <- setNames(patients[[cv]][prow], names(pid))
    d[[paste0(".cov_", cv)]] <- v[d$sample_id]
  }

  rows <- list()
  for (pop in sort(unique(d$population))) {
    dp <- d[d$population == pop, , drop = FALSE]
    cvc <- paste0(".cov_", covariates)
    dp <- dp[stats::complete.cases(dp[, c(value_col, "group", cvc), drop = FALSE]), ,
             drop = FALSE]
    if (nrow(dp) < 2L * min_n) next
    dp$.y <- rank(dp[[value_col]])
    dp$group <- factor(dp$group)
    if (nlevels(dp$group) < 2L) next
    # Coerce each covariate to the narrowest type it genuinely is: a numeric-like
    # column read as character would otherwise become a factor with one level per
    # donor and consume every degree of freedom.
    use <- character(0)
    for (cv in covariates) {
      k <- paste0(".cov_", cv)
      num <- suppressWarnings(as.numeric(as.character(dp[[k]])))
      if (all(is.finite(num))) dp[[k]] <- num
      else {
        dp[[k]] <- factor(as.character(dp[[k]]))
        if (nlevels(dp[[k]]) < 2L) next
      }
      use <- c(use, k)
    }
    if (!length(use)) next

    full <- try(stats::lm(stats::reformulate(c(use, "group"), ".y"), data = dp),
                silent = TRUE)
    red  <- try(stats::lm(stats::reformulate(use, ".y"), data = dp), silent = TRUE)
    if (inherits(full, "try-error") || inherits(red, "try-error")) next
    rdf <- stats::df.residual(full)
    if (!is.finite(rdf) || rdf < min_resid_df) {
      rows[[length(rows) + 1L]] <- data.frame(
        population = pop, measure = value_col,
        covariates = paste(covariates, collapse = "+"),
        n = nrow(dp), residual_df = rdf, F_group = NA_real_, p_group_adjusted = NA_real_,
        status = paste0("NOT FITTED - only ", rdf, " residual df (need ",
                        min_resid_df, "); too few samples to adjust"),
        stringsAsFactors = FALSE)
      next
    }
    an <- try(stats::anova(red, full), silent = TRUE)
    if (inherits(an, "try-error") || nrow(an) < 2L) next
    rows[[length(rows) + 1L]] <- data.frame(
      population = pop, measure = value_col,
      covariates = paste(covariates, collapse = "+"),
      n = nrow(dp), residual_df = rdf,
      F_group = round(an$F[2], 3), p_group_adjusted = an[["Pr(>F)"]][2],
      status = "EXPLORATORY - rank ANCOVA; compare with the unadjusted Wilcoxon",
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj_BH <- stats::p.adjust(out$p_group_adjusted, method = "BH")
  out[order(out$p_group_adjusted), , drop = FALSE]
}

#' Paired / repeated-measures abundance test
#'
#' WHAT: when samples pair up (pre/post, longitudinal, matched donors), compares
#' WITHIN pair instead of between groups.
#'
#' WHY IT CANNOT BE INFERRED AND MUST BE DECLARED: nothing in an FCS file or a
#' filename says two tubes came from the same donor at two timepoints. Treating
#' paired samples as independent -- which is what the baseline necessarily does,
#' having no notion of pairing -- discards the pairing's variance reduction AND
#' violates the independence assumption the rank test rests on. So this runs only
#' when --paired-column names a column, and it counts incomplete pairs rather
#' than silently dropping them, because a pairing that quietly halves n is
#' exactly the failure this function exists to prevent.
#'
#' WHY WILCOXON SIGNED-RANK / FRIEDMAN: the paired counterparts of the rank tests
#' used everywhere else here. Friedman generalises to 3+ conditions and is the
#' conventional choice for repeated-measures frequency data at this scale.
#'
#' @param pair_of named vector sample_id -> pairing unit (donor)
#' @param cond_of named vector sample_id -> condition (timepoint / arm)
#' @param freq Population frequency table, one row per sample x population.
#' @param reference The group every other group is compared against.
#' @param value_col Column to test or plot instead of the default measure.
#' @param min_pairs Minimum complete pairs before a paired test is attempted. Default `3L`.
#' @export
stats_paired_comparison <- function(freq, pair_of, cond_of, reference = NULL,
                                    value_col = NULL, min_pairs = 3L) {
  if (is.null(freq) || !nrow(freq)) return(NULL)
  if (is.null(value_col)) value_col <- abundance_measure(freq)$col
  if (!value_col %in% names(freq)) return(NULL)
  d <- qc_pass_rows(freq)
  d <- d[is.finite(d[[value_col]]), , drop = FALSE]
  d$pair <- unname(pair_of[d$sample_id])
  d$cond <- unname(cond_of[d$sample_id])
  d <- d[!is.na(d$pair) & !is.na(d$cond), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  conds <- sort(unique(as.character(d$cond)))
  if (length(conds) < 2L) return(NULL)
  if (is.null(reference) || !reference %in% conds) reference <- conds[1]

  rows <- list()
  for (pop in sort(unique(d$population))) {
    dp <- d[d$population == pop, c("pair", "cond", value_col), drop = FALSE]
    # Average replicate tubes of the same (pair, condition) before reshaping: two
    # rows for one donor-timepoint is a duplicate measurement, not a second pair,
    # and reshape() would otherwise silently keep only one of them.
    dp <- stats::aggregate(stats::reformulate(c("pair", "cond"), value_col),
                           dp, mean)
    w <- stats::reshape(dp, idvar = "pair", timevar = "cond", direction = "wide")
    vn <- paste0(value_col, ".", conds)
    have <- intersect(vn, names(w))
    if (length(have) < 2L) next
    complete <- stats::complete.cases(w[, have, drop = FALSE])
    n_drop <- sum(!complete)
    w <- w[complete, , drop = FALSE]
    if (nrow(w) < min_pairs) next

    p_omni <- NA_real_
    if (length(conds) > 2L && length(have) == length(conds)) {
      fr <- try(stats::friedman.test(as.matrix(w[, have, drop = FALSE])), silent = TRUE)
      if (!inherits(fr, "try-error")) p_omni <- fr$p.value
    }
    rcol <- paste0(value_col, ".", reference)
    if (!rcol %in% names(w)) next
    for (cd in setdiff(conds, reference)) {
      ccol <- paste0(value_col, ".", cd)
      if (!ccol %in% names(w)) next
      a <- w[[rcol]]; b <- w[[ccol]]
      wt <- try(stats::wilcox.test(a, b, paired = TRUE, exact = FALSE), silent = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        population = pop, measure = value_col,
        reference_condition = reference, comparison_condition = cd,
        n_pairs = nrow(w), n_incomplete_pairs_dropped = n_drop,
        median_reference = stats::median(a), median_comparison = stats::median(b),
        median_paired_difference = stats::median(b - a),
        # Matched-pairs rank-biserial: the share of pairs moving up minus the
        # share moving down. Bounded \code{[-1,1]} and reads directly as "how consistent
        # is the DIRECTION of change", which is the paired question.
        rank_biserial = round((sum(b > a) - sum(b < a)) / max(1L, sum(b != a)), 3),
        test = "Wilcoxon signed-rank (paired)",
        p_value = if (inherits(wt, "try-error")) NA_real_ else wt$p.value,
        p_omnibus_friedman = p_omni, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj_BH <- stats::p.adjust(out$p_value, method = "BH")
  out$significant_raw <- !is.na(out$p_value) & out$p_value < 0.05
  out$significant_BH  <- !is.na(out$p_adj_BH) & out$p_adj_BH < 0.05
  out[order(out$p_value), , drop = FALSE]
}
