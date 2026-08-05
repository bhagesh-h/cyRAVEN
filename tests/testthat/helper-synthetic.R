# Synthetic fixtures with PLANTED structure.
#
# Every generator below builds data whose answer is known in advance, so the
# tests can assert that a function finds that specific thing rather than merely
# returning a data frame of the right shape. A test that only checks nrow() > 0
# passes just as happily when the statistic is wrong.

MARKERS <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14")
POPS    <- c("CD4 T cells", "CD8 T cells", "B cells", "NK cells",
             "Classical monocytes")

#' 16 samples in three cohorts, 6 of them the reference.
synth_groups <- function() {
  g <- c(rep("HC", 6), rep("S1", 5), rep("S2", 5))
  stats::setNames(g, sprintf("S%02d", seq_along(g)))
}

#' Per-sample x population x marker medians.
#'
#' PLANTED: CD4 in "CD4 T cells" is +1.5 arcsinh units in cohort S1 only.
#' Nothing else differs, so any function that reports a different top hit is
#' finding noise.
synth_mfi <- function(seed = 7, n_cells = 500L) {
  withr::local_seed(seed)
  g <- synth_groups()
  do.call(rbind, lapply(names(g), function(s) {
    do.call(rbind, lapply(POPS, function(p) {
      eff <- if (p == "CD4 T cells" && g[[s]] == "S1") 1.5 else 0
      data.frame(
        sample_id = s, panel = "P1", population = p, n_cells = n_cells,
        marker = MARKERS,
        median_asinh = stats::rnorm(length(MARKERS), 2, 0.25) +
          c(0, eff, 0, 0, 0, 0),
        pct_positive = pmin(100, pmax(0, stats::rnorm(length(MARKERS), 40, 6) +
                                        c(0, eff * 12, 0, 0, 0, 0))),
        is_control = FALSE, qc_status = "pass", stringsAsFactors = FALSE)
    }))
  }))
}

#' Population frequencies that genuinely sum to 100 per sample.
#'
#' PLANTED: CD4 T cells are truly depressed in both patient cohorts. Because the
#' parts are closed, the OTHER populations will appear to rise — which is exactly
#' the compositional artefact the CLR test exists to separate out.
synth_freq <- function(seed = 11) {
  withr::local_seed(seed)
  g <- synth_groups()
  do.call(rbind, lapply(names(g), function(s) {
    base <- c(30, 20, 12, 10, 28)
    if (g[[s]] != "HC") base <- base * c(0.55, 1, 1, 1, 1)
    v <- base * exp(stats::rnorm(length(base), 0, 0.08))
    v <- 100 * v / sum(v)
    data.frame(sample_id = s, panel = "P1", population = POPS,
               pct_of_cd45_pos = v, count = round(v * 100),
               is_control = FALSE, qc_status = "pass", stringsAsFactors = FALSE)
  }))
}

#' Per-sample thresholds.
#' PLANTED: the CD4 cut sits ~0.9 units higher in both patient cohorts, so the
#' population's DEFINITION moves with cohort.
synth_thresholds <- function(seed = 13) {
  withr::local_seed(seed)
  g <- synth_groups()
  do.call(rbind, lapply(names(g), function(s) data.frame(
    sample_id = s, panel = "P1", marker = MARKERS,
    threshold = stats::rnorm(length(MARKERS), 1.5, 0.02) +
      c(0, if (g[[s]] == "HC") 0 else 0.9, 0, 0, 0, 0),
    source = "valley", is_control = FALSE, qc_status = "pass",
    stringsAsFactors = FALSE)))
}

#' Embedded cells with well-separated populations along UMAP-1.
#' PLANTED: batch is perfectly confounded with cohort (HC on day1, patients on
#' day2), which is the situation where no batch correction can be safe.
synth_cells <- function(seed = 17, per_sample = 120L) {
  withr::local_seed(seed)
  g <- synth_groups()
  labs <- c(POPS, "Other CD45+")
  ctr <- stats::setNames(seq_along(labs) * 4, labs)
  do.call(rbind, lapply(names(g), function(s) {
    pl <- sample(labs, per_sample, TRUE)
    d <- data.frame(
      sample_id = s, population_label = pl, cohort = g[[s]],
      batch = if (g[[s]] == "HC") "day1" else "day2",
      event_index = seq_len(per_sample),
      umap_1 = ctr[pl] + stats::rnorm(per_sample, 0, 0.8),
      umap_2 = stats::rnorm(per_sample, 0, 1.2), stringsAsFactors = FALSE)
    for (m in MARKERS)
      d[[m]] <- stats::rnorm(per_sample,
                             ifelse(pl == "CD4 T cells" & m == "CD4", 4, 1), 0.5)
    d
  }))
}

synth_patients <- function(seed = 19) {
  withr::local_seed(seed)
  ids <- names(synth_groups())
  # PLANTED: age is strongly unbalanced across cohorts (adult controls, paediatric
  # patients) — the classic confounding this design cannot escape.
  data.frame(patient_id = ids,
             age = c(stats::rnorm(6, 40, 5), stats::rnorm(10, 12, 3)),
             sex = sample(c("male", "female"), length(ids), TRUE),
             stringsAsFactors = FALSE)
}

synth_smap <- function() {
  ids <- names(synth_groups())
  data.frame(sample_id = ids, patient_id = ids,
             donor = rep(sprintf("D%d", seq_len(8)), each = 2),
             timepoint = rep(c("pre", "post"), 8), stringsAsFactors = FALSE)
}
