# Rank ANCOVA: cohort effect after adjusting for covariates (EXPLORATORY)

WHAT: fits lm(rank(y) ~ covariates + group) per population and reports
the partial F test for `group`. Ranking the response first is the
Conover-Iman rank-transform approach: it keeps the robustness of a rank
test while allowing the design matrix a parametric test cannot otherwise
get from these n.

## Usage

``` r
stats_rank_ancova(
  freq,
  patients,
  smap,
  group_of,
  covariates = c("age", "sex"),
  value_col = NULL,
  min_resid_df = 5L,
  min_n = 3L
)
```

## Arguments

- freq:

  Population frequency table, one row per sample x population.

- patients:

  Patient metadata table, keyed by patient_id.

- smap:

  Sample map joining sample_id to patient_id.

- group_of:

  Named character vector mapping sample_id to group label.

- covariates:

  Column names in the patient table to screen as confounders. Default
  `c("age", "sex")`.

- value_col:

  Column to test or plot instead of the default measure.

- min_resid_df:

  Minimum residual degrees of freedom before a model is fitted. Default
  `5L`.

- min_n:

  Minimum samples per group before a test is attempted. Default `3L`.

## Details

WHY IT IS LABELLED EXPLORATORY IN THE OUTPUT ITSELF, not merely in the
docs: the rank transform is not exact in the presence of covariates, its
Type I error drifts at small n, and – the binding constraint here – if a
covariate is near-collinear with cohort the adjusted estimate is
extrapolation. The `min_resid_df` guard refuses to fit rather than
return a number that looks like the others in the results folder but is
not comparable to them.

READ IT ALONGSIDE stats_confounding(): if that table shows no covariate
is both unbalanced and associated, this adjustment should barely move
the unadjusted p-value, and agreement between the two is the
reassurance. A large disagreement is a finding about the DESIGN, not
about the biology.
