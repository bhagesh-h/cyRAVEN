# Is a covariate confounded with cohort, and does it track the outcome?

WHY SPEARMAN / KRUSKAL / FISHER: the same reasoning as everywhere else
here – small n, no free normality assumption. Numeric covariates get
Spearman (monotone, rank-based). Categorical covariates get
Kruskal-Wallis against the outcome and a Fisher exact test against
cohort, which stays exact at the cell counts a 2x3 sex-by-cohort table
produces where chi-squared does not.

## Usage

``` r
stats_confounding(
  freq,
  patients,
  smap,
  group_of,
  covariates = c("age", "sex"),
  value_col = NULL,
  min_n = 3L
)
```

## Arguments

- freq:

  population_frequencies table

- patients:

  patient table (patient_id + the covariates)

- smap:

  sample map (sample_id -\> patient_id)

- group_of:

  named vector sample_id -\> cohort

- covariates:

  column names in `patients` to examine

- value_col:

  outcome column in `freq`; defaults to the abundance measure

- min_n:

  Minimum samples per group before a test is attempted. Default `3L`.

## Value

data.frame; rows with population = NA carry the covariate-vs-cohort
balance test, the rest carry covariate-vs-outcome association
