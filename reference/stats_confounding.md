# Is a covariate confounded with cohort, and does it track the outcome?

WHY RANK TESTS AND FISHER: the same reasoning as everywhere else here –
small n, no free normality assumption. A Fisher exact test asks whether
a categorical covariate is balanced across cohorts, which stays exact at
the cell counts a 2x3 sex-by-cohort table produces where chi-squared
does not; a numeric one gets Kruskal-Wallis.

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

## Details

THE OUTCOME ASSOCIATION IS TAKEN WITHIN COHORTS. A confounder is
associated with the exposure and with the outcome conditional on the
exposure (Rothman, Greenland & Lash 2008, Modern Epidemiology, 3rd ed.).
Tested across the pooled cohorts, a real group difference makes any
unbalanced covariate track the outcome by itself: with two cohorts of
ten, ages 30 and 55 and a group effect on the outcome with none of age,
the pooled Spearman test called age a HIGH confounder risk in 199 of 200
simulated cohorts. So the association is a partial Spearman correlation
controlling for cohort – the ranks with each cohort's mean rank removed,
tested on n - 2 - (k - 1) degrees of freedom, as ppcor does – and, for a
categorical covariate, a rank F test of the covariate after the cohort
term. Where age and cohort barely overlap this has little power, which
is the honest answer: the design cannot separate them, and the balance
row already says the covariate is unbalanced.
