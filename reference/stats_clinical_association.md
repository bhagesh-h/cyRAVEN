# Clinical variables against population frequencies and marker intensities

One question, asked per population and per marker: does this quantity
move with that clinical variable? Numeric variables are tested with
Spearman's rho, two-level variables with Wilcoxon rank-sum and Cliff's
delta, and variables with more levels with Kruskal-Wallis with
epsilon-squared. P-values are adjusted with Benjamini-Hochberg within
each variable, because each variable is its own question asked of every
population.

## Usage

``` r
stats_clinical_association(
  freq,
  mfi = NULL,
  clin = list(),
  patient_of = NULL,
  value_col = NULL,
  min_n = 4L
)
```

## Arguments

- freq:

  population_frequencies table.

- mfi:

  population_marker_mfi table, optional.

- clin:

  named list of clinical variables, sample_id -\> value.

- patient_of:

  named vector sample_id -\> patient. When given, a variable that is
  constant within a patient is tested on one value per patient rather
  than one per sample.

- value_col:

  the frequency column to use.

- min_n:

  fewest observations for a test.

## Value

list(populations, markers)

## The unit of analysis

Every test runs on one value per sample or per patient, never on pooled
cells. Which of the two depends on the variable, and the choice is
recorded in the `unit` column.

A variable that is constant within a patient – 28-day survival, BMI, age
– asks a between-subject question, and repeated samples from one patient
are repeated measurements of one person. Testing them as independent
observations is pseudoreplication: it inflates the apparent sample size,
narrows the interval and makes significance more likely without one
extra patient having been recruited. Those variables are collapsed to
one value per patient, and `n` is the patient count.

A variable that moves within a patient – a severity score taken at each
draw – asks a within-subject question, and collapsing would delete the
variation being asked about. The per-sample test is kept,
`repeated_measures` is set, and `n_patients` is reported beside `n`. The
dedicated method for that question is a repeated-measures correlation,
discussed in
[`vignette("statistics", package = "cyRAVEN")`](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.md).

With no `patient_of` every sample is treated as its own subject, which
is correct for a cross-sectional cohort and is all that can be assumed
when the sheet does not say otherwise.

## See also

[`fig_clinical_forest`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_forest.md),
[`fig_clinical_correlogram`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_correlogram.md),
[`fig_clinical_landscape`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_landscape.md).
