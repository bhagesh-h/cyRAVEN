# The cohort as one picture: populations, samples and every clinical variable

One column per sample, ordered along the clinical variable named in
`order_by`. Below the strips, one row per population, filled by that
sample's value expressed as a z-score within the population, so a row is
readable whether the population is 40% of CD45+ or 0.4%. Above them, one
strip per clinical variable.

## Usage

``` r
fig_clinical_landscape(
  freq,
  clin,
  outfile,
  order_by = NULL,
  value_col = NULL,
  max_vars = 6L,
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- freq:

  population_frequencies table.

- clin:

  named list of clinical variables, sample_id -\> value.

- outfile:

  path.

- order_by:

  which clinical variable orders the columns; the first numeric one by
  default.

- value_col:

  frequency column.

- max_vars:

  most clinical strips to draw, so the key cannot outgrow the figure.

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The patchwork composition, invisibly.

## Details

The z-score is taken within a row and not across the matrix. Populations
differ by two orders of magnitude in abundance; on a shared scale the
frequent lineages saturate the palette and every rare one is the same
colour, so the figure would show which populations are large – which is
already known – instead of which samples are unusual.

Nothing here is a test. A gradient that reads clearly across twelve
ordered columns can still be the ordering the eye was given, and the
tests in `clinical_association.csv` are what decide.

## See also

[`stats_clinical_association`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_clinical_association.md)
for the tests,
[`fig_clinical_correlogram`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_correlogram.md)
for whether the strips are independent of each other.
