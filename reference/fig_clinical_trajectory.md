# Per-patient trajectories across timepoints, split by outcome

Thin lines are patients, the thick line is the median at each timepoint
within each outcome group.

## Usage

``` r
fig_clinical_trajectory(
  freq,
  time_of,
  patient_of,
  outfile,
  outcome_of = NULL,
  outcome_name = "outcome",
  value_col = NULL,
  ncol = 4L,
  time_levels = NULL,
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- freq:

  population_frequencies table.

- time_of:

  named vector sample_id -\> timepoint.

- patient_of:

  named vector sample_id -\> patient.

- outfile:

  path.

- outcome_of:

  optional named vector sample_id -\> outcome, used to colour the lines
  and to split the median.

- outcome_name:

  label for the outcome in the key.

- value_col:

  frequency column.

- ncol:

  panels per row.

- time_levels:

  order of the timepoints; sorted unique values by default.

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.

## Details

Lines rather than a box per timepoint, because in a repeated-measures
cohort the signal is usually the direction each patient moves, and a box
per timepoint averages exactly that away: two patients rising and two
falling look like a flat cohort. In the sepsis literature this is the
figure that separates survivors from non-survivors, which differ in
trajectory rather than at any one timepoint.

Descriptive. The test for this design is the paired comparison in
`paired_comparison_stats.csv`, which uses only patients present at both
timepoints. A line that stops early is a patient who was not sampled
again.

## See also

[`stats_paired_comparison`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_paired_comparison.md)
for the test on the same design.
