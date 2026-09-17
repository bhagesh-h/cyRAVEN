# Per-patient trajectories across timepoints, split by outcome

WHY LINES AND NOT A BOX PER TIMEPOINT. In a repeated-measures cohort the
signal is usually the DIRECTION each patient moves, and a box per
timepoint averages exactly that away: two patients rising and two
falling look like a flat cohort. In the sepsis literature this is the
figure that separates survivors from non-survivors, because they differ
in trajectory rather than at any one timepoint.

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
  and to split the mean.

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

  palette.

## Details

WHAT IT IS NOT. Descriptive. The test for a repeated-measures design is
the paired comparison in `paired_comparison_stats.csv`; a trajectory
figure shows the shape of the data that test was run on.
