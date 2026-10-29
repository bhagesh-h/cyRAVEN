# Cohort-level MFI per population and marker, overall and per study group

One row per panel, population (or cluster), marker and group, where the
group is "All samples" plus each study group when there are two or more.
Controls and staining-QC failures are dropped first, as for every group
comparison.

## Usage

``` r
summarise_mfi(
  mfi,
  transforms = NULL,
  group_of = NULL,
  label_col = "population",
  value_col = "median_asinh",
  pos_col = "pct_positive",
  count_col = "n_cells",
  scatter = character(0)
)
```

## Arguments

- mfi:

  Per-sample table: `sample_id`, the label column, `marker`, the
  analysis-scale median, and optionally `panel`, `pct_positive`,
  `n_cells`.

- transforms:

  Either one transform or a list of them named by panel.

- group_of:

  Optional named vector, sample_id -\> group.

- label_col:

  Name of the population or cluster column.

- value_col:

  Name of the analysis-scale median column.

- pos_col:

  Name of the percent-positive column, if any.

- count_col:

  Name of the cell-count column, if any.

- scatter:

  Markers on a log10 scatter scale. See
  [`mfi_back_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/mfi_back_transform.md).

## Value

data.frame, or NULL when there is nothing to summarise.
