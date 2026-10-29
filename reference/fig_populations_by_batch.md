# Population abundance per batch

WHY THIS EXISTS SEPARATELY from the batch diagnostics. Those ask whether
the EMBEDDING separates by batch, which is a question about the cells.
This asks whether the reported numbers differ by batch, which is a
question about the result: a population that shifts with acquisition
batch is a population whose between-group difference may be an
acquisition difference.

## Usage

``` r
fig_populations_by_batch(
  freq,
  batch_of,
  outfile,
  value_col = NULL,
  ncol = 4L,
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- freq:

  population_frequencies table.

- batch_of:

  named vector sample_id -\> batch.

- outfile:

  path.

- value_col:

  frequency column.

- ncol:

  panels per row.

- dpi:

  resolution.

- colors:

  palette.
