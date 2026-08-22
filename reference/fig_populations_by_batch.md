# Population abundance per acquisition batch

One box per batch per population, with a point per sample.

Separate from the batch diagnostics on purpose. Those ask whether the
*embedding* separates by batch, which is a question about the cells.
This asks whether the *reported numbers* differ by batch, which is a
question about the result: a population that steps with acquisition
batch is one whose between-group difference may be an acquisition
difference.

Descriptive only. `batch_group_confounding.csv` reports whether batch
and study group can be told apart at all; where they cannot, no
correction here or elsewhere can separate them.

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

  named vector, sample_id -\> batch.

- outfile:

  path.

- value_col:

  frequency column; defaults to the run's abundance measure.

- ncol:

  panels per row.

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.
