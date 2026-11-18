# Per-marker distributional drift between acquisition batches

For each marker, the largest Earth Mover's distance between any two
batches, reported in analysis units and scaled by the marker's pooled
MAD so markers can be compared with each other.

## Usage

``` r
marker_batch_drift(
  cells,
  batch_col,
  markers = NULL,
  max_cells = 20000L,
  min_cells = 200L,
  flag_at = 0.5,
  seed = 42L
)
```

## Arguments

- cells:

  cell-level table carrying `sample_id`, the batch column, and one
  column per marker on the analysis scale

- batch_col:

  name of the batch column in `cells`

- markers:

  marker columns to test; defaults to every numeric column that is not
  structural

- max_cells:

  cells sampled per batch before the distances are computed; the
  quantile function stops moving long before the whole batch is used

- min_cells:

  smallest batch worth comparing

- flag_at:

  value of `emd_over_mad` at or above which a marker is flagged

- seed:

  seed for the local RNG stream, which is restored on exit

## Value

a data.frame ordered by `emd_over_mad`, or NULL

## Details

WHAT A FLAG MEANS. A marker whose distribution differs between batches
by an appreciable fraction of its own spread was measured differently in
them. That is a statement about the assay, not about the donors, unless
batch and study group coincide – which is the question
[`batch_mixing_report()`](https://bhagesh-h.github.io/cyRAVEN/reference/batch_mixing_report.md)
answers with Cramer's V, and which should be read first. Where the two
overlap, this table cannot separate a reagent lot from the biology
either.
