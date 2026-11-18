# Learn and validate a gate for each externally supplied label

Learn and validate a gate for each externally supplied label

## Usage

``` r
explain_external_labels(
  cells,
  features,
  min_cells = 200L,
  max_labels = 6L,
  donor_col = "sample_id",
  max_donors = 8L,
  transfer_max_cells = 20000L,
  seed = 42L,
  ...
)
```

## Arguments

- cells:

  embedding cell table carrying `external_label`

- features:

  marker columns to gate on

- min_cells:

  smallest label worth gating

- max_labels:

  ceiling on how many labels to process

- donor_col:

  column holding the donor identity for the transferability split

- max_donors, transfer_max_cells:

  caps passed to
  [`gate_transferability()`](https://bhagesh-h.github.io/cyRAVEN/reference/gate_transferability.md)

- seed:

  seed for the local RNG stream

- ...:

  passed to
  [`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)

## Value

list(summary, polygons, transfer, transfer_summary, strategies), or NULL
