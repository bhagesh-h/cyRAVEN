# Write learned gating strategies as a Gating-ML 2.0 document

Each strategy becomes a chain of polygon gates, level *n* declared as
the child of level *n-1*, which is what makes the file a gating STRATEGY
rather than a bag of regions: software reading it applies them in the
order they were learned, on the cells the previous gate kept.

## Usage

``` r
write_gating_ml(
  polygons,
  path,
  transform,
  id_col = "label",
  x_col = "x_transformed",
  y_col = "y_transformed",
  n_per_edge = 24L,
  channel_map = NULL
)
```

## Arguments

- polygons:

  a polygons table from
  [`explain_external_labels()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_external_labels.md)
  or
  [`explain_unmatched_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_unmatched_clusters.md)

- path:

  destination `.xml`

- transform:

  the transform object the run used

- id_col:

  column naming the strategy each row belongs to

- x_col, y_col:

  columns holding the vertex coordinates

- n_per_edge:

  edge subdivision passed to
  [`polygon_to_linear()`](https://bhagesh-h.github.io/cyRAVEN/reference/polygon_to_linear.md)

- channel_map:

  optional named character vector, marker symbol to channel

## Value

`path`, invisibly, or NULL when there was nothing to write

## Details

COORDINATES ARE LINEAR AND UNCOMPENSATED-REFERENCED. The gates are
emitted in the units the FCS file stores, so no transformation element
is needed and there is nothing for the reader to get wrong.
`compensation-ref` is declared as `uncompensated` because cyRAVEN
applies the acquisition spillover matrix from the file's own keyword
block before gating; a reader that applies its own compensation on top
would be compensating twice.

DIMENSION NAMES ARE MARKER SYMBOLS. cyRAVEN resolves channels through
`$PnS`, so a gate names CD3 rather than the detector. Software keyed on
detector names needs `channel_map` to translate.
