# The same polygons as a plain table in linear units

The format anyone can read without an XML parser, and the one to use
when redrawing the gate by hand at the instrument. Carries both scales
side by side, so a vertex can be checked against the figure it came
from.

## Usage

``` r
polygons_linear_table(
  polygons,
  transform,
  id_col = "label",
  x_col = "x_transformed",
  y_col = "y_transformed",
  n_per_edge = 24L
)
```

## Arguments

- polygons:

  a polygons table from
  [`explain_external_labels()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_external_labels.md)
  or
  [`explain_unmatched_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_unmatched_clusters.md)

- transform:

  the transform object the run used

- id_col:

  column naming the strategy each row belongs to

- x_col, y_col:

  columns holding the vertex coordinates

- n_per_edge:

  edge subdivision passed to
  [`polygon_to_linear()`](https://bhagesh-h.github.io/cyRAVEN/reference/polygon_to_linear.md)

## Value

a data.frame, or NULL
