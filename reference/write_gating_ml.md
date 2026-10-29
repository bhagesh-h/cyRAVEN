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
  channel_map = NULL,
  space = c("transformed", "linear"),
  compensation_ref = c("FCS", "uncompensated", "supplied"),
  spillover = NULL,
  M = 4.5
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
  [`polygon_to_linear()`](https://bhagesh-h.github.io/cyRAVEN/reference/polygon_to_linear.md),
  for `space = "linear"`

- channel_map:

  optional named character vector, marker symbol to the files' `$PnN`;
  required with `compensation_ref = "supplied"`

- space:

  "transformed" (default where the transform has a Gating-ML form) or
  "linear"

- compensation_ref:

  "FCS" (default), "uncompensated" or "supplied"; see above

- spillover:

  for `compensation_ref = "supplied"`, the spillover matrix the files
  were compensated with, the files' `$PnN` as its column names (see
  [`gating_ml_spillover()`](https://bhagesh-h.github.io/cyRAVEN/reference/gating_ml_spillover.md))

- M:

  decades for the Gating-ML asinh parameters. Any positive value states
  the same transform; 4.5 matches the logicle default.

## Value

`path`, invisibly, or NULL when there was nothing to write

## Details

COORDINATES ARE ON THE ANALYSIS SCALE by default: one Gating-ML
transformation per marker (`fasinh` for arcsinh, `logicle` for logicle)
and each vertex divided into its units, which states the fitted region
exactly. `space = "linear"` writes the region in the units the FCS file
stores instead, edges subdivided before inversion.

COMPENSATION. A gate drawn on data compensated with the file's own
`$SPILLOVER` is declared `compensation-ref="FCS"`, which tells the
reader to apply that same matrix; data with no matrix in the file,
already unmixed by the acquisition software, are declared
`uncompensated`. The earlier export declared every gate `uncompensated`,
which told readers to gate the stored, uncompensated values for gates
cyRAVEN had placed on compensated ones. Data compensated with a matrix
from `--spillover` (`"supplied"`) carry that matrix in the document as a
`spectrumMatrix`, and each gate dimension it compensated references its
identifier, as Gating-ML 2.0 section 7 defines: `"FCS"` would tell the
reader to apply each file's own `$SPILLOVER`, the matrix the run
replaced. A dimension whose detector the matrix does not hold was never
compensated and is declared `uncompensated`.

READING "FCS" WITH flowCore. A reader applying an `"FCS"` gate fetches
each file's `$SPILLOVER` itself. flowCore 2.18 to 2.22, every release up
to Bioconductor 3.22, stops there on R 4.2 or later ("the condition has
length \> 1": `getSpilloverFromFlowFrame()` tests
`class(mySpill) != "matrix"`, and a matrix has two classes since R 4.0);
the fix (`!is(mySpill, "matrix")`, upstream commit 903b5b3) is in
flowCore 2.25 (Bioconductor devel), and the cyRAVEN image installs
2.18.0 with it (`inst/scripts/patch_flowcore.R`). Checked with flowUtils
1.63.0: on the patched 2.18 the three cases apply and count exactly the
events computed directly.

DIMENSION NAMES. Gating-ML 2.0 names a dimension of an FCS file by its
`$PnN` (section 4.1), also under `"FCS"` (4.2.2), so `channel_map`
should map each marker symbol to the files' `$PnN`
([`gating_ml_markers()`](https://bhagesh-h.github.io/cyRAVEN/reference/gating_ml_markers.md));
without it the marker symbols are written, which no FCS file holds as a
column name. Under a `spectrumMatrix` a compensated dimension is named
by its fluorochrome instead (7.3), here the marker symbol, prefixed
`Comp-` where it would equal a detector name, which section 7.2 forbids;
a detector carrying no marker of the run is `Comp-` and its `$PnN`.
