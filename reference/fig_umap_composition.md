# The shared embedding coloured by label, with a legend of composition

Beside the map, one row per label: its colour, its name, its defining
markers, a bar of its share of the embedded cells, and its share within
each study group (up to four; the table carries any number). The
catch-all is drawn in its own colour beneath everything else and listed
last.

## Usage

``` r
fig_umap_composition(
  cells,
  comp,
  outfile,
  label_col = "population_label",
  title = "Composition of the shared embedding",
  subtitle = NULL,
  header = "Cell type and defining markers",
  xlim = NULL,
  ylim = NULL,
  dpi = 220
)
```

## Arguments

- cells:

  Embedded cell table with `umap_1`, `umap_2` and the label column.

- comp:

  Output of
  [`composition_table()`](https://bhagesh-h.github.io/cyRAVEN/reference/composition_table.md),
  with a `colour` column.

- outfile:

  Destination PNG.

- label_col:

  Label column.

- title, subtitle:

  Figure title and subtitle.

- header:

  Heading over the legend's name column.

- xlim, ylim:

  Optional axis limits for the map, so figures of subsets of one
  embedding (a timepoint each) are drawn on the same extent.

- dpi:

  Resolution.

## Value

The plot, invisibly, or NULL when there was nothing to draw.

## Details

To the right of that, where `comp` carries `pct_pos_<marker>` columns, a
grid of the percent of each cell type positive for each marker, on the
same rows: the composition by cell type and, across the row, by marker.
A cell type's defining markers read near 0 or 100 by construction; the
others are what the figure adds.
