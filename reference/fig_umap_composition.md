# The shared embedding coloured by label, with a legend of composition

Beside the map, one row per label: its colour, its name, its defining
markers, a bar of its share of the embedded cells, and that share split
by study group (up to four; the table carries any number), with a total
row. The catch-all is drawn in its own colour beneath everything else
and listed last.

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
grid of the cells of each type positive for each marker. Every number in
the figure is a percent of ALL embedded cells, so they add up: a row's
group columns to its share, a group column to that group's share of the
cells, a marker tile to at most its row's share, and a marker column to
the share of the embedded cells positive for it. The within-group and
within-type percentages are `pct_<group>` and `pct_pos_<marker>` in the
table.
