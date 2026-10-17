# Panel borders and gutters for any faceted figure

WHY THIS IS NOT IN
[`theme_cyto()`](https://bhagesh-h.github.io/cyRAVEN/reference/theme_cyto.md):
putting a border on every panel would box the single-panel figures too,
where there is nothing to divide and the frame is just ink. This is for
figures whose panels are meant to be COMPARED.

## Usage

``` r
theme_panel_borders(colors = fcs_colors(), spacing = 0.7)
```

## Arguments

- colors:

  Named list of colours; defaults to the package palette. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

- spacing:

  Gutter between panels, in lines. Default `0.7`.

## Details

WHY A BORDER AT ALL. Two scatters side by side with nothing between them
read as one wide cloud: the eye has to find the boundary from a gap in
the point density, and where a facet's cells reach the edge of their
panel there is no gap to find. A thin border and a wider gutter make the
division explicit at a glance. `umap_<marker>_by_<category>.png` has
always done this; every other multi-panel figure in the package now does
it through this one function, so a report cannot contain framed and
unframed versions of the same kind of figure.
