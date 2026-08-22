# One UMAP per marker, faceted by study group

WHY: the two existing by-group figures answer different questions.
`umap_markers.png` colours by marker intensity but pools the groups, so
a shift that happens in one group only is averaged away.
`umap_overview_by_group.png` splits by group but colours by population
or covariate, never by intensity. Neither answers "is this marker
brighter, or somewhere else, in one group", which is the comparison a
case-control or longitudinal design is usually about.

ONE FILE PER MARKER rather than one markers x groups grid. A grid is
unreadable past a handful of each, and a single marker's comparison is
the unit that goes into a slide or a figure panel.

The colour scale is SHARED across the facets within a file, and clipped
to the 1st-99th percentile of the cells drawn, so a colour difference
between two panels of one file is a real intensity difference. It is not
shared between files, for the reason given in
[`fig_marker_grid()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_marker_grid.md).

Panel density is comparable only because the embedding draws the same
number of cells from every sample; see
[`plan_subsample()`](https://bhagesh-h.github.io/cyRAVEN/reference/plan_subsample.md).

TWO FILES PER MARKER ON A GROUPED RUN, not one. `umap_<marker>.png`
pools every sample; `umap_<marker>_by_<group>.png` splits by group. They
answer different questions – where the marker is at all, against whether
it sits differently in one group – and the first has no other home. Each
facet of the split figure holds only a subset of the cells, and
`umap_markers.png` shrinks every marker into one cell of a grid, so
without the pooled panel no figure shows one marker over the whole
cohort at full size. With no group column resolved, only the pooled
panel is written.

## Usage

``` r
fig_marker_umaps_by_group(
  cells,
  markers,
  outdir,
  group_col = NULL,
  panel_label = "",
  dpi = 170,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  Data frame of embedded cells with `umap_1`, `umap_2` and the marker
  columns, plus the group column when there is one.

- markers:

  Character vector of marker columns to draw.

- outdir:

  Directory to write into; created if absent.

- group_col:

  Name of the column to facet by, or NULL for one unfaceted panel per
  marker. A column resolving to a single level is treated as NULL.

- panel_label:

  Marker-panel name added to each title; empty for none. Default `""`.

- dpi:

  Resolution in dots per inch. Default `170`.

- colors:

  Named list of colours; defaults to the package palette. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

Character vector of the files written, invisibly.
