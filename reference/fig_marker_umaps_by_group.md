# One UMAP per marker, pooled and split by every category present

WHY: the two existing by-group figures answer different questions.
`umap_markers.png` colours by marker intensity but pools the groups, so
a shift that happens in one group only is averaged away.
`umap_overview_by_group.png` splits by group but colours by population
or covariate, never by intensity. Neither answers "is this marker
brighter, or somewhere else, in one group", which is the comparison a
case-control or longitudinal design is usually about.

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

  Data frame of embedded cells with `umap_1`, `umap_2`, the marker
  columns and any category columns.

- markers:

  Character vector of marker columns to draw.

- outdir:

  Directory to write into; created if absent.

- group_col:

  Category column to draw first, normally the run's `--group-column`.
  Further categories are detected automatically. NULL is fine: detection
  still runs. A column resolving to one level is skipped.

- panel_label:

  Marker-panel name added to each title; empty for none. Default `""`.

- dpi:

  Resolution in dots per inch. Default `170`.

- colors:

  Named list of colours; defaults to the package palette. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

Character vector of the files written, invisibly.

## Details

     ONE FILE PER MARKER rather than one markers x groups grid. A grid is
     unreadable past a handful of each, and a single marker's comparison is
     the unit that goes into a slide or a figure panel.

     The colour scale is SHARED across the facets within a file, and clipped
     to the 1st-99th percentile of the cells drawn, so a colour difference
     between two panels of one file is a real intensity difference. It is not
     shared between files, for the reason given in [fig_marker_grid()].

     Panel density is comparable only because the embedding draws the same
     number of cells from every sample; see [plan_subsample()].

     ONE POOLED FILE PLUS ONE PER CATEGORY. `umap_<marker>.png` pools every
     sample; `umap_<marker>_by_<category>.png` splits it. They answer
     different questions -- where the marker is at all, against whether it
     sits differently between categories -- and the pooled one has no other
     home, because each facet of a split figure holds only a subset of the
     cells and `umap_markers.png` shrinks every marker into one cell of a
     grid.

     EVERY category present is drawn, not only the statistical group column.
     A cohort usually carries several -- timepoint, infection focus,
     phenotype -- and naming one of them as `--group-column` says nothing
     about which are worth seeing. See [marker_facet_cols()] for what counts
     as a category and what is excluded.
