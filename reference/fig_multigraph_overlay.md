# FlowJo-style Multigraph Overlay: every cluster x every marker, one peak per group

WHAT THIS IS: a port of FlowJo's Layout Editor "Make Multigraph Overlay
-\> Histograms"
(docs.flowjo.com/flowjo/graphical-reports/graph-options-and-
annotation/le-mgo/), which draws a histogram of every parameter with
each overlaid subset kept in its own colour. Here the grid is CLUSTER
(row) x MARKER (column) and the overlaid subsets are the STUDY GROUPS,
so every panel asks one question: inside this population, does this
marker sit at the same intensity in every cohort?

## Usage

``` r
fig_multigraph_overlay(
  cells,
  outfile,
  markers,
  group_col = "cohort",
  min_cells = 20L,
  panel_label = "",
  dpi = 200,
  subcluster = NULL,
  reference = NULL,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  embedding cell table (population_label, group_col, marker cols)

- outfile:

  Path to write the figure to.

- markers:

  marker columns to draw, normally the UMAP feature set

- group_col:

  Name of the column holding the biological grouping (cohort). Default
  `"cohort"`.

- min_cells:

  a group needs at least this many cells in a population before its
  curve is drawn; below it a KDE is noise shaped like a result.

- panel_label:

  Marker-panel name added to the figure title; empty for none. Default
  `""`.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `200`.

- subcluster:

  Integer vector of subcluster assignments, one per cell.

- reference:

  The group every other group is compared against.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHY IT IS WORTH A FIGURE OF ITS OWN: umap_overview_by_group.png shows
WHERE a population sits and how big it is; this shows WHAT it expresses.
A cohort can carry a normal-sized cluster in a normal position whose
marker intensity has shifted, and no abundance figure in this pipeline
can surface that – they all count cells rather than reading them.

WHY THE PEAKS ARE MODE-NORMALISED (each curve scaled so its own maximum
is 100, which is FlowJo's "Count (%)" convention): the groups differ
several-fold in size, and within a rare population they differ far more.
Plotted on a shared density axis the largest cohort's curve simply
towers over the others and the figure silently becomes a cell-count
comparison – which population_frequencies.png and group_comparison.png
already do properly, with statistics. Normalising to the mode puts every
group on the same footing so the eye compares SHAPE and POSITION, the
only thing this figure is for. Abundance is deliberately somebody else's
job.

MEMORY: the curves are precomputed with stats::density() rather than
handing ggplot the raw per-cell values, for exactly the reason
documented above fig_gating_qc – a (population x marker x group) grid of
raw cells rebuilds the multi-million-row frame that OOM-killed this
pipeline. Do not "simplify" this into geom_density() either.
