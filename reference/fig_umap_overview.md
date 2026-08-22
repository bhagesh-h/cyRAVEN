# Population + sample + covariate colouring panel set

WHAT: writes the main UMAP figure (population identity, sample, and any
patient covariates present) for one panel group. WHY: one shared
embedding viewed many ways is what makes populations comparable across
patients; separate embeddings per sample would give coordinates that
cannot be compared.

## Usage

``` r
fig_umap_overview(
  cells,
  outfile,
  panel_label = "",
  covariates = NULL,
  feature_cols = character(0),
  width = 15,
  height = 11,
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  data.frame with umap_1, umap_2, population_label, sample_id and
  optionally covariate columns.

- outfile:

  Path to write the figure to.

- panel_label:

  Marker-panel name added to the figure title; empty for none. Default
  `""`.

- covariates:

  NULL (default) to DISCOVER plottable covariates from the joined table,
  or an explicit character vector to fix the set.

  WHY discovery is the default: a fixed list can only plot the
  covariates someone anticipated. The single most important variable in
  a case-control study is the group column – cohort, genotype, treatment
  arm – and its name is study-specific, so a hardcoded list silently
  omits exactly the panel the experiment was designed to produce.
  Anything that arrives from the sample map or patient table and varies
  across cells is plottable; the exclusion list below is structural
  (coordinates, gate bookkeeping, marker intensities), because those are
  plotted by other figures or are not covariates at all.

- feature_cols:

  The feature cols. Default `character(0)`.

- width:

  Figure width in inches. Default `15`.

- height:

  Figure height in inches. Default `11`.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `300`.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

the composed patchwork object (also written to `outfile`).
