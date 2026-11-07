# Per-sample (or per-group) density comparison over the shared embedding

WHY: side-by-side 2D density of the SAME embedding shows where each
sample or group gains or loses cells – the comparison the shared
embedding exists for. Called with `facet_by = "sample_id"` (the default)
for per-sample QC, and again with the resolved group column when 2+
groups exist, so cohorts can be compared the same way sample_id already
is.

## Usage

``` r
fig_density_by_sample(
  cells,
  outfile,
  panel_label = "",
  facet_by = "sample_id",
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`,
  `sample_id`, `population_label` and one column per marker.

- outfile:

  Path to write the figure to.

- panel_label:

  Marker-panel name added to the figure title; empty for none. Default
  `""`.

- facet_by:

  any column present in `cells` – sample_id, cohort, or any other
  covariate the patient table carried through.

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
