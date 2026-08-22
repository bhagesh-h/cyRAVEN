# Heatmap of clinical association across populations or markers

The signed effect on each tile – Spearman's rho, or Cliff's delta – so
the colour carries direction as well as strength. A Kruskal-Wallis row
has no signed effect and is drawn grey with its p-value, because a
three-level variable has no single direction to show.

Significance is marked on the tile rather than encoded in the colour, so
a large effect that did not survive correction still reads as large. On
a small cohort that distinction is most of the message.

## Usage

``` r
fig_clinical_heatmap(
  assoc,
  key_col,
  outfile,
  title = "Clinical association",
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- assoc:

  data frame of associations, as returned inside
  [`stats_clinical_association`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_clinical_association.md).

- key_col:

  "population" or "marker".

- outfile:

  path.

- title:

  figure title.

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.
