# Heatmap of clinical association across populations or markers

WHAT IS ON THE TILE. The signed effect where one exists – Spearman's
rho, or Cliff's delta – so the colour carries direction as well as
strength. A Kruskal-Wallis row has no signed effect and is drawn grey
with its p-value, because a three-level variable has no single direction
to show.

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

  data frame from
  [`clin_associate()`](https://bhagesh-h.github.io/cyRAVEN/reference/clin_associate.md).

- key_col:

  "population" or "marker".

- outfile:

  path.

- title:

  figure title.

- dpi:

  resolution.

- colors:

  palette.

## Details

Significance is marked on the tile rather than encoded in the colour, so
a large effect that did not survive correction still reads as large. On
a cohort this size that distinction is most of the message.
