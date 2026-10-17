# QC heatmap for –absolute-counts: one tile per sample x population

WHAT: log10(cells/uL) as tile colour, one row per matched sample, one
column per population, ordered by group where available. Blank tiles are
values genuinely absent from the source file for that sample x
population, not zero.

## Usage

``` r
fig_absolute_counts_qc(
  ac,
  outfile,
  group_of = NULL,
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- ac:

  The ac.

- outfile:

  Path to write the figure to.

- group_of:

  Named character vector mapping sample_id to group label.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `200`.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHY THIS COMES BEFORE THE GROUP-COMPARISON FIGURE (absolute_counts.png):
this is externally supplied data this pipeline did not measure – a
transcription slip, a stray order-of-magnitude error, or a sample that
silently failed to match (see the NOTE lines the loader prints) will not
look anomalous in a single boxplot the way it jumps out in a full grid.
Controls and QC-failed samples are labelled rather than dropped, on the
same principle fig_gating_qc() uses: a QC figure that hides the excluded
rows cannot be used to catch a problem with them.
