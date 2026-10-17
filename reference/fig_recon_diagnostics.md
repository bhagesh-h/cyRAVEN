# Initial reconnaissance / QC diagnostic figure – run on EVERY batch

WHAT: a three-row diagnostic across all input files: row 1 log10 FSC-A
vs log10 SSC-A occupancy with the derived leukocyte gate row 2 CD45
(asinh) vs log10 SSC-A with the derived CD45 cutoff row 3 overlaid
density of the lineage markers within the CD45+ parent WHY: this is the
figure that catches the three failure modes that silently ruin a run:
(a) gating debris instead of leukocytes – visible as a dense low-FSC
blob sitting at CD45 background, (b) a mis-placed CD45 cutoff, and (c)
unstained or failed-staining files, visible as marker densities
collapsed into a single background peak. It must be inspected before any
downstream result is trusted.

## Usage

``` r
fig_recon_diagnostics(
  recon,
  outfile,
  max_points = 40000L,
  dpi = 300,
  dens_markers = NULL,
  colors = fcs_colors()
)
```

## Arguments

- recon:

  list of per-file diagnostic records, each a list with: sample_id,
  panel, fsc_log10, ssc_log10 (numeric vectors, may be subsampled), cd45
  (asinh vector or NULL), gate (list with fsc_lo/fsc_hi/ssc_lo/ssc_hi),
  cd45_threshold (numeric or NA), marker_densities (named list of asinh
  vectors within the CD45+ parent), verdict (character, from staining
  QC).

- outfile:

  Path to write the figure to.

- max_points:

  The max points. Default `40000L`.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `300`.

- dens_markers:

  The dens markers.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

the composed patchwork object (also written to `outfile`).
