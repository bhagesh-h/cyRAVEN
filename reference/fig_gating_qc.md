# Per-gate QC figure: every applied gate shown per sample with its threshold

WHY THIS IS A DELIVERABLE: a gated analysis is only as trustworthy as
its gates, and a threshold that is numerically plausible can still sit
in the wrong place. Showing every cut on the distribution it was derived
from is the only way a reviewer can check the gating without re-running
anything. DO NOT "SIMPLIFY" THIS BACK TO geom_density() ON RAW CELLS. It
was written that way and it OOM-killed the container mid-run.

## Usage

``` r
fig_gating_qc(recon, outfile, dpi = 200, colors = fcs_colors())
```

## Arguments

- recon:

  The recon.

- outfile:

  Path to write the figure to.

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

The obvious implementation hands every CD45+ cell to geom_density() and
lets ggplot compute the curve. That builds one row per CELL per MARKER:
with 25 samples x 13 markers x ~813k CD45+ cells (–max-events-per-file
300000) it is 10,572,757 rows through do.call(rbind, ...), on top of the
~1.9 GB of reads\$exprs + pops\$tmat that main() is still holding at
this point. Measured under –memory=6g: SIGKILL, exit 137. The failure is
SILENT – the OOM killer leaves no R error, no traceback, nothing in the
log. The run simply stops after STEP 4 with recon_diagnostics.png
written and gating_qc.png missing. It survived at –max-events-per-file
120000 purely by luck of scale.

The cure is NOT to subsample the cells. Measured on this study's own
data, thinning to 5,000 cells per panel distorts the drawn curve by up
to 21% of its peak height and makes a 2%-of-CD45+ population – a
perfectly ordinary size for Vd1/Vd2, NKT or dendritic cells here –
invisible in 19 of 20 draws. A gating-QC figure that silently drops the
rare populations whose gates most need checking is worse than no figure.

So the curve is computed here, from EVERY cell, and only the curve is
passed to ggplot: stats::density() reduces each (sample, marker) to 512
points, which is what geom_density would have drawn anyway (same default
bw.nrd0 bandwidth, same 512-point grid). 10.6M rows become 325 x 512 =
166k, a ~64x reduction, with NO loss of fidelity and no subsampling. It
is also RNG-free, so unlike a sampling-based fix it cannot shift the
.Random.seed stream that STEP 6's UMAP cell selection (sample(), line
~3559) draws from – the embedding is bit-identical to before this
change.
