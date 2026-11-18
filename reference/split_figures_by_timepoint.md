# Re-emit the whole figure set once per timepoint

Re-emit the whole figure set once per timepoint

## Usage

``` r
split_figures_by_timepoint(ctx, outdir)
```

## Arguments

- ctx:

  Shared artefacts; see
  [`emit_timepoint_figures()`](https://bhagesh-h.github.io/cyRAVEN/reference/emit_timepoint_figures.md).
  Must also carry `smap`.

- outdir:

  Run output directory; figures go to `by_timepoint/<level>/`.

## Value

Character vector of files written, invisibly.
