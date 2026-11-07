# Save the complete analysis state to an .RData file

WHY: reading and gating a multi-hundred-megabyte batch is the expensive
part of this pipeline. Saving the state means re-plotting,
re-thresholding or re-embedding never requires re-reading the FCS files.

## Usage

``` r
save_session(path, state, keep_exprs = FALSE)
```

## Arguments

- path:

  File path.

- state:

  The state.

- keep_exprs:

  The keep exprs. Default `FALSE`.

## Details

The raw expression matrices are dropped by default (they are large and
re-readable from the FCS files); everything derived – masks, thresholds,
coordinates, tables, gate geometry – is kept. Pass keep_exprs = TRUE to
include them for a fully self-contained but much larger file.
