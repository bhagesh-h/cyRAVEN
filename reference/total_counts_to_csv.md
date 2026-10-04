# Normalize any –total-counts input down to one CSV on disk

Mirrors
[`absolute_counts_to_csv()`](https://bhagesh-h.github.io/cyRAVEN/reference/absolute_counts_to_csv.md)
deliberately: same reasons (one irregular-layout scan, written once,
into outdir so the user can open exactly what was parsed), and the same
CSV convention every other tabular input to this pipeline already uses.

## Usage

``` r
total_counts_to_csv(path, outdir)
```

## Arguments

- path:

  File path.

- outdir:

  Directory to write outputs to.
