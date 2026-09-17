# Write the declared specification plus a draft entry per uncovered cluster

Write the declared specification plus a draft entry per uncovered
cluster

## Usage

``` r
write_suggested_config(cfg, ex, outfile)
```

## Arguments

- cfg:

  The parsed run config, whose `populations` block is copied through
  unchanged.

- ex:

  What
  [`run_explore()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_explore.md)
  returned: `dir` and `gaps`.

- outfile:

  Destination YAML.

## Value

`list(n_declared, n_added, path)`, or NULL when there is nothing to add.
