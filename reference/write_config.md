# Write the derived values out as an editable YAML config

Every threshold carries its source and a review flag. The config is an
OVERRIDE mechanism: the script derives everything at runtime and runs
without it. Values under `derived_from_batch` document what THIS batch
produced – editing them changes future runs; deleting the file changes
nothing.

## Usage

``` r
write_config(
  path,
  derived,
  cofactors,
  spec = default_population_spec(),
  blocks = default_functional_blocks(),
  colmap = default_column_map(),
  valmap = default_value_map(),
  colors = default_colors()
)
```

## Arguments

- path:

  File path.

- derived:

  The derived.

- cofactors:

  The cofactors.

- spec:

  Population specification mapping population name to marker directions.
  See
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).
  Default
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).

- blocks:

  Functional-marker blocks. See
  [`default_functional_blocks()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_functional_blocks.md).
  Default
  [`default_functional_blocks()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_functional_blocks.md).

- colmap:

  The colmap. Default
  [`default_column_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_column_map.md).

- valmap:

  The valmap. Default
  [`default_value_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_value_map.md).

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`default_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_colors.md).
