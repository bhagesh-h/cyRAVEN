# Resolve a functional block's scope to concrete population names Precedence: explicit `populations`, else `require`, else all minus `exclude`.

Resolve a functional block's scope to concrete population names
Precedence: explicit `populations`, else `require`, else all minus
`exclude`.

## Usage

``` r
resolve_block_populations(block, spec, available)
```

## Arguments

- block:

  The block.

- spec:

  Population specification mapping population name to marker directions.
  See
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).

- available:

  The available.
