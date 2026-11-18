# The live share of single cells, where a viability gate was applied

What the gate kept, so a sample whose dead cells form the dye's tallest
mode reads as mostly live (see
[`live_mode()`](https://bhagesh-h.github.io/cyRAVEN/reference/live_mode.md)).
The gate's `live_pct_if_dominant_dead` is the second reading for that
case.

## Usage

``` r
live_share(gate)
```

## Arguments

- gate:

  a gate result from
  [`apply_gate_hierarchy()`](https://bhagesh-h.github.io/cyRAVEN/reference/apply_gate_hierarchy.md).

## Value

percent, or NA where no viability gate was applied.
