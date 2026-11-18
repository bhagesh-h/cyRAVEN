# Decide whether a hierarchy gate rests on evidence or on the fallback constant

Separated from
[`apply_gate_hierarchy()`](https://bhagesh-h.github.io/cyRAVEN/reference/apply_gate_hierarchy.md)
so the decision is testable on its own and so the reason travels with it
as text, rather than being implied by a threshold silently becoming NA.

## Usage

``` r
gate_needs_autofix(
  source,
  enabled = FALSE,
  kept = NA_real_,
  parent = NA_real_,
  min_retention = 0.05
)
```

## Arguments

- source:

  The `source` string
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
  returned.

- enabled:

  Whether `--auto-fix-gates` is set.

## Value

`list(skip, reason)`. `skip = TRUE` means apply no cut at this gate.
