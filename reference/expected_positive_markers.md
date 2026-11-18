# Markers each population's gate definition requires to be POSITIVE

Feeds the heatmap's gate-audit outline. Reads the same `populations:`
spec the scoring uses, so the audit cannot drift from the definitions it
audits. `any_of` members are included: any one of them being bright
satisfies the gate, so all are legitimate places to look.

## Usage

``` r
expected_positive_markers(spec)
```

## Arguments

- spec:

  Population specification mapping population name to marker directions.
  See
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).
