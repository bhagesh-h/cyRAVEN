# Derive the upper bound of an "intermediate" gate

WHAT: among the marker-positive cells, find the valley separating the
intermediate mode from the bright mode; that valley is the upper bound.
WHY: "CD14 int CD16 int" is a genuine three-level gate (negative /
intermediate / bright), and monocyte subsetting depends on it. Returns
NA when the positive fraction has no second valley – the population is
then reported UNAVAILABLE rather than being silently merged with the
bright subset.

## Usage

``` r
derive_intermediate_bounds(tmat, thr, parent, spec)
```

## Arguments

- tmat:

  The tmat.

- thr:

  The thr.

- parent:

  The parent.

- spec:

  Population specification mapping population name to marker directions.
  See
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).
