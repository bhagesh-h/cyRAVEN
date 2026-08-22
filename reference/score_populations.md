# Evaluate the population spec against one file's thresholded markers

Evaluate the population spec against one file's thresholded markers

## Usage

``` r
score_populations(
  tmat,
  thr,
  parent,
  spec = default_population_spec(),
  hi_thr = NULL
)
```

## Arguments

- tmat:

  matrix of transformed marker values (cells x markers, named cols)

- thr:

  named numeric vector of thresholds

- parent:

  logical mask (the CD45+ parent)

- spec:

  Population specification mapping population name to marker directions.
  See
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).
  Default
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).

- hi_thr:

  named numeric: upper bound for "intermediate" requirements. Derived by
  derive_intermediate_bounds(); a marker requested as "intermediate"
  without an upper bound makes its population UNAVAILABLE rather than
  silently collapsing to "above".

## Value

list(masks, labels, unavailable)
