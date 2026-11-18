# Sample ids the sample sheet declares as controls

Used wherever a cohort-wide quantity is pooled over samples – the
transform, and the peer cuts and peer comparison of thresholds – so that
a control tube is gated like a sample but never counted as one of the
cohort.

## Usage

``` r
declared_controls(smap = NULL)
```

## Arguments

- smap:

  the sample map, or NULL.

## Value

character vector of sample ids, empty when none are declared.
