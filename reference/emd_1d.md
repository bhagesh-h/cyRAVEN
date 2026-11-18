# One-dimensional Earth Mover's distance between two samples

The mean absolute difference between the two empirical quantile
functions, evaluated on a common grid. Equivalent to the area between
the two cumulative distributions, and to the Wasserstein-1 distance.

## Usage

``` r
emd_1d(a, b, n = 512L)
```

## Arguments

- a, b:

  numeric vectors

- n:

  quantile grid resolution

## Value

the distance in the units of `a` and `b`, or NA if either side is too
small to describe a distribution
