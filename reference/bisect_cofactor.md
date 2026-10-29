# Bisect one marker's cofactor to the background-IQR target

The search
[`derive_cofactor()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_cofactor.md)
runs per marker, on the 5th to 60th percentile of `x`.

## Usage

``` r
bisect_cofactor(x, target_iqr = 0.5, lo = 1, hi = 5000)
```

## Arguments

- x:

  numeric vector, linear units.

- target_iqr, lo, hi:

  see
  [`derive_cofactor()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_cofactor.md).

## Value

the cofactor, or NA when fewer than 100 values fall in the band.
