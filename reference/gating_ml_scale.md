# The Gating-ML 2.0 form of a cyRAVEN transform

The Gating-ML 2.0 form of a cyRAVEN transform

## Usage

``` r
gating_ml_scale(transform, M = 4.5)
```

## Arguments

- transform:

  the transform object from
  [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md).

- M:

  decades for the asinh parameters.

## Value

a function of the marker name returning `list(element, params, unit)`,
where a value on the analysis scale divided by `unit` is the Gating-ML
transformed value; or NULL when the transform has no Gating-ML form
here.
