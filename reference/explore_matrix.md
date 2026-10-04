# Transformed matrix over every eligible channel for one sample

Mirrors how the declared path builds `tmat`, with one difference: it is
not restricted to the channels a specification happens to reference.
Fluorescence goes through the panel transform; scatter is log10, because
it is strictly positive and spans decades and the arcsinh cofactor
derived for fluorescence is meaningless on it.

## Usage

``` r
explore_matrix(rd, feats, tr, rows = NULL)
```

## Arguments

- rd:

  One element of `reads`.

- feats:

  Feature names wanted.

- tr:

  Transform object from
  [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md).

## Value

Numeric matrix, cells x features.
