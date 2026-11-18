# Rebuild a transform from its stored method and parameters

The inverse of keeping `list(method, params)` from a
[`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md)
object, which is what a baseline and a saved UMAP model store so a later
run can carry its values onto their scale.

## Usage

``` r
rebuild_transform(tr)
```

## Arguments

- tr:

  list with `method` and `params`, or NULL.

## Value

a transform object, or NULL when it cannot be rebuilt.
