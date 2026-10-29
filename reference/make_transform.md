# Build the intensity transform the whole run uses

WHAT IT RETURNS: an object with `fn(x, marker)`, so every call site
applies the transform the same way and none of them needs to know which
method is in force. Adding a method here reaches the gating, scoring,
embedding and figure code without touching any of them.

## Usage

``` r
make_transform(
  method = c("arcsinh", "logicle", "none"),
  cofactor = NULL,
  logicle = NULL,
  per_marker = NULL
)
```

## Arguments

- method:

  "arcsinh", "logicle" or "none"

- cofactor:

  arcsinh cofactor (required for method = "arcsinh")

- logicle:

  named list of per-marker parameters from
  [`derive_logicle_pooled()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_logicle_pooled.md)
  (required for method = "logicle")

- per_marker:

  named numeric vector of arcsinh cofactors by marker, or NULL. A marker
  named here is transformed with its own cofactor; any other, and
  scatter, with `cofactor`. See
  [`derive_cofactor()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_cofactor.md).

## Value

list(method, fn, params, label)
