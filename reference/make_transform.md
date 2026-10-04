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
  logicle = NULL
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

## Value

list(method, fn, params, label)
