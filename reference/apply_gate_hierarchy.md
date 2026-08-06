# Apply the full gate hierarchy to one file

Chain: all events -\> leukocytes (scatter) -\> single cells -\> live -\>
CD45+. Each step degrades gracefully: a missing viability dye or missing
CD45 skips that step with a warning rather than aborting, because a
panel that lacks one is still analysable.

## Usage

``` r
apply_gate_hierarchy(
  rd,
  cofactor,
  cfg = list(),
  control_ref = NULL,
  singlet_k = 3,
  viability_name = NULL,
  cd45_name = "CD45",
  transform = NULL
)
```

## Arguments

- rd:

  The rd.

- cofactor:

  Arcsinh cofactor. See
  [`derive_cofactor()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_cofactor.md).

- cfg:

  The cfg. Default [`list()`](https://rdrr.io/r/base/list.html).

- control_ref:

  The control ref.

- singlet_k:

  The singlet k. Default `3`.

- viability_name:

  The viability name.

- cd45_name:

  The cd45 name. Default `"CD45"`.

- transform:

  Intensity transform from
  [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md).
  Defaults to arcsinh with `cofactor`, which is what every caller got
  before the transform became selectable.

## Value

list of masks, derived geometry, thresholds, and a tidy counts table
