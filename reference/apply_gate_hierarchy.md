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
  transform = NULL,
  overrides = NULL,
  autofix = FALSE,
  adaptive = FALSE
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

- overrides:

  Optional per-marker override entries for THIS sample, as returned by
  indexing the config's `sample_overrides` block by sample id. See
  [`sample_override()`](https://bhagesh-h.github.io/cyRAVEN/reference/sample_override.md).
  Absent by default, in which case every threshold is derived exactly as
  before.

- autofix:

  Skip a hierarchy gate whose threshold came from `quantile_fallback`,
  carrying its parent through unchanged. Such a gate has its retention
  decided by the fallback constant rather than by the data, so it keeps
  a fixed share of the parent whatever the sample holds. `FALSE` by
  default, because skipping a gate changes every count in an affected
  file. See
  [`gate_needs_autofix()`](https://bhagesh-h.github.io/cyRAVEN/reference/gate_needs_autofix.md).

- adaptive:

  Sweep the kernel bandwidth and try the tail and Otsu rules before
  falling back to a quantile, for both hierarchy gates. See
  [`best_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/best_threshold.md).
  `FALSE` by default, because it moves thresholds.

## Value

list of masks, derived geometry, thresholds, and a tidy counts table
