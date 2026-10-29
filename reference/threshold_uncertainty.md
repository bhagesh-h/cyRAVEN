# Uncertainty in one per-sample threshold

Returns the standard uncertainty of the cut, in the units of `x`, split
into the component from finite sampling and the component from the
placement settings, plus their quadrature sum.

## Usage

``` r
threshold_uncertainty(
  x,
  source = "valley",
  B = 100L,
  seed = 42L,
  max_events = 20000L,
  fallback_q = 0.9,
  grid = valley_setting_grid(),
  floor = -Inf,
  prefer = "deepest",
  tail = "upper"
)
```

## Arguments

- x:

  parent-gate values for one marker in one sample, on the analysis scale

- source:

  derivation string from
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)

- B:

  bootstrap replicates

- seed:

  seed for the local RNG stream

- max_events:

  cap on events per bootstrap replicate; the cut is a histogram feature
  and stops moving long before the full parent gate is used

- fallback_q:

  the quantile
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
  falls back to

- grid:

  setting combinations for the method component

- floor:

  the floor
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
  placed the cut under, so every resampled cut obeys the same rule.
  Default `-Inf`.

- prefer:

  the valley preference the cut was placed with; see
  [`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md).
  Default `"deepest"`.

- tail:

  the side a `tailgate` cut was placed on, as
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
  recorded it. Default `"upper"`.

## Value

list(u_sampling, u_method, u_combined, rel_depth, valley_rate, n_events,
basis)

## Details

HOW EACH SOURCE IS TREATED. `source` is the string
[`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
recorded for this threshold, and it determines what "uncertain" even
means:

- `config`: the value was declared, not estimated. Uncertainty is zero
  by construction, and saying so is different from saying it is small.

- `valley`: both components are computed as described at the top of this
  file.

- `quantile_fallback`: no valley was found, so there is no valley to
  resample. Sampling uncertainty is the bootstrap spread of the
  quantile. The method component sweeps which quantile, because the
  choice of 0.90 is arbitrary and the spread it produces is normally
  large. That is the correct answer: a fallback threshold is barely
  determined by the data, and this is the number that says so.

- `control_q995` and `control_q995_valley_rejected`: the cut came from a
  separate control tube not passed here, so neither component can be
  computed from `x`. Both are NA with a reason, rather than a small
  number that would read as confidence.

- `tailgate`: no valley, so the tail rule placed the cut. Sampling
  uncertainty is the bootstrap spread of the tail cut; the method
  component sweeps the multiple of the robust SD it is placed at.

- `peer_median`: `--auto-fix-gates` took the cut from the density minima
  of the other samples (see
  [`peer_cut_replacements()`](https://bhagesh-h.github.io/cyRAVEN/reference/peer_cut_replacements.md)),
  so it was not placed on `x` and there is nothing in `x` to resample.
  NA with a reason.
