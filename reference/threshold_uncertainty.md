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
  grid = valley_setting_grid()
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
