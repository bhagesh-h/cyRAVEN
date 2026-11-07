# Choose a threshold by trying every method and scoring the cuts alike

Choose a threshold by trying every method and scoring the cuts alike

## Usage

``` r
best_threshold(
  x,
  ladder = c(4, 6, 8, 12, 16, 24, 32, 48, 64),
  bins = 400L,
  min_quality = 0.05,
  fallback_q = 0.9,
  floor = -Inf,
  prefer = c("deepest", "lowest"),
  tail = c("upper", "lower")
)
```

## Arguments

- x:

  transformed values for one marker within the parent gate.

- ladder:

  smoothing windows for the mindensity sweep, widest last.

- bins:

  histogram resolution for the mindensity sweep. Default `400L`.

- min_quality:

  a candidate must sit in a trough at least this deep to be preferred
  over the tail rule. Default `0.05`.

- fallback_q:

  quantile used only when nothing else is usable.

- floor:

  no candidate at or below this value is considered; see
  [`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md).
  Default `-Inf`.

- prefer:

  `"deepest"` (default) scores every candidate by the trough it sits in.
  `"lowest"`, for a gate that removes only the negative mode, takes the
  lowest density valley that sits in a real trough, and falls back to
  the scored choice when no valley does; see
  [`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md).

- tail:

  `"upper"` (default, also meaning nothing was declared) offers a tail
  gate on either side of the dominant mode and lets the trough decide;
  `"lower"`, the analyst's declaration that the dominant mode is
  positive (`thresholds: <marker>: {tail: lower}`), offers only the
  lower one and takes it as the last resort too, as the default path
  does. Without this the declaration was ignored under
  `--adaptive-gates`.

## Value

`list(threshold, source, quality)`, where `source` names the method that
won so a reader can see how each cut was placed.
