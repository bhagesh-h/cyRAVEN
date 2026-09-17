# Choose a threshold by trying every method and scoring the cuts alike

Choose a threshold by trying every method and scoring the cuts alike

## Usage

``` r
best_threshold(
  x,
  ladder = c(4, 6, 8, 12, 16, 24, 32, 48, 64),
  bins = 400L,
  min_quality = 0.05,
  fallback_q = 0.9
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

## Value

`list(threshold, source, quality)`, where `source` names the method that
won so a reader can see how each cut was placed.
