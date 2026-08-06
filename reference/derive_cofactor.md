# Derive the arcsinh cofactor from the data

WHAT: for each marker, bisect the cofactor so the inter-quartile range
of the transformed background bulk equals `target_iqr`; the panel
cofactor is the median across markers. WHY: the cofactor sets how much
the near-zero region is expanded. Spectrally unmixed data has a wide
negative/background band; too small a cofactor (the template's default
of 5) expands that noise until it dominates all distances. On the test
batch this recovers ~150 for the baseline panel – derived, not assumed,
so a different instrument or gain setting gets its own value.

## Usage

``` r
derive_cofactor(
  ex,
  marker_cols,
  target_iqr = 0.5,
  lo = 1,
  hi = 5000,
  max_markers = 40L
)
```

## Arguments

- ex:

  Numeric expression matrix, events x channels.

- marker_cols:

  Named integer vector mapping marker symbol to column index.

- target_iqr:

  Interquartile range the transformed background is bisected towards.
  Default `0.5`.

- lo:

  Lower bound of the search interval. Default `1`.

- hi:

  Upper bound of the search interval. Default `5000`.

- max_markers:

  Ceiling on the number of markers used. Default `40L`.
