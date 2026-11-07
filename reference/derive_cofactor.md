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
  max_markers = 40L,
  rows = NULL
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

- rows:

  Logical or integer index of the events to estimate from, or NULL for
  every event. The pipeline passes the scatter- and singlet-gated cells,
  so debris does not set the background.

## Details

WHAT THE CRITERION AMOUNTS TO. For a marker whose events are mostly
negative, N(0, sigma), the background band is the 5th to 60th percentile
and its IQR is about 0.79 sigma, so an IQR of 0.5 after asinh(x / c)
puts c near 1.6 sigma: the linear region of the transform just wider
than the negative population. That is the condition that keeps a
symmetric negative population one mode, which it stays exactly while c
\>= sigma, and the same idea as the logicle width rule (Parks et al.
2006). It is a heuristic, not a published criterion.

PER MARKER, FROM THE NEGATIVE POPULATION. Detectors differ many-fold in
how wide their negatives are: on one spectral panel the 5th percentile
of the values below zero ran from 147 (CCR5) to 2,330 (CCR9). One
cofactor then splits the widest negatives into two peaks and compresses
the narrowest, and the published recommendation is a cofactor per
channel (den Braanker et al. 2021; flowVS, Azad et al. 2016). The bulk
criterion above cannot be applied per marker as it stands, because for a
marker most cells carry (CD45, CD3) the 5th-60th percentile band is the
POSITIVE population, which gave CD45 a cofactor of 3,256 on the same
panel. So each marker's own negative population is reconstructed from
its values below zero, mirrored about zero, which positive cells never
reach; the same bisection is run on that. Returned as attribute
`per_marker_negative`, NA where fewer than 100 values sit below zero
(uncompensated or linearised log data, for instance).
