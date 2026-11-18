# Fraction of each cluster's cells positive for each marker

The heart of what this module adds over a standalone clusterer.
Positivity is judged against each cell's OWN SAMPLE'S threshold, so a
brightly stained sample and a dim one are not compared against a common
cut. Falls back to the pooled per-cluster median when a sample has no
threshold for that marker.

## Usage

``` r
explore_positivity(X, cluster, sample_id, thr_by_sample = list())
```

## Arguments

- X:

  Transformed matrix, cells x features.

- cluster:

  Integer cluster per cell.

- sample_id:

  Character sample id per cell.

- thr_by_sample:

  Named list: sample id -\> named numeric vector of thresholds on the
  transformed scale. May be empty.

## Value

Numeric matrix, clusters x features, values in `[0, 1]`, with attribute
"source" naming how each feature was called.
