# Test each marker's per-sample threshold for a cohort difference

Test each marker's per-sample threshold for a cohort difference

## Usage

``` r
stats_threshold_drift(thr, group_of, spec = NULL, min_n = 3L)
```

## Arguments

- thr:

  thresholds_used table (sample_id, marker, threshold, source)

- group_of:

  named vector sample_id -\> cohort

- spec:

  population spec, used to name which populations each flagged marker
  feeds into

- min_n:

  Minimum samples per group before a test is attempted. Default `3L`.
