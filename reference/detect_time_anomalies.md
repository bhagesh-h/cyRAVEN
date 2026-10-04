# Flag acquisition intervals whose rate or signal departs from the file

Flag acquisition intervals whose rate or signal departs from the file

## Usage

``` r
detect_time_anomalies(
  time,
  mat,
  n_bins = 40L,
  mad_k = 5,
  min_events_per_bin = 50L
)
```

## Arguments

- time:

  numeric acquisition time per event, in the file's own units

- mat:

  matrix of channel values on the analysis scale, one column per
  channel, same number of rows as `time`

- n_bins:

  number of equal-width time intervals

- mad_k:

  robust z above which a bin is flagged

- min_events_per_bin:

  bins holding fewer events than this are not judged on their channel
  medians, because a median over a handful of events is noise; they are
  still judged on their rate, which is the point

## Value

list(bins = data.frame, flagged = logical per event, reason)
