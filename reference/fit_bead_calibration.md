# Fit a channel-units to assigned-units calibration from a bead acquisition

Fit a channel-units to assigned-units calibration from a bead
acquisition

## Usage

``` r
fit_bead_calibration(bead_exprs, channel_cols, assigned, min_r2 = 0.98)
```

## Arguments

- bead_exprs:

  the bead file's expression matrix, linear units

- channel_cols:

  named integer vector, marker or channel name to column

- assigned:

  data.frame with a `marker` column and one column per bead population,
  or a long form carrying `marker`, `peak` and `value`

- min_r2:

  fit quality below which a channel is reported as not calibrated

## Value

a data.frame, one row per channel
