# Fit a channel-units to assigned-units calibration from a bead acquisition

A BEAD FILE IS NAMED BY DETECTOR. Beads carry no antibody, so a bead
file usually has no `$PnS` and its channels resolve to detector names
(`FITC-A`) where the samples resolve to markers (`CD3`). The
assigned-value table may therefore name either: a marker symbol, or a
detector as `$PnN` writes it. With `channel_names` supplied, a detector
name is looked up there, and the detector each row was fitted on is
returned in `channel`, so the fit is applied to the samples by detector
(`$PnN` names the detector, Spidlen et al. 2010) rather than by a stain
name the bead file does not have.

## Usage

``` r
fit_bead_calibration(
  bead_exprs,
  channel_cols,
  assigned,
  min_r2 = 0.98,
  channel_names = NULL
)
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

- channel_names:

  the bead file's `$PnN` values, as read, or NULL.

## Value

a data.frame, one row per channel
