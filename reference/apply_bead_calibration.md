# Apply a fitted calibration to a sample's linear channel values

Only channels whose fit was accepted are converted. A channel that was
not calibrated is returned untouched, and the caller must not present
the two as the same units.

## Usage

``` r
apply_bead_calibration(exprs, channel_cols, calib, channel_names = NULL)
```

## Arguments

- exprs:

  a sample's expression matrix, linear units

- channel_cols:

  named integer vector, marker to column

- calib:

  the table from
  [`fit_bead_calibration()`](https://bhagesh-h.github.io/cyRAVEN/reference/fit_bead_calibration.md)

- channel_names:

  the sample's `$PnN` values, as read, or NULL. Where given and `calib`
  records the detector each fit came from, the fit is applied to the
  sample's channel of that detector.

## Value

list(exprs = converted matrix, applied = character vector of markers)
