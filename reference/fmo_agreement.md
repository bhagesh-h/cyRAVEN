# Distance between a derived cut and its FMO-anchored equivalent

The diagnostic the feature exists for. Positive means the sample's own
density put the cut above the control; negative means below.

## Usage

``` r
fmo_agreement(
  thr_all,
  fmo_thresholds,
  unc = NULL,
  agree_at = 1,
  disagree_at = 3
)
```

## Arguments

- thr_all:

  the thresholds table, carrying `sample_id`, `marker`, `threshold` and
  `source`

- fmo_thresholds:

  data.frame of `sample_id`, `marker`, `fmo_threshold`, `fmo_sample`,
  and optionally `fmo_verdict` from
  [`control_cut()`](https://bhagesh-h.github.io/cyRAVEN/reference/control_cut.md).
  A row whose `fmo_threshold` is `NA` is a control that was refused, and
  is reported as such rather than compared.

- unc:

  optional `thresholds` element of
  [`run_gate_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_gate_uncertainty.md),
  used to scale the distance. A cut the control supplied has no
  uncertainty of its own, because the events it came from are in a
  separate tube, so those rows carry `distance` but not `distance_in_u`.
  The scaled verdict is available for the rows where the sample's own
  cut was the one applied.

- agree_at:

  distance in uncertainties within which the two agree

- disagree_at:

  distance beyond which they are reported as disagreeing

## Value

a data.frame, or NULL

## Details

HOW TO READ `distance_in_u`. It is the gap expressed in units of that
threshold's own standard uncertainty, from `threshold_uncertainty.csv`.
Within about one, the two methods agree to the precision either can
claim, and the derived cut is corroborated by an independent experiment.
Beyond about three, they disagree by more than either can explain and
one of them is wrong: a derived cut far above the FMO is discarding real
signal, and one far below is calling spillover positive.
