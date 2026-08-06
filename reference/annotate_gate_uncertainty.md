# Compare a between-group difference against the gate uncertainty behind it

Adds two columns to a group comparison table: the typical within-sample
uncertainty on the quantity being compared, and the ratio of the
observed difference to it.

## Usage

``` r
annotate_gate_uncertainty(gstats, ufreq)
```

## Arguments

- gstats:

  a table from
  [`stats_group_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_group_comparison.md)

- ufreq:

  the `frequencies` element of
  [`run_gate_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_gate_uncertainty.md)

## Value

`gstats` with `gate_u_pct_points` and `difference_over_gate_u` appended,
unchanged if the uncertainty table is missing

## Details

READ IT AS A SCREEN, NOT A TEST. The gate uncertainty is partly common
to every sample in a run, since they share a panel, a transform and a
placement rule, so it cancels to some degree in a difference and the
ratio is conservative. A ratio below 1 says the groups differ by less
than the typical distance the cut itself can move, which is a reason to
look at threshold_uncertainty.csv before interpreting the result, not a
p-value.
