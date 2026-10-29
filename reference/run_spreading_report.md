# Spreading across a cohort, and which markers it is costing

Aggregates
[`spreading_pairs()`](https://bhagesh-h.github.io/cyRAVEN/reference/spreading_pairs.md)
over samples and joins the result to how each receiver's threshold was
obtained. A marker that receives a great deal of spreading AND falls
back to a quantile is the actionable finding: its cut is unresolved for
a reason the panel can fix.

## Usage

``` r
run_spreading_report(
  pops,
  gates,
  thr_all = NULL,
  panel_of = NULL,
  max_samples = 8L,
  flag_ratio = 1.25,
  trend_flag = 0.3
)
```

## Arguments

- pops:

  per-sample scoring results, carrying `tmat` and `thresholds`

- gates:

  per-sample gate objects

- thr_all:

  the thresholds table, used for the fallback rate per marker

- panel_of:

  named character vector mapping sample to panel

- max_samples:

  cap on samples contributing, since the ranking stabilises long before
  every sample is used

- flag_ratio:

  spreading ratio at or above which a pair is called substantial

- trend_flag:

  absolute rank correlation at or above which a pair is flagged as a
  possible compensation or unmixing error. Default 0.3, a conventional
  moderate correlation, not a validated cut. A negative trend on a pair
  whose receiver is `flag_ratio` times wider at a given source level is
  reported as `trend_unresolved` instead, since truncation at the
  receiver's cut produces one from spreading alone.

## Value

list(pairs, receivers), each a data.frame, or NULL
