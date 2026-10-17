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
  flag_ratio = 1.25
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

## Value

list(pairs, receivers), each a data.frame, or NULL
