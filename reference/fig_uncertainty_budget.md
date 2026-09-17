# Where each population's uncertainty comes from

Stacked contribution per term, median across samples, in percentage
points. The terms are the markers the population's definition reads,
plus the CD45 parent gate, which enters every population at once because
it sets the denominator.

## Usage

``` r
fig_uncertainty_budget(budget, outfile, dpi = 200, colors = fcs_colors())
```

## Arguments

- budget:

  the `budget` element of
  [`run_gate_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_gate_uncertainty.md)

- outfile:

  Path to write the figure to.

- dpi:

  Resolution in dots per inch. Default `200`.

- colors:

  Named list of colours; defaults to the package palette. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHY IT IS WORTH A FIGURE OF ITS OWN. It says which gate to fix. A
population whose budget is dominated by one marker has one threshold
worth inspecting in gating_qc.png; a cohort whose budgets are dominated
by the parent term has a CD45 gate problem that no amount of attention
to the downstream markers will improve. Operator studies of manual
gating find the same concentration at the first gate, which is the one
place this figure and that literature can be read against each other.
