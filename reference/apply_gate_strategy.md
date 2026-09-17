# Apply a learned gating strategy to new cells

Runs the polygons in order, each on the cells the previous one kept.
This is what makes a proposal executable: the same function scores the
cells the gate was fitted on, a held-out donor, and next year's cohort.

## Usage

``` r
apply_gate_strategy(X, strategy, depth = NULL)
```

## Arguments

- X:

  marker matrix carrying at least the columns each level names

- strategy:

  the list returned by
  [`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)

- depth:

  stop after this many levels; NULL uses all of them

## Value

logical vector, TRUE for cells inside every gate
