# Learn a convex two-marker gate for a labelled set of cells

WHAT: fits an intersection of half-planes that separates `y == 1` cells
from the rest in the plane of two markers, then reports how well it does
on cells held out of the fit.

## Usage

``` r
learn_convex_gate(
  X,
  y,
  n_planes = 8L,
  s_schedule = 800,
  lambda_grid = c(0, 0.5, 1, 2, 4, 8),
  holdout = 0.3,
  n_target_total = NULL,
  max_cells = 20000L,
  seed = 42L
)
```

## Arguments

- X:

  n x 2 numeric matrix of marker values (arcsinh scale)

- y:

  0/1 or logical target indicator, length n

- n_planes:

  total half-planes; 4 are PC-aligned and fixed, the rest free

- s_schedule:

  sigmoid steepness. Supplying several increasing values anneals: each
  fit warm-starts the next. The default is a single value, because
  annealing measured no reliable improvement (+0.002 to +0.007 F1 on the
  cases tried, inside the seed-to-seed spread) for roughly twice the
  runtime. It is left available for a dataset where the fit does get
  stuck

- lambda_grid:

  tightness values to search

- holdout:

  fraction of cells reserved for evaluation, stratified by label

- n_target_total:

  recall denominator; defaults to the targets supplied

- max_cells:

  subsample ceiling for the fit

- seed:

  RNG seed; the stream is restored on exit

## Value

list with `polygon`, `hull`, `metrics`, `metrics_hull`,
`metrics_insample`, `lambda`, `W`, `b`, `markers`, or NULL if no gate
could be fitted

## Details

HOW THE TIGHTNESS IS CHOSEN: `lambda` is not a tuning parameter the
caller has to guess. The function fits at each value of `lambda_grid`,
scores every fit on held-out cells, and keeps the best; it then refines
once around the winner. Searching against held-out F1 rather than
training F1 is what stops the search from always selecting the loosest
gate, which is what fits the training cells best and generalises worst.

WHY IT RETURNS TWO POLYGONS. `polygon` is the fitted half-plane
intersection. `hull` is the convex hull of the target cells that fell
inside it – strictly tighter, usually more precise, and much more prone
to memorising the training cells, since it is drawn around them. Both
are returned with their own held-out metrics so the trade is visible
instead of being made silently.
