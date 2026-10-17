# Does the gate work on a donor it was not fitted to

Refits the whole strategy with one donor withheld and scores it on that
donor, once per donor. Returns the per-donor scores and their spread.

## Usage

``` r
gate_transferability(
  X,
  y,
  donor,
  min_donors = 3L,
  max_donors = 8L,
  max_cells = 20000L,
  seed = 42L,
  ...
)
```

## Arguments

- X:

  marker matrix

- y:

  0/1 target indicator

- donor:

  grouping vector, one entry per row of `X`

- min_donors:

  refuse below this many donors, where the statistic would be an
  anecdote

- max_donors:

  ceiling on folds

- max_cells:

  ceiling on training rows per fold

- seed:

  seed for the local RNG stream

- ...:

  passed to
  [`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)

## Value

list(per_donor = data.frame, summary = data.frame), or NULL

## Details

READ THE MINIMUM, NOT THE MEAN. A gate that transfers is one whose worst
donor is acceptable. The mean hides exactly the failure this function
exists to expose.

COST, AND THE TWO CAPS THAT BOUND IT. Every fold refits the whole
strategy, so the work is one fit per donor per label and grows with
both. On a cohort of twenty-odd donors and half a dozen labels an
uncapped version runs for hours, which is not a validation statistic
anyone will wait for.

`max_donors` bounds the number of folds. The subset is drawn at random
from the fixed seed rather than taken as the largest donors, because a
worst-donor statistic computed only on the best-represented donors is
optimistic in exactly the direction that matters. `max_cells` subsamples
each fold's TRAINING rows, stratified by label; the held-out donor is
always scored in full. Both caps are logged when they bind, so a summary
over eight donors is never mistaken for one over twenty.
