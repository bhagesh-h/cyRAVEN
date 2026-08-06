# Precision, recall and F1 for a binary gate

Recall is measured against `n_target_total` rather than against the
targets present, so that in a hierarchy it stays CUMULATIVE: a level-3
gate is judged on the fraction of the ORIGINAL population it still
holds, not on the fraction of whatever survived level 2. A hierarchy
that discards half the population at each level and reports 100% recall
three times is exactly the failure this guards against.

## Usage

``` r
gate_metrics(inside, y, n_target_total = sum(y == 1))
```

## Arguments

- inside:

  logical, gate membership

- y:

  0/1 target indicator

- n_target_total:

  denominator for recall

## Value

named numeric vector
