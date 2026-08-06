# Rank markers by how differently they are distributed in targets and the rest

WHY THE (median, p1, p99) TRIPLE AND NOT A MEAN DIFFERENCE: a marker can
separate a population by having a different SPREAD rather than a
different centre – a subset that is uniformly dim on a marker everything
else is bimodal on, for instance – and a difference of means scores that
at zero. Comparing three order statistics catches shifts in location and
in both tails, and costs one pass over the column.

## Usage

``` r
rank_gate_markers(X, y)
```

## Arguments

- X:

  n x m matrix of marker values

- y:

  0/1 target indicator

## Value

data.frame of marker and score, best first

## Details

Values are scaled to the unit interval first, so a marker on a wider
axis cannot outrank a more informative one purely by having larger
numbers.
