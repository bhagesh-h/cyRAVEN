# Align one marker's distribution across batches

WHAT: for each batch, build the monotone map that carries that batch's
quantiles onto the pooled reference quantiles, and apply it.

## Usage

``` r
align_quantiles(x, batch, probs = seq(0, 1, length.out = 101))
```

## Arguments

- x:

  numeric vector for one marker, across all cells

- batch:

  batch label per cell

- probs:

  quantile grid used to build the map

## Value

numeric vector, corrected

## Details

WHY MONOTONE INTERPOLATION AND NOT A SHIFT OR A Z-SCORE: a shift
corrects the location of a distribution and leaves its shape, so a batch
whose negative population is wider stays wider and its gate still lands
somewhere else. A z-score assumes the spread is meaningful on both sides
of zero, which is false for a bimodal marker where the two modes have
different variances. Matching quantiles corrects the shape too, and
because the map is monotone it preserves the ORDER of cells within a
batch – the property that stops a correction from inventing structure.
