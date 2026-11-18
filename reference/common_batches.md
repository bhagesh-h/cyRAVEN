# Batches holding cells of every study group

Batches holding cells of every study group

## Usage

``` r
common_batches(batch, group, min_cells = 10L)
```

## Arguments

- batch:

  batch label per cell.

- group:

  study group per cell.

- min_cells:

  cells of a group a batch must hold to count, the size
  [`align_quantiles()`](https://bhagesh-h.github.io/cyRAVEN/reference/align_quantiles.md)
  needs to build a batch's quantiles.

## Value

character vector of batch labels, possibly empty.
