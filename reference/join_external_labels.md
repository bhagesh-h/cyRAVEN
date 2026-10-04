# Attach external labels to the embedding cell table

Reports the overlap rather than assuming it. A join that matched a
handful of cells is the normal consequence of two tools having
subsampled differently, and it has to be visible before anything is
fitted to the result.

## Usage

``` r
join_external_labels(cells, labels)
```

## Arguments

- cells:

  embedding cell table with sample_id and event_index

- labels:

  the table from
  [`read_external_labels()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_external_labels.md)

## Value

`cells` with an `external_label` column, or NULL when nothing matched
