# Local Inverse Simpson's Index of a label over a neighbourhood graph

Local Inverse Simpson's Index of a label over a neighbourhood graph

## Usage

``` r
lisi_score(coords, labels, k = 30L)
```

## Arguments

- coords:

  numeric matrix, cells x dimensions (the embedding, or scaled marker
  space)

- labels:

  character/factor vector, one per row of `coords`

- k:

  neighbourhood size

## Value

numeric vector of per-cell iLISI
