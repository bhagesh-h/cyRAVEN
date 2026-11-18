# Spreading received by each marker from each other marker

Spreading received by each marker from each other marker

## Usage

``` r
spreading_pairs(tmat, thr, parent = NULL, min_cells = 200L)
```

## Arguments

- tmat:

  transformed marker matrix for one sample

- thr:

  named threshold vector for the same markers

- parent:

  logical mask of the cells to use, normally the parent gate

- min_cells:

  smallest group of cells that can support a spread estimate

## Value

a data.frame with one row per ordered (source, receiver) pair, or NULL
