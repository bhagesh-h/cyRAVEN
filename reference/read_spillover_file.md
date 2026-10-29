# Read a spillover matrix from a CSV file

The layout FlowJo and AutoSpill export and the `$SPILLOVER` keyword
encodes: a square matrix whose columns are the detectors, row i giving
the fraction of detector i's fluorochrome seen in each detector, with
ones on the diagonal. A first column holding the row names is accepted
and dropped.

## Usage

``` r
read_spillover_file(path)
```

## Arguments

- path:

  CSV path, or NULL.

## Value

numeric matrix with detector column names, or NULL when `path` is NULL
or empty.
