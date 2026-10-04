# Locate the Time channel in a read FCS file

The reading stage keeps Time out of `marker_cols` and `scatter_cols` on
purpose, so it cannot drive a gate or an embedding. This finds it in the
raw matrix for the one use it does have.

## Usage

``` r
find_time_column(rd)
```

## Arguments

- rd:

  a list from
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)

## Value

the column index, or NA when the file carries no Time channel
