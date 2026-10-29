# The supplied spillover matrix as the run applied it

[`maybe_compensate()`](https://bhagesh-h.github.io/cyRAVEN/reference/maybe_compensate.md)
inverts only the detectors a file holds, so the matrix a reader must
apply is the supplied one restricted to the detectors the run's files
hold. Its columns are renamed to the files' `$PnN`, since
[`read_spillover_file()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_spillover_file.md)
also accepts a header in `make.names` spelling, which no file holds.

## Usage

``` r
gating_ml_spillover(reads, spillover)
```

## Arguments

- reads:

  the run's reads, each carrying `channel_names` and `keywords`.

- spillover:

  the matrix from
  [`read_spillover_file()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_spillover_file.md),
  or NULL.

## Value

the restricted matrix, or NULL.
