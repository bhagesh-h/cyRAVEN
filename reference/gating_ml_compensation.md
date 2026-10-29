# The compensation a Gating-ML gate should declare for this run

"supplied" where the files were compensated with the `--spillover`
matrix, which the document must then carry
([`gating_ml_spillover()`](https://bhagesh-h.github.io/cyRAVEN/reference/gating_ml_spillover.md));
"FCS" where they were compensated with their own `$SPILLOVER`;
"uncompensated" where none had a matrix and the stored values were gated
as they are.

## Usage

``` r
gating_ml_compensation(reads)
```

## Arguments

- reads:

  the run's reads, each carrying `compensation`.

## Value

"supplied", "FCS" or "uncompensated".
