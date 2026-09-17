# Apply the spillover matrix if one is present WHY: spectral instruments (e.g. Sony ID7000) write already-unmixed data and carry no \$SPILLOVER keyword. Compensating twice, or erroring because the matrix is absent, are both wrong – detect and report.

Apply the spillover matrix if one is present WHY: spectral instruments
(e.g. Sony ID7000) write already-unmixed data and carry no \$SPILLOVER
keyword. Compensating twice, or erroring because the matrix is absent,
are both wrong – detect and report.

## Usage

``` r
maybe_compensate(ff_exprs, keywords)
```

## Arguments

- ff_exprs:

  The ff exprs.

- keywords:

  FCS keyword list.
