# Apply the spillover matrix if one is present

WHY: spectral instruments (e.g. Sony ID7000) write already-unmixed data
and carry no \$SPILLOVER keyword. Compensating twice, or erroring
because the matrix is absent, are both wrong – detect and report.

## Usage

``` r
maybe_compensate(ff_exprs, keywords, external = NULL)
```

## Arguments

- ff_exprs:

  The ff exprs.

- keywords:

  FCS keyword list.

- external:

  A spillover matrix from
  [`read_spillover_file()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_spillover_file.md),
  or NULL.

## Value

`ff_exprs`, compensated where a matrix applied, with attribute
`compensation` stating what was done.

## Details

WHAT AN ABSENT KEYWORD MEANS. Nothing, on its own. `$SPILLOVER` is the
only standard carrier of a compensation matrix (FCS 3.1, Spidlen et al.
2010), so a file without it is either unmixed or compensated already, or
conventional data nobody compensated. The keyword cannot tell those
apart, and the result says the status is unknown rather than guessing
the reassuring answer. FCS 3.2 files can say it per parameter in
`$PnTYPE` ("Unmixed_Fluorescence"), and that is read where it is
written.

A MATRIX FROM ELSEWHERE. `external`, from `--spillover`, replaces the
file's own matrix: a matrix recomputed from single-stain controls, in
AutoSpill or in the acquisition software, is applied to every file
instead of whatever the files carry.
