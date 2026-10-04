# One `$Pn*` keyword, or "" when the file does not carry it

read.FCSHeader returns a NAMED CHARACTER VECTOR rather than a list, and
`[[` on one raises "subscript out of bounds" for a name that is absent,
where the same call on a list returns NULL. An `if (is.null(v))` guard
therefore never runs: the error is raised before it is reached.

## Usage

``` r
kw_value(kw, i, suffix)
```

## Arguments

- kw:

  Keyword vector from
  [`flowCore::read.FCSheader()`](https://rdrr.io/pkg/flowCore/man/read.FCSheader.html).

- i:

  Parameter index.

- suffix:

  "N" or "S".

## Value

Trimmed character scalar, "" when absent or NA.

## Details

A parameter with no `$PnS` is ordinary rather than exceptional. Time
channels routinely carry none, and some instruments omit it for scatter.
Reading a cohort containing one used to abort the whole command.
