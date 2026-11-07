# Which pulse feature each channel records

FCS 3.2 states it in `$PnFEATURE` (Area, Height or Width). Older files
say it only through the channel name, by the `-A`/`-H`/`-W` suffix every
acquisition package writes, or through a `$PnS` ending in "- Area" as
some spectral instruments write it. The keyword wins where it is
present.

## Usage

``` r
pulse_feature(nm, desc, features = NULL)
```

## Arguments

- nm:

  `$PnN` values, as read (`alter.names` applied).

- desc:

  `$PnS` values.

- features:

  `$PnFEATURE` values, `NA` where absent. Default none.

## Value

character vector of `"A"`, `"H"`, `"W"` or `NA`.
