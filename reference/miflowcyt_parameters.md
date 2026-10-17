# Per-parameter detector table for one file

`$PnN` is the detector and is always present; `$PnS` is the marker and
is optional, which is why the reading stage resolves one to the other.
Voltage, range and amplification are recorded where the instrument wrote
them.

## Usage

``` r
miflowcyt_parameters(kw)
```

## Arguments

- kw:

  keyword list for one file

## Value

a data.frame, one row per parameter
