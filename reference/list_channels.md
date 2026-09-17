# List the acquisition parameters and the names the run resolves them to

One row per parameter: index, the detector name (`$PnN`), the
description (`$PnS`), the symbol the run resolves it to, and how the run
will use it. The resolved symbol is the string a population
specification has to name.

## Usage

``` r
list_channels(fcs)
```

## Arguments

- fcs:

  Paths to the FCS files.

## Value

invisible data.frame of the reference file's parameters.

## Details

Reads headers only, so it costs nothing on large files, and reads every
file rather than the first so a panel that differs between files is
reported.
