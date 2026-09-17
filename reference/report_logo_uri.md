# The package logo as a data URI, or NA when it cannot be found

An installed package keeps `man/figures/` at `help/figures/`, while a
source tree loaded with pkgload keeps the original path, so both are
tried rather than assumed. Returns NA rather than failing: a report is
worth writing without its logo, and this runs at the end of an analysis
that may have taken an hour.

## Usage

``` r
report_logo_uri()
```
