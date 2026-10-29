# The marker name cyRAVEN resolves from an FCS description

Spectral instruments write `$PnS` as `"CD45 : SparkUV-387 - Area"`, the
antibody, the detector, and the pulse statistic in one string. The name
a population specification has to use is the antibody alone, so
everything from `" : "` onward is stripped. Where `$PnS` is empty the
channel name `$PnN` stands in.

## Usage

``` r
marker_symbol(desc, nm = NULL)
```

## Arguments

- desc:

  `$PnS` values.

- nm:

  `$PnN` values, used where `desc` is empty.

## Value

Character vector of resolved marker names.

## Details

EXTRACTED SO IT CANNOT DIVERGE.
[`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)
and
[`report_input_check()`](https://bhagesh-h.github.io/cyRAVEN/reference/report_input_check.md)
both need it, and for a while only the first had it: `--check` then
reported markers as the full `$PnS` string and told the user that every
population named a marker no file contained, on data where the run
itself would have resolved all of them. A false alarm from the one
command whose whole purpose is to catch real ones.
