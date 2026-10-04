# Size above which a table is named rather than embedded

Set with `options(cyRAVEN.report_table_max_mb = )`. The default of 8 MB
sits above every summary table a run writes and below the per-cell
exports, whose row count is the number of cells rather than the number
of samples.

## Usage

``` r
report_table_max_bytes()
```
