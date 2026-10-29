# Parse a whole column of dates, resolving day/month order from the column

WHY NOT try formats one at a time until one sticks: R's `%Y` happily
accepts a two-digit year, so `as.Date("12/11/58", "%d/%m/%Y")` returns
0058-11-12 rather than failing. A trial-and-error loop therefore locks
onto the first permissive format and silently produces dates ~1900 years
wrong. Two-digit years must be parsed as such, and the field order
decided from the column as a whole.

## Usage

``` r
parse_dates_column(v, reference_date = Sys.Date(), label = "date")
```

## Arguments

- v:

  character vector of raw date strings.

- reference_date:

  Date used to place two-digit years and to flag values in the future.

- label:

  column name, used in messages.

## Value

Date vector, NA where unparseable.

## Details

WHY column-wise: "12/11/58" is ambiguous in isolation (12 Nov or 11
Dec?), but a column of clinical dates almost always contains at least
one value whose first component exceeds 12, which settles the order for
every row. Deciding per-row would silently mix conventions within one
column.
