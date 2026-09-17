# Resolve one canonical column name against a sheet's headers

Accepts the canonical spelling and every alias in `column_map`, so a
sheet exported in German resolves without being translated first.

## Usage

``` r
match_sheet_column(nms, canon, column_map = default_column_map())
```

## Arguments

- nms:

  character vector of the sheet's column names.

- canon:

  canonical name to look for.

- column_map:

  The column map. Default
  [`default_column_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_column_map.md).

## Value

index into `nms`, or NA.
