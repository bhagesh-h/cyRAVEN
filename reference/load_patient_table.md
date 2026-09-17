# Read, filter, rename and translate the patient table

Detects the delimiter and the encoding (German exports are often
latin1), trims whitespace, tolerates umlaut spelling variants, converts
decimal commas, and PASSES THROUGH any value not in the dictionary while
reporting it – silently mangling an unmapped clinical value would be
worse than leaving it in German for the user to map.

## Usage

``` r
load_patient_table(
  path,
  column_map = default_column_map(),
  value_map = default_value_map(),
  reference_date = Sys.Date()
)
```

## Arguments

- path:

  File path.

- column_map:

  The column map. Default
  [`default_column_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_column_map.md).

- value_map:

  The value map. Default
  [`default_value_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_value_map.md).

- reference_date:

  The reference date. Default
  [`Sys.Date()`](https://rdrr.io/r/base/Sys.time.html).
