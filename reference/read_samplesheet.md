# Read the unified sample sheet

One row per acquired file. Splits into the three structures the pipeline
consumes, applying the same coercions as the separate loaders so that a
study analysed either way produces the same numbers.

## Usage

``` r
read_samplesheet(
  path,
  fcs_files,
  column_map = default_column_map(),
  value_map = default_value_map(),
  reference_date = Sys.Date(),
  count_unit = "cells/uL"
)
```

## Arguments

- path:

  Path to the sheet (CSV or TSV; delimiter and encoding detected).

- fcs_files:

  The fcs files being analysed, for the coverage check.

- column_map:

  The column map. Default
  [`default_column_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_column_map.md).

- value_map:

  The value map. Default
  [`default_value_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_value_map.md).

- reference_date:

  The reference date. Default
  [`Sys.Date()`](https://rdrr.io/r/base/Sys.time.html).

- count_unit:

  "cells/uL" or "cells/mL"; the latter is divided by 1000.

## Value

list(smap, patients, counts, study_columns)
