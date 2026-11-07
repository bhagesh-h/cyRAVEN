# Parse the normalized CSV into long format: sample_raw, population, cells_per_ul

Auto-detects which columns are real population data (header text AND at
least one numeric value somewhere below it) versus a units-label spacer
column (header text, never any numeric value below it – e.g. "in 1 mL
blood"), and uses that spacer text, if found, to convert per-mL counts
to the cells_per_ul convention this pipeline's own wbc_per_ul route
uses. Falls back to assuming cells/uL already, LOUDLY, when no unit is
found – silently guessing a 1000x factor wrong is worse than asking the
user to check.

## Usage

``` r
parse_absolute_counts_csv(csv_path)
```

## Arguments

- csv_path:

  The csv path.
