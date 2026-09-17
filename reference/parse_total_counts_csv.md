# Parse the normalized CSV into long format: patient_raw, timepoint, total_cells

Blocks are located by their header cell rather than by a fixed column
offset, because the spacer columns between blocks are a formatting
choice that varies between two exports of the same sheet.

## Usage

``` r
parse_total_counts_csv(csv_path)
```

## Arguments

- csv_path:

  The csv path.
