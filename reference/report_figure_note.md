# One line saying what a figure shows

The figure counterpart of
[`report_table_note()`](https://bhagesh-h.github.io/cyRAVEN/reference/report_table_note.md).
It is printed under each figure and is what the report's search bar
matches when a reader types what they are looking for ("composition",
"MFI", "threshold") rather than a file name they have never seen.
Figures whose names carry a marker, a variable or a timepoint are
matched by pattern, and a panel or timepoint suffix is stripped before
the exact lookup. An unknown figure gets no line.

## Usage

``` r
report_figure_note(file)
```

## Arguments

- file:

  figure file name, with or without a directory.

## Value

A single string, or "" when the figure is not known.
