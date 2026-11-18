# Legacy HTML table renderer

Retained because the interactive tables are rendered from JSON in the
browser; this is the static fallback used when a table is too small for
searching to be worth the controls.

## Usage

``` r
html_table(d, max_rows = 20L)
```

## Arguments

- d:

  a data.frame

- max_rows:

  rows beyond which the table is truncated, with a note
