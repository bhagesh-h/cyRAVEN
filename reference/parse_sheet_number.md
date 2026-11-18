# Read a number written in a spreadsheet's own locale

Spreadsheets export what a cell displays, thousands separators included,
and a count in cells/uL is routinely above 1,000. Replacing every comma
by a point read "1,500" as 1.5, and "1.500" is 1.5 to as.numeric(): a
thousandfold error either way. The delimiter tells the conventions
apart, as the export that wrote the file followed one: a
semicolon-separated file is written with decimal commas and points
between thousands, a comma-separated one with decimal points and commas
between thousands. Only a separator in thousands position (groups of
exactly three digits) is removed.

## Usage

``` r
parse_sheet_number(x, delim = ",")
```

## Arguments

- x:

  values as read.

- delim:

  the sheet's delimiter.

## Value

numeric vector.
