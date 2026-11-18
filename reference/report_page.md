# Assemble a report page around its sections

Shared by the run report and the explore report, so both carry the same
search bar, numbering, contents and behaviour.

## Usage

``` r
report_page(secs, title, masthead, banner = "")
```

## Arguments

- secs:

  Sections from
  [`report_section()`](https://bhagesh-h.github.io/cyRAVEN/reference/report_section.md),
  in reading order.

- title:

  The page title.

- masthead:

  HTML for the masthead.

- banner:

  HTML placed under the masthead.

## Value

A character vector of HTML lines.
