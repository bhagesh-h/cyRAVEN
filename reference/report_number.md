# Number the sections, and resolve cross-references between them

Called once the sections are in reading order and the empty ones are
gone, which is the first point at which a section's number is known.
Replaces the `@N@` token report_section() leaves in each section, and
resolves `@REF:<id>@` in descriptions to "section N", so prose that
points at another section cannot go stale when a section is added or a
run omits one.

## Usage

``` r
report_number(secs)
```

## Arguments

- secs:

  List of sections from
  [`report_section()`](https://bhagesh-h.github.io/cyRAVEN/reference/report_section.md).

## Value

The non-empty sections, numbered.
