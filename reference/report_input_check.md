# Report what a run would do, from headers and metadata alone

Report what a run would do, from headers and metadata alone

## Usage

``` r
report_input_check(fcs, smap, sheet, spec, opt)
```

## Arguments

- fcs:

  Paths to the FCS files.

- smap:

  Sample map or sheet-derived equivalent; may be NULL.

- sheet:

  Result of
  [`read_samplesheet()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_samplesheet.md),
  or NULL on the three-file route.

- spec:

  Population specification.

- opt:

  Parsed options.

## Value

invisible TRUE when nothing fatal was found, FALSE otherwise.
