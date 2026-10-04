# Check that a subject attribute does not disagree between a patient's rows

Check that a subject attribute does not disagree between a patient's
rows

## Usage

``` r
assert_subject_consistent(df, key, cols)
```

## Arguments

- df:

  the sheet.

- key:

  character vector of patient identifiers, one per row.

- cols:

  subject columns present in the sheet.

## Value

invisible TRUE, or stops naming every conflict.
