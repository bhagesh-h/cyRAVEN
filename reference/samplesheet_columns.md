# Reserved column names of the unified sample sheet

`acquisition` names the file and how it was acquired; `subject` names
properties of the patient, which repeat across that patient's rows. Any
column not listed here and not carrying the count prefix is kept as a
study variable, usable as `--group-column`, `--batch-column` or a
covariate.

## Usage

``` r
samplesheet_columns()
```

## Value

list(acquisition, subject, required, count_prefix)
