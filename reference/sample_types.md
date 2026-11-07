# Each sample's specimen type, where the sheet records one

Read from a `sample_type`, `specimen`, `material` or `sample_material`
column of the sample map, or of the patient table through `patient_id`.

## Usage

``` r
sample_types(smap, patients = NULL)
```

## Arguments

- smap:

  sample map, or NULL.

- patients:

  patient table, or NULL.

## Value

named character vector, sample id -\> type, or NULL when no such column
exists.
