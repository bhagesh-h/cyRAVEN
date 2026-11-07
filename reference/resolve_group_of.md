# Resolve sample_id -\> an arbitrary patient-table column, via the sample map

Shared by the between-group figures (grouping column, e.g. "cohort") and
by –absolute-counts sample matching (see below): both need patient_id as
the join key between the patient table and the sample map's per-file
rows.

## Usage

``` r
resolve_group_of(patients, smap, gcol)
```

## Arguments

- patients:

  Patient metadata table, keyed by patient_id.

- smap:

  Sample map joining sample_id to patient_id.

- gcol:

  The gcol.
