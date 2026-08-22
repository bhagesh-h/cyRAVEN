# Parse the `fmo_for` column into a marker-to-file map

Parse the `fmo_for` column into a marker-to-file map

## Usage

``` r
parse_fmo_map(smap)
```

## Arguments

- smap:

  the sample map, or NULL

## Value

a data.frame with one row per (file, marker) the file controls for,
carrying `sample_id`, `marker` and `control_group`; NULL when no FMO is
declared
