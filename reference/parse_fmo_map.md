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
carrying `sample_id`, `marker`, `control_group`, and `control_kind`
(`minus_one` where the file leaves out a single reagent,
`minus_multiple` where it leaves out several) with `n_markers_in_file`;
NULL when no FMO is declared
