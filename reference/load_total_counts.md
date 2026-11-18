# Load, normalize, parse and sample-match a –total-counts input

Load, normalize, parse and sample-match a –total-counts input

## Usage

``` r
load_total_counts(path, smap, outdir)
```

## Arguments

- path:

  File path.

- smap:

  Sample map joining sample_id to patient_id.

- outdir:

  Directory to write outputs to.

## Value

data.frame(sample_id, patient_raw, timepoint, total_cells), or NULL
