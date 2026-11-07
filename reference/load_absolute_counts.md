# Load, normalize, parse and sample-match an –absolute-counts input

Load, normalize, parse and sample-match an –absolute-counts input

## Usage

``` r
load_absolute_counts(path, smap, outdir)
```

## Arguments

- path:

  File path.

- smap:

  Sample map joining sample_id to patient_id.

- outdir:

  Directory to write outputs to.

## Value

data.frame(sample_id, sample_raw, population, cells_per_ul) or NULL
