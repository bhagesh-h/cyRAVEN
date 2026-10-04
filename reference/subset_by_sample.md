# Restrict any sample-keyed table to a set of samples

Returns NULL rather than an empty frame when nothing is left, because
every figure function in this package treats NULL as "nothing to draw"
and an empty frame as "draw an empty figure".

## Usage

``` r
subset_by_sample(x, ids)
```

## Arguments

- x:

  A data.frame carrying `sample_id`, or NULL.

- ids:

  Sample ids to keep.

## Value

The subset, or NULL.
