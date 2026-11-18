# Count-model differential abundance, as diffcyt-DA-edgeR

Count-model differential abundance, as diffcyt-DA-edgeR

## Usage

``` r
stats_edger_abundance(
  counts,
  group_of,
  reference = NULL,
  min_cells = 3L,
  min_samples = NULL,
  min_n = 3L
)
```

## Arguments

- counts:

  long table with `sample_id`, `population`, `count` (exclusive labels),
  and optionally `is_control`/`qc_status`.

- group_of:

  named vector, sample id to group.

- reference:

  reference group, or NULL for the first.

- min_cells, min_samples:

  a population is tested when it holds at least `min_cells` cells in at
  least `min_samples` samples (default half of them).

- min_n:

  fewest samples per group for a group to be tested.

## Value

data.frame, one row per population and comparison, or NULL.
