# Choose k per population from the reference group's cells

Choose k per population from the reference group's cells

## Usage

``` r
choose_subcluster_k(
  cells,
  markers,
  group_col = "cohort",
  reference = NULL,
  k_range = 2:5,
  min_ref = 150L,
  seed = 42L
)
```

## Arguments

- cells:

  Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`,
  `sample_id`, `population_label` and one column per marker.

- markers:

  Character vector of marker names to use.

- group_col:

  Name of the column holding the biological grouping (cohort). Default
  `"cohort"`.

- reference:

  The group every other group is compared against.

- k_range:

  Candidate values of k to score. Default `2:5`.

- min_ref:

  Minimum reference-group cells required to fit. Default `150L`.

- seed:

  Random seed. The RNG stream is restored afterwards, so this cannot
  perturb later sampling. Default `42L`.

## Value

list(k = named integer population -\> k, curve = data.frame of scores)
