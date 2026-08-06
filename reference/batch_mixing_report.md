# Batch-mixing diagnostic with a permutation null and a confounding check

Batch-mixing diagnostic with a permutation null and a confounding check

## Usage

``` r
batch_mixing_report(
  cells,
  batch_col,
  group_col = "cohort",
  k = 30L,
  max_cells = 4000L,
  n_perm = 20L,
  seed = 42L
)
```

## Arguments

- cells:

  embedding cell table with umap_1/umap_2

- batch_col:

  column naming the batch (acquisition date, run, operator)

- group_col:

  biological grouping, used for the confounding check

- k:

  neighbourhood size

- max_cells:

  subsample ceiling; LISI is a local statistic and converges quickly, so
  a few thousand cells is ample and keeps this from dominating the
  runtime of the whole pipeline

- n_perm:

  permutations for the null

- seed:

  Random seed. The RNG stream is restored afterwards, so this cannot
  perturb later sampling. Default `42L`.

## Value

list(summary = data.frame, per_cell = numeric, confounding = data.frame)
