# Correct batch effects in the transformed marker matrix

WHAT IT DOES: aligns each marker's distribution across batches, after
checking that batch and biological group are separable enough for that
to be meaningful.

## Usage

``` r
correct_batch(
  tmat,
  batch,
  group = NULL,
  markers = NULL,
  max_cramers_v = 0.6,
  force = FALSE,
  method = c("quantile", "cluster", "cytonorm"),
  cluster_k = 10L,
  seed = 42L
)
```

## Arguments

- tmat:

  numeric matrix, cells x markers, already transformed

- batch:

  batch label per cell

- group:

  biological group per cell, for the confounding check

- markers:

  markers to correct; defaults to all columns

- max_cramers_v:

  refusal threshold

- force:

  proceed despite confounding

- method:

  `"quantile"` for one map per marker over the whole file, or
  `"cluster"` (synonym `"cytonorm"`) for one map per marker per cell
  type

- cluster_k:

  number of cell-type clusters when `method = "cluster"`

- seed:

  RNG seed for the clustering, so a run reproduces

## Value

list(tmat, cramers_v, corrected, markers, reason, method)

## Details

WHAT IT REFUSES TO DO: correct when Cramer's V between batch and group
exceeds `max_cramers_v`. At that point the two are close to the same
variable, and any correction removes the effect being looked for.
Passing `force = TRUE` proceeds and records the decision in the returned
object so it reaches the run manifest.

This is descriptive of intensity only. It does not touch the gate
hierarchy, which is derived per sample and is already batch-local by
construction.

TWO METHODS. `"quantile"` fits one map per marker over the whole file.
`"cluster"` fits one map per marker per cell type, which is what
CytoNorm does and what a whole-file map cannot do. `"cytonorm"` is
accepted as a synonym for `"cluster"`, because that is the name the
method is known by.

THE METHOD IS CHOSEN AFTER THE REFUSAL, NEVER BEFORE IT. A better
alignment algorithm does not make a confounded design correctable. Where
Cramer's V says batch and group are close to the same variable, both
methods are refused identically, and `force = TRUE` remains the only way
past.
