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
  force = FALSE
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

## Value

list(tmat, cramers_v, corrected, markers, reason)

## Details

WHAT IT REFUSES TO DO: correct when Cramer's V between batch and group
exceeds `max_cramers_v`. At that point the two are close to the same
variable, and any correction removes the effect being looked for.
Passing `force = TRUE` proceeds and records the decision in the returned
object so it reaches the run manifest.

This is descriptive of intensity only. It does not touch the gate
hierarchy, which is derived per sample and is already batch-local by
construction.
