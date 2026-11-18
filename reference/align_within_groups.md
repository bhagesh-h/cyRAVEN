# Align each batch towards a goal that respects the design

WHY NOT ONE POOLED GOAL. Mapping each batch's distribution onto the
pooled distribution of every sample assumes the batches hold the same
mixture of biology. Where one batch holds mostly patients and another
mostly controls, the pooled goal is a blend of the two groups, each
batch is pulled towards it, and a real difference between the groups is
reduced along with the batch effect. That homogeneity is exactly what
anchor samples exist to avoid assuming (Schuyler et al. 2019). Without
anchors, CytoNorm 2.0 tailors its goal distribution to the experimental
design (Quintelier et al. 2025).

## Usage

``` r
align_within_groups(x, batch, group = NULL, goal = "batch_mean")
```

## Arguments

- x:

  numeric vector for one marker.

- batch:

  batch label per cell.

- group:

  study group per cell, or NULL for one goal over all cells.

- goal:

  passed to
  [`align_quantiles()`](https://bhagesh-h.github.io/cyRAVEN/reference/align_quantiles.md).

## Value

numeric vector, corrected.

## Details

WHAT IS DONE. Each study group is aligned separately: a batch's cells of
group g are mapped onto group g's goal, the mean of its quantiles over
the batches EVERY group was measured in (see
[`common_batches()`](https://bhagesh-h.github.io/cyRAVEN/reference/common_batches.md)
and
[`align_quantiles()`](https://bhagesh-h.github.io/cyRAVEN/reference/align_quantiles.md)).
A batch effect common to the groups is removed in each of them, and the
groups keep whatever separates them. Cells with no group are aligned
towards the same batches.

WHY THE GOAL USES ONLY THE SHARED BATCHES. Each group's goal carries the
mean batch effect of the batches it is built from. Built from each
group's own batches, two groups measured in different batches end on
different scales: with batch 2 reading one unit high, a group measured
in batch 1 alone stayed where it was while a group spread over both
moved half a unit up, so a correction meant to remove a batch effect put
a 0.5 difference between groups of identical cells. Built from the same
batches for every group, every goal carries the same batch effect. Where
no batch holds every group there is no shared scale to align to, and
nothing is moved.
