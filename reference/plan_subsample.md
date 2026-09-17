# Draw a size-balanced subsample of gated cells across samples

WHAT: returns row indices, per sample, capped so no sample dominates.
WHY: event counts differ by orders of magnitude between files. Embedding
all cells lets the largest file dictate the manifold and makes "sample"
structure indistinguishable from acquisition depth. Equal-N per sample
(or the smallest sample's N, whichever is smaller) makes cross-sample
comparison of the shared embedding meaningful.

## Usage

``` r
plan_subsample(n_per_sample, cap = 20000L, equalise = TRUE, total_cap = NULL)
```

## Arguments

- n_per_sample:

  named integer vector: sample_id -\> number of gated cells.

- cap:

  maximum cells to take from any one sample.

- equalise:

  if TRUE, take min(cap, smallest sample's N) from every sample so all
  samples contribute equally; if FALSE, take min(cap, N) per sample.

- total_cap:

  optional overall ceiling on embedded cells; the per-sample allowance
  is reduced proportionally if the total would exceed it.

## Value

named integer vector: sample_id -\> number of cells to draw.
