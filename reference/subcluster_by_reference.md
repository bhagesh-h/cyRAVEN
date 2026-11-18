# Split each cluster into subclusters DEFINED ON THE REFERENCE GROUP

WHY THE REFERENCE DEFINES THE SUBCLUSTERS, AND NOT ALL CELLS POOLED: the
question is "what inside this cluster differs from healthy". Clustering
the pooled cells would let the patient cells help draw the very
boundaries they are then tested against – a disease-specific subset
would carve out its own subcluster and come back looking like a normal
compartment that simply contains patient cells. Fitting k-means on the
HEALTHY cells alone makes each subcluster a piece of normal biology, and
assigning every study's cells to the nearest healthy centroid asks the
honest question: given the compartments a healthy immune system has,
where do patient cells fall, and do they look the same once they get
there?

## Usage

``` r
subcluster_by_reference(
  cells,
  markers,
  group_col = "cohort",
  reference = NULL,
  k = 3L,
  min_ref = 150L,
  min_cell = 20L,
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

  the group whose cells define the subclusters (e.g. controls)

- k:

  maximum subclusters per cluster; reduced when the reference is small

- min_ref:

  Minimum reference-group cells required to fit. Default `150L`.

- min_cell:

  Minimum cells a population must have before it is split. Default
  `20L`.

- seed:

  Random seed. The RNG stream is restored afterwards, so this cannot
  perturb later sampling. Default `42L`.

## Value

integer vector, one subcluster id per row of `cells`

## Details

That yields the two ways a study can deviate, and the figure and table
downstream report both: OCCUPANCY a subcluster holds a different share
of the cluster's cells INTENSITY the cells that are there express a
marker differently

A cluster whose reference has too few cells to subcluster safely is left
whole (subcluster 1), which is the honest answer rather than splitting
noise.
