# Explore mode: MFI per cluster, and the composition of the embedding

The explore counterparts of population_mfi_summary.csv and
umap_composition.png, built by the same functions so the two modes read
alike: `explore_cluster_mfi.csv` with its heatmap and dot plot, and
`explore_cluster_composition.csv` with `explore_umap_composition.png`.

## Usage

``` r
explore_mfi_composition(
  cells,
  feats,
  prof,
  subs = NULL,
  tr = NULL,
  thr_by_sample = list(),
  scatter = character(0),
  ex_dir,
  tag = ""
)
```

## Arguments

- cells:

  Explore cell table: `sample_id`, `cluster`, `umap_1`, `umap_2`, the
  feature columns and optionally `group`.

- feats:

  Feature names.

- prof:

  Cluster profile, with `phenotype` and `frac_pos.*` columns.

- subs:

  Optional subset annotation from annotate_clusters_with_subsets().

- tr:

  The panel's transform.

- thr_by_sample:

  Named list: sample -\> named thresholds on the analysis scale.

- scatter:

  Feature names that are log10 scatter.

- ex_dir:

  Output directory.

- tag:

  Panel suffix.

## Value

Character vector of file names written.

## Details

Per-cluster medians are taken per SAMPLE first and summarised over
samples, as for declared populations, so the donor with the most events
does not set a cluster's intensity. A sample contributes to a cluster
only where it holds at least ten of its cells. Percent positive is
against each sample's own threshold where the declared run supplied one;
otherwise the "All samples" row falls back to the cluster's pooled call,
the basis explore_cluster_profile already reports.
