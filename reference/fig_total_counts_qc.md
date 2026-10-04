# Sanity figure for externally supplied total cell counts

WHY A LOG SCALE. Yields across a clinical cohort routinely span an order
of magnitude, and the failure this figure exists to catch – a value
entered in the wrong unit – lands three or six decades away. On a linear
axis one such value flattens every real one into the same pixel row and
hides exactly the problem being looked for.

## Usage

``` r
fig_total_counts_qc(
  tc,
  outfile,
  group_of = NULL,
  all_samples = NULL,
  panel_label = ""
)
```

## Arguments

- tc:

  Matched total counts from
  [`load_total_counts()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_total_counts.md).

- outfile:

  Destination PNG.

- group_of:

  Optional named vector, sample_id -\> group.

- all_samples:

  Optional character vector of every sample in the run, so the ones with
  no external total can be named rather than silently absent.

- panel_label:

  Marker-panel name added to the title; empty for none. A multi-panel
  run writes one of these per panel, and without the label in the image
  the two are distinguishable only by filename – which fails the moment
  either is pasted into a slide or a report. Same convention as
  [`fig_umap_overview()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_umap_overview.md)
  and
  [`fig_group_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_group_comparison.md).
