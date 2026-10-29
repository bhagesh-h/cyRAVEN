# The per-marker embeddings, pooled first and then split by category

WHY A HELPER. A run writes one figure per marker and then one per marker
per category, so the count is a property of the cohort rather than
something a section can name in advance. Returned as paths relative to
`outdir`, which is what
[`report_section()`](https://bhagesh-h.github.io/cyRAVEN/reference/report_section.md)
and the catch-all both work in.

## Usage

``` r
marker_umap_files(outdir)
```

## Arguments

- outdir:

  The run directory.

## Details

Ordering is deliberate: the pooled view of a marker comes before its
splits, so the reader sees the whole embedding before any comparison
drawn on it.
