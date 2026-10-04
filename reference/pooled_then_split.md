# A figure and its split-out variants, pooled one first

WHY THIS EXISTS RATHER THAN ONE list.files() PATTERN. Sorting cannot be
left to the locale here. R collates through the platform's collation
order, which on the usual locales ignores punctuation when comparing –
so "umap_markers_d0.png" and "umap_markers.png" compare as
"umapmarkersd0png" against "umapmarkerspng", and the day figure sorts
FIRST. The section descriptions tell the reader to take the pooled view
before the split ones, and a tab strip that opens on d0 contradicts that
on the first click.

## Usage

``` r
pooled_then_split(outdir, stem, suffix = "_[A-Za-z0-9]+")
```

## Arguments

- outdir:

  the results directory.

- stem:

  filename stem, e.g. "umap_markers".

- suffix:

  regex for the split variants' suffix, after the stem and any
  `_panel_N`. Defaults to any single alphanumeric run.

## Value

filenames present in `outdir`: pooled first, then the split ones.

## Details

The split figures are then ordered by the number inside the suffix
rather than as text, so a d10 visit follows d7 instead of landing
between d0 and d3.
