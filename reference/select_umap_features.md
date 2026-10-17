# Choose the markers that define the embedding space

WHAT: returns the subset of available markers to use as UMAP input
features. WHY: the embedding must be driven by lineage/phenotype
identity, not by (a) scatter channels, which encode size/granularity and
are already used for gating, (b) the viability dye, whose variation is
technical, (c) height/width channels, which are redundant duplicates of
area and would double-weight every marker, or (d) time. The original
template selected features by channel-name regex and so retained 26
redundant height channels – the single largest cause of its
uninformative UMAP.

## Usage

``` r
select_umap_features(
  marker_cols,
  prefer = NULL,
  exclude = NULL,
  lineage_only = TRUE
)
```

## Arguments

- marker_cols:

  named integer vector: marker symbol -\> column index (AREA channels
  only, as produced by the reading module).

- prefer:

  character vector of preferred lineage markers; those present are used.
  If none of them are present, falls back to all eligible markers.

- exclude:

  character vector of marker symbols to drop unconditionally (e.g. the
  detected viability marker).

- lineage_only:

  logical; if FALSE, every eligible marker is used.

## Value

character vector of marker symbols, in stable sorted order.
