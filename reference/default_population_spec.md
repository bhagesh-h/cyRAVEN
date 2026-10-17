# Built-in population spec, derived from the supplied gating strategy

This is DATA, not code: each population is a list of (marker, direction)
requirements evaluated against the resolved thresholds, all within the
CD45+ parent. Users add or edit populations in the YAML config without
touching R. A population whose markers are absent from the panel is
reported UNAVAILABLE. This is a LITERAL transcription of the 15
populations in the gating strategy document, in document order. Each
entry is (marker, direction) with direction one of "above", "below" or
"intermediate". Nothing is added, renamed or inferred: populations the
document does not list (e.g. dendritic cells, double-negative T cells,
NKT cells) are deliberately absent, and gates the document specifies as
double-positive are kept double-positive.

## Usage

``` r
default_population_spec()
```

## Details

`lineage` marks the three monocyte subsets, because the document scopes
the functional marker blocks by monocyte vs non-monocyte membership.
