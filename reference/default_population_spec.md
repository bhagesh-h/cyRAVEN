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

ONE DEPARTURE IN HOW, NOT WHAT. The regulatory T cell population keeps
the document's name and markers, but its CD25 and CD127 cuts are derived
within the CD3+ CD4+ T cells (`parent:`) rather than across every CD45+
cell. A CD25 cut across all leukocytes falls wherever the CD25+ cells of
every lineage put it, and above it sit activated conventional T cells as
well as regulatory ones; the consensus gate is CD25+ CD127lo within CD4+
T cells (Santegoets et al. 2015; Liu et al. 2006; Seddiki et al. 2006).
