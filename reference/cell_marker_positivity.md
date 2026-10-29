# Is each embedded cell above its own sample's threshold, per marker

Read from the scored populations rather than from the embedding's own
columns, so every marker with a threshold is covered and not only the
ones the UMAP was built on. Scatter channels are left out: "positive"
for side scatter is a gate, not a marker.

## Usage

``` r
cell_marker_positivity(cells, pops, markers)
```

## Arguments

- cells:

  Embedded cell table with `sample_id` and `event_index`.

- pops:

  Scored populations per sample, carrying `tmat` and `thresholds`.

- markers:

  Markers to call.

## Value

Logical matrix, cells x markers.
