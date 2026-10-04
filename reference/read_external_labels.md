# Read cell labels produced by another tool

Accepts a CSV carrying one row per cell with a sample identifier, the
event's index within that sample's file, and a label. Column names are
matched loosely, so a table exported from a clustering package usually
needs no editing.

## Usage

``` r
read_external_labels(
  path,
  sample_col = NULL,
  event_col = NULL,
  label_col = NULL
)
```

## Arguments

- path:

  CSV path

- sample_col, event_col, label_col:

  column names; NULL matches by convention

## Value

data.frame(sample_id, event_index, label)

## Details

WHY THE KEY MUST BE EXPLICIT AND NOT POSITIONAL. Tools subsample.
cyCONDOR takes `max_cell` events per file and cyRAVEN takes its own cap
for the embedding, so row *i* of one table is not row *i* of the other
and a positional join silently mislabels every cell. Requiring the event
index makes the mismatch an error instead of a result.
