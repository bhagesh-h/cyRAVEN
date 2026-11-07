# Percent of each label positive for each marker, from per-cell calls

Percent of each label positive for each marker, from per-cell calls

## Usage

``` r
composition_marker_pct(labels, pos, w = NULL)
```

## Arguments

- labels:

  Character vector, one label per cell.

- pos:

  Logical matrix, cells x markers: is the cell above its own sample's
  threshold for the marker. NA where no threshold exists.

- w:

  Optional per-cell weights (`sampling_weight`).

## Value

Numeric matrix, labels x markers, in percent; NA where a label has no
cell with a call for the marker.
