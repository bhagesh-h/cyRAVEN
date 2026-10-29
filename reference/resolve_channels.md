# Decide which channels are markers and which are scatter

ONE CHANNEL PER DETECTOR: area where the detector records it, height
where it records only height, width never. Area is preferred because it
is the integrated signal, and one channel per detector because a marker
recorded as both area and height would otherwise be counted twice in the
embedding.

## Usage

``` r
resolve_channels(nm, desc, features = NULL)
```

## Arguments

- nm:

  `$PnN` values, as read.

- desc:

  `$PnS` values.

- features:

  `$PnFEATURE` values, or NULL.

## Value

list(markers = named integer, marker symbol -\> column; scatter = named
integer keyed `FSC-A`, `FSC-H`, `SSC-A`, `SSC-H`; height_markers =
markers taken from pulse height; height_only = TRUE when every marker
was)

## Details

The choice is made per detector, not per file. A file with area for one
detector and height for the others, which is how flowCore's FACSCalibur
example and many clinical archives look, used to keep the one area
channel and silently drop every height-only marker. FCS 3.2 records the
feature per parameter in `$PnFEATURE` for this reason (Spidlen et al.
2021).
