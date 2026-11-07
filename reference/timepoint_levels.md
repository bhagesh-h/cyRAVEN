# The timepoint levels present in a sample sheet, in collection order

The timepoint levels present in a sample sheet, in collection order

## Usage

``` r
timepoint_levels(smap, min_samples = 2L)
```

## Arguments

- smap:

  Sample map.

- min_samples:

  Levels with fewer samples than this are dropped: a figure set for a
  single acquisition is a set of one-column boxplots and one-sample
  UMAPs, which is noise rather than a view.

## Value

Character vector, or `NULL` when there is nothing to split by.
