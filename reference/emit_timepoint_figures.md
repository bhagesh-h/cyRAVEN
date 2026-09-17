# Write the figure set for one timepoint

Every call is guarded individually: a figure that cannot be drawn from
one visit's samples – a population with no cells at that visit, a marker
absent from those acquisitions – must not take the rest of the set down
with it.

## Usage

``` r
emit_timepoint_figures(ctx, dir, label)
```

## Arguments

- ctx:

  Named list of the shared artefacts: `freq`, `mfi`, `fx`, `rt`,
  `cells`, `tc`, `ufreq`, `markers`, `covariates`, `feature_cols`,
  `panel_label`, `colors`.

- dir:

  Destination directory; created if absent.

- label:

  Timepoint label, added to every figure title.

## Value

Character vector of the files written.
