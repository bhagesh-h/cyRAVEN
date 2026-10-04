# Record one gate adjustment for gate_adjustments.csv

Record one gate adjustment for gate_adjustments.csv

## Usage

``` r
gate_adjustment_row(
  sample_id,
  gate,
  marker,
  source,
  threshold,
  action,
  reason,
  n_parent,
  n_kept
)
```

## Arguments

- sample_id:

  Sample identifier.

- gate:

  Gate name.

- marker:

  Marker the gate is placed on.

- source:

  Threshold source.

- threshold:

  The threshold that would have been applied.

- action:

  "skipped" or "kept".

- reason:

  Free text.

- n_parent:

  Events entering the gate.

- n_kept:

  Events the gate retained.
