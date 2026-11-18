# A CD45-negative reference from the events the scatter gate rejected

Events well below the scatter gate's FSC boundary are sub-cellular
debris and carry no CD45, so the upper tail of their CD45 values is
where CD45-negative ends in this tube, measured under the same staining
and the same settings.

## Usage

``` r
cd45_negative_reference(cd45_x, in_gate, q = 0.995, min_events = 500L)
```

## Arguments

- cd45_x:

  Transformed CD45 values for every event.

- in_gate:

  Logical, `TRUE` for events NOT to use as debris. The caller passes
  everything except the debris core, the events well below the FSC
  boundary.

- q:

  Quantile taken. Default `0.995`.

- min_events:

  Fewer rejected events than this gives `NA`. Default `500L`.

## Value

The reference value, or `NA_real_`.
