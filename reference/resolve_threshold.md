# Resolve a threshold for one marker, recording where the value came from

Priority (per the agreed contract):

1.  explicit config value – always wins

2.  density valley within parent – bimodal markers

3.  unstained-control quantile – unimodal markers with a control present

4.  quantile fallback – flagged needs_review A valley falling BELOW the
    control-derived value is rejected: it lies inside the background
    distribution and would call noise positive.

## Usage

``` r
resolve_threshold(
  marker,
  x_parent,
  cfg_value = NULL,
  control_x = NULL,
  control_q = 0.995,
  fallback_q = 0.9
)
```

## Arguments

- marker:

  Marker name.

- x_parent:

  The x parent.

- cfg_value:

  The cfg value.

- control_x:

  The control x.

- control_q:

  The control q. Default `0.995`.

- fallback_q:

  The fallback q. Default `0.9`.
