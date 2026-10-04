# Resolve a threshold for one marker, recording where the value came from

Priority (per the agreed contract): 0. per-sample manual override – wins
over everything, and is recorded

1.  explicit config value – applies to every sample

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
  fallback_q = 0.9,
  override = NULL,
  control_kind = c("control_q995", "fmo_q995", "fmm_q995"),
  control_threshold = NULL,
  adaptive = FALSE
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

- override:

  Optional list for THIS sample and marker carrying `threshold`, and
  optionally `reason` and `set_by`. See
  [`sample_override()`](https://bhagesh-h.github.io/cyRAVEN/reference/sample_override.md).

- control_kind:

  What `control_x` is, which decides the `source` string recorded.
  `"control_q995"` for an unstained tube, `"fmo_q995"` for a
  fluorescence-minus-one control, `"fmm_q995"` for a tube missing
  several reagents at once, whose negative is narrower than a true
  FMO's, so the cut it anchors is permissive. The arithmetic is
  identical in all three cases; they are named apart because they are
  different experiments and support different claims. See
  [`parse_fmo_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/parse_fmo_map.md).

- control_threshold:

  A cut already placed on the control by the caller, used in place of
  `control_q` of `control_x`.
  [`control_cut()`](https://bhagesh-h.github.io/cyRAVEN/reference/control_cut.md)
  produces one and refuses a control that is brighter than the sample,
  which a minus control cannot be. `NULL` by default, in which case the
  quantile is taken here as before.

- adaptive:

  Sweep the kernel bandwidth and try the tail and Otsu rules before
  falling back to a quantile. See
  [`best_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/best_threshold.md).
  `FALSE` by default, because it moves thresholds and a threshold that
  moves moves every frequency beneath it.

## Details

WHY LEVEL 0 EXISTS AND WHY IT IS DISTINCT FROM LEVEL 1. When this
package flags a threshold for review, the only previous responses were
to accept it or to pin that marker across the whole run through
`thresholds:`. Pinning is the worse of the two: it applies one number to
every sample and so reintroduces exactly the fixed-coordinate bias
described in README section 1, in order to correct one tube. A
per-sample override corrects the tube.

It is a separate `source` value rather than reusing `config` because the
two are different claims. `config` says the assay declares this cut;
`manual` says a named person moved this one sample's cut for a stated
reason on a stated run. An override that is recorded and attributed is
more defensible than an automated cut nobody was permitted to correct,
and less defensible than one nobody needed to.
