# Resolve a threshold for one marker, recording where the value came from

Priority (per the agreed contract): 0. per-sample manual override – wins
over everything, and is recorded

1.  explicit config value – applies to every sample

2.  density valley within parent – bimodal markers

3.  control quantile – unimodal markers with a control present

4.  tail gate – unimodal markers without a control; flagged needs_review

5.  quantile fallback – only where the tail gate cannot be placed;
    flagged needs_review A valley falling BELOW the control-derived
    value is rejected: it lies inside the background distribution and
    would call noise positive.

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
  adaptive = FALSE,
  floor = -Inf,
  prefer = c("deepest", "lowest"),
  fmo_sets_cut = FALSE,
  tail = c("upper", "lower")
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

- floor:

  No derived cut is accepted at or below this value; see
  [`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md).
  Pass `transform$fn(0, marker)` for a fluorescence marker. Default
  `-Inf`. A manual, config or control cut is not floored: each is a
  decision someone made, not one the density made.

- prefer:

  Passed to
  [`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md)
  and
  [`best_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/best_threshold.md):
  `"lowest"` for a gate that should remove only the negative mode, such
  as CD45. Default `"deepest"`.

- fmo_sets_cut:

  When `TRUE` and the control is a true FMO
  (`control_kind = "fmo_q995"`), the FMO's cut is used even where the
  sample's own valley sits higher, and the source reads
  `fmo_q995_valley_lowered`. The FMO is the published reference for
  where positivity starts (Maecker & Trotter 2006; Cossarizza et al.
  2019), so events between its edge and a higher valley are dim
  positives, not negatives. `FALSE` by default: an FMO argues a cut down
  validly only where its negative is the sample's own, the same cells
  stained the same day, and one FMO tube serves every sample of its
  panel or `control_group`, which this package cannot check against each
  sample; a valley is evidence from the sample itself. Raising a valley
  that sits below the FMO's edge is applied by default, since that errs
  towards fewer false positives. Never applied to an unstained or
  minus-multiple tube: their negatives lack the panel's spreading and
  sit too low to argue a cut down.

- tail:

  which side of the dominant mode a tail gate goes, where one is needed:
  `"upper"` (default) takes the mode to be the negative population;
  `"lower"` is the analyst's declaration that it is the positive one,
  read from `thresholds: <marker>: {tail: lower}`. See
  [`tail_gate_cut()`](https://bhagesh-h.github.io/cyRAVEN/reference/tail_gate_cut.md).

## Details

WHY THE TAIL GATE COMES BEFORE THE QUANTILE. A marker with one mode and
no control used to be cut at its 90th percentile, which reports exactly
10% positive whatever the sample holds: a constant, not a measurement,
and no source supports it. The tail gate places the cut a fixed number
of robust standard deviations above the dominant mode, so a sample with
no positive cells reports almost none and a sample with a bright tail
reports the tail. It answers the case openCyto's `gate_tail` exists for,
a single mode with a positive tail (Finak et al. 2014), by a different
rule: see
[`tail_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/tail_threshold.md).
It assumes the dominant mode is the negative one, which a single mode
cannot confirm (see
[`tail_gate_cut()`](https://bhagesh-h.github.io/cyRAVEN/reference/tail_gate_cut.md)
for why the side is not guessed), so the cut stays flagged for review.
The quantile is kept only for the case the tail rule cannot handle, a
distribution with no spread to measure.

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
