# Derive the cuts a specification places within parent populations

For each parent named in the specification's `cut_scope` (see
[`resolve_population_parents()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_population_parents.md)),
the parent is scored in this sample and every marker a child places
within it is cut on the parent's events alone, with the same rule the
CD45+ cuts use (`resolve`). Bounds for `bright` and `intermediate`
requirements are derived within the parent too. Parents are processed
outermost first, so a parent that is itself scoped is scored with its
own scoped cuts.

## Usage

``` r
derive_scoped_cuts(
  tmat,
  thr,
  parent,
  spec,
  hi_thr = NULL,
  resolve = NULL,
  min_events = 200L
)
```

## Arguments

- tmat:

  transformed marker matrix.

- thr:

  named CD45+ cuts.

- parent:

  the CD45+ mask.

- spec:

  resolved specification.

- hi_thr:

  CD45+ upper bounds, from
  [`derive_intermediate_bounds()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_intermediate_bounds.md).

- resolve:

  `function(marker, x)` returning a
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
  result for values `x`. Default:
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
  with its defaults.

- min_events:

  fewest parent events a cut is derived from. Default 200, the least the
  tail rule accepts.

## Value

list keyed by parent population, each
`list(thr, hi, details, n_parent)`, or NULL when the specification
scopes nothing.

## Details

A parent with fewer than `min_events` events gives no cut: its children
are reported unavailable in that sample rather than cut on a handful of
cells.
