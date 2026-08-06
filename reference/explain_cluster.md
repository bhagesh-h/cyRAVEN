# Learn a hierarchical gating strategy for one labelled population

WHAT: repeatedly picks the two most discriminating markers among the
cells that survive so far, fits a convex gate in that plane, keeps the
cells inside it, and goes again – which is exactly the shape of a manual
gating strategy, and is why the output can be executed by hand.

## Usage

``` r
explain_cluster(X, y, max_depth = 4L, min_gain = 0.02, markers = NULL, ...)
```

## Arguments

- X:

  n x m matrix of marker values (arcsinh scale)

- y:

  0/1 or logical target indicator

- max_depth:

  maximum gates in the strategy

- min_gain:

  minimum held-out F1 improvement to justify another level

- markers:

  optional restriction of the marker set

- ...:

  passed to
  [`learn_convex_gate()`](https://bhagesh-h.github.io/cyRAVEN/reference/learn_convex_gate.md)

## Value

list(levels = list of gates, summary = data.frame, best_depth), or NULL

## Details

WHY IT STOPS EARLY. Each additional level can only remove cells, so
recall falls monotonically while precision usually rises. The strategy
is worth extending only while F1 improves; `min_gain` sets how much
improvement counts as worth another gate for whoever has to draw it. It
also stops when a level fails to fit, rather than skipping it – a gap in
the middle of a hierarchy is not a gating strategy, and returning the
levels that did work is the honest truncation.
