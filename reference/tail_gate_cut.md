# The tail gate of a marker with no density minimum, where it can be placed

Above the dominant mode, openCyto's default side, never below zero
signal (a negative population ending below zero puts the cut at zero),
and only where it leaves under 60% of the parent above it (more means
the dominant mode was not a negative one). NOTHING above the cut is an
answer rather than a failure: a viability dye in a sample with no dead
cells, a marker no cell in this sample carries. Shared by the default
path of
[`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
and the last resort of
[`best_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/best_threshold.md),
so the two cannot disagree on where a tail goes.

## Usage

``` r
tail_gate_cut(x, floor = -Inf, side = c("upper", "lower"))
```

## Arguments

- x:

  transformed values within the parent gate.

- floor:

  the transformed value of zero signal, `-Inf` for scatter.

- side:

  `"upper"` (default), or `"lower"` where the analyst has declared the
  dominant mode positive; the cut is then the mode minus 3.5 SD and must
  leave over 40% of the parent above it.

## Value

the cut, or NA where no tail gate can be placed.

## Details

WHY THE SIDE IS NOT GUESSED. A single mode does not say whether it is
the marker's negative or positive population. Its distance from zero
signal does not either: granulocyte autofluorescence and the spillover
spread of bright cells put a NEGATIVE mode well above zero, and a
per-marker cofactor small enough to resolve the negatives stretches that
offset into many standard deviations. A version that cut below any mode
lying 3.5 of its SDs above zero took the lower side 17 times in two runs
on one 24-sample whole-blood cohort scored against manual gating, and
was wrong every time it could be checked: Vd1, Vd2 and CD11c, all rare,
called positive on 90-99% of CD45+ cells, emptying CD4 and CD8 T cells
(defined Vd1- Vd2-), and a CD16 cut below the bulk that took classical
monocytes to 0% and non-classical to 12% (manual 6.4% and 0.1%).
openCyto takes the side from the analyst for this reason. A marker
carried by most of its parent with no density minimum gets a cut above
its mode, flagged `needs_review`, unless the analyst declares the other
side, as openCyto's `side` argument does:
`thresholds: <marker>: {tail: lower}` in the configuration. A control, a
fixed cut or a `parent:` in which the marker separates are the
alternatives.
