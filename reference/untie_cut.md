# Move a cut off a block of identical values

Populations are scored with strict inequalities, `>` for above and `<`
for below, so an event whose value EQUALS the cut is in neither. On
continuous data that never happens; on a pile of identical values it
decides the result. A quantile of data whose top 15% sit on one value IS
that value, and `> cut` then selects nothing: on the demo cohort that
emptied every population of 11 samples. Integer channel data, and the
capped values of `--keep-margin-events`, are piles of exactly this kind.

## Usage

``` r
untie_cut(x, cut, target_above = NA_real_)
```

## Arguments

- x:

  the values the cut is applied to.

- cut:

  the proposed cut.

- target_above:

  intended fraction of `x` above the cut, or NA.

## Value

the cut, unchanged unless some value of `x` equals it.

## Details

The cut is moved to the midpoint between the tied value and its nearest
neighbour on one side, so the block falls wholly on the other. With a
`target_above` (a quantile's intended positive fraction) the side giving
the fraction nearer the target is taken; without one the block is kept
below the cut, which leaves a density-derived cut where it was to within
one step.
