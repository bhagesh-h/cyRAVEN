# Estimate logicle parameters for one marker

WHAT: the standard automatic-logicle rule (Parks, Roederer & Moore 2006;
Moore & Parks 2012). `t` is the top of scale, `m` the number of decades,
and `w` the width of the linear region, chosen so that the linear region
just covers the spread of the NEGATIVE population:

## Usage

``` r
logicle_params_one(x, m = 4.5, a = 0, neg_q = 0.05)
```

## Arguments

- x:

  numeric vector of raw intensities for one marker

- m:

  decades on the display scale

- a:

  additional negative decades to display

- neg_q:

  quantile of the negative values used as `r`

## Value

list(w, t, m, a)

## Details

    w = (m - log10(t / |r|)) / 2

where `r` is a low quantile of the values below zero.

WHY A QUANTILE OF THE NEGATIVES AND NOT THE MINIMUM: the minimum is one
event. On a few hundred thousand events it is whatever the noisiest cell
in the tube did, and `w` derived from it stretches the linear region
until the positive decades are squashed. The 5th percentile of the
negative values describes the bulk of that population instead of its
worst member.

Returns `w = 0` when a channel has no negative values at all – with
nothing below zero there is nothing for a linear region to rescue, and
logicle degenerates to a log transform, which is the right answer rather
than a failure.
