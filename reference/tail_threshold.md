# Cut a robust spread above (or below) the dominant mode

openCyto's `tailgate`. The spread is estimated by mirroring the half of
the distribution the positive tail cannot reach, so a long tail does not
inflate the width of the mode it is being measured against.

## Usage

``` r
tail_threshold(x, k = 3.5, side = c("upper", "lower"), adjust = 2)
```

## Arguments

- x:

  transformed values.

- k:

  how many robust standard deviations from the mode. Default `3.5`.

- side:

  `"upper"` for a minority-positive marker, `"lower"` for a marker whose
  dominant mode is itself the positive population.

- adjust:

  kernel bandwidth multiplier. Default `2`.

## Value

the cut.
