# Where the live cells sit on a viability dye

The tallest mode of the dye's density, which is the live cells in any
sample with more live cells than dead ones. The lowest substantial mode
was tried instead, to cover a mostly dead sample, and failed on real
whole blood: there the dye often shows a second live mode near zero, as
tall as a third of the main one, and the cut fell between the two live
modes. A sample with more dead than live cells is the case this misses:
its live/dead valley lies below the dominant mode and is not taken, and
the tail gate placed above that mode instead is flagged for review.

## Usage

``` r
live_mode(x)
```

## Arguments

- x:

  numeric vector.

## Value

the mode, or NA with fewer than 200 finite values.
