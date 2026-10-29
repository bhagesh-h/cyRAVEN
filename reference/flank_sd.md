# Standard deviation of a density mode from its half-maximum width on one side

The distance from the mode to where the density falls to half its
height, interpolated between grid points, is the half-width at half
maximum; for a normal population it is sqrt(2 ln 2) standard deviations.
A kernel density estimate of a normal population with SD `s` is normal
with variance `s^2 + bw^2`, so the kernel's variance is subtracted.

## Usage

``` r
flank_sd(d, i, side = c("upper", "lower"))
```

## Arguments

- d:

  a [`stats::density()`](https://rdrr.io/r/stats/density.html) result.

- i:

  index of the mode in `d$y`.

- side:

  `"upper"` measures the flank above the mode, `"lower"` below.

## Value

the SD, or NA where the density never falls to half on that side.
