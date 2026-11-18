# Put every channel on a linear scale and find the top of its range

WHY. An analog instrument with log amplifiers (a FACSCalibur, for one)
stores integer channel numbers 0 to `$PnR` - 1 on a log scale, and
`$PnE f1,f2` declares it: f1 decades, offset f2. Every transform this
package fits assumes linear intensity, so an arcsinh or logicle laid
over channel numbers puts a log on top of a log and moves every valley.
FCS 3.1 allows log storage only for integer data (Spidlen et al. 2010).

## Usage

``` r
linearise_channels(ex, kw, label = "")
```

## Arguments

- ex:

  Expression matrix, events x channels, in file order.

- kw:

  FCS keywords.

- label:

  File name for the log line.

## Value

list(exprs, scale), `scale` a data.frame with one row per channel: `f1`,
`f2`, `stored` ("linear", "log, linearised" or "log declared on floats,
used as stored"), `range` (top of scale in the units now held, or NA)
and `margin_at` (the value at and above which an event is at the margin,
or NA).

## Details

WHAT. Integer channels with f1 \> 0 are converted to f2 \* 10^(f1 \* x /
`$PnR`), with f2 = 0 read as 1, which is exactly what flowCore's
`read.FCS(transformation = "linearize")` computes. Their range becomes
f2 \* 10^f1.

FLOATING-POINT DATA WITH A NON-ZERO `$PnE` breaks the standard, and the
values are not converted, because flowCore does not convert them either
and nothing says what scale they are on. flowCore's own GvHD example,
written back out with `write.FCS()`, is like this: linear floats up to
10,000 under the FACSCalibur keywords `$PnE 4,0` and `$PnR 1024`. The
`$PnR` of such a channel describes the original channel numbers rather
than the stored values, so it is not used as the channel's range.
Capping those values at 1024 is what put a pile of identical values in
the middle of the demo cohort's CD45 data.
