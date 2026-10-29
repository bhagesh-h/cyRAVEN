# Find the events sitting at the top of a channel's range

WHAT: an event is at the margin when, in any of the channels `used`, it
is at or above that channel's highest value (`margin_at` from
[`linearise_channels()`](https://bhagesh-h.github.io/cyRAVEN/reference/linearise_channels.md)).

## Usage

``` r
margin_events(ex, scale, used, warn_frac = 0.01, label = "")
```

## Arguments

- ex:

  Expression matrix, linear units.

- scale:

  The `scale` table from
  [`linearise_channels()`](https://bhagesh-h.github.io/cyRAVEN/reference/linearise_channels.md).

- used:

  Column indices judged: the marker and scatter channels.

- warn_frac:

  Warn when more than this fraction of events is at a margin. Default
  `0.01`.

- label:

  File name for the log line.

## Value

list(event = logical per event, n = count, by_channel = named integer
count per judged channel with any)

## Details

WHY THEY ARE REMOVED RATHER THAN KEPT. A saturated event has no true
value in the saturated channel: it is somewhere above the top, by an
unknown amount. flowAI removes "margin events in the upper limit of the
dynamic range" (Monaco et al. 2016) and PeacoQC's `RemoveMargins()` does
the same (Emmaneel et al. 2022). Capping them instead, which flowCore's
`read.FCS()` and earlier versions of this package did, keeps the event
and invents a value for it: a pile of identical values at the top of the
scale. Floating-point files can also carry values far above the declared
range, up to 6.8e8 under a declared 262,144 on one cohort, which are not
measurements at all.

On the 30 spectral files of one sepsis cohort the margin events were
0.01% to 0.7% of each file. A file losing more than 1% is named in a
warning, because at that rate the voltage was set too high or the
declared range is wrong.
