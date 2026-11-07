# Cap values at the top of each channel's range

The behaviour of earlier versions, kept for `--keep-margin-events`:
every value above a channel's range is set to that range, the clamping
flowCore's `read.FCS()` does by default. Only channels whose range is
known are capped, so a floating-point channel whose `$PnR` describes
some other scale (see
[`linearise_channels()`](https://bhagesh-h.github.io/cyRAVEN/reference/linearise_channels.md))
is left alone rather than cut in the middle of its data.

## Usage

``` r
clamp_to_declared_range(
  ex,
  scale,
  used = seq_len(ncol(ex)),
  warn_frac = 0.01,
  label = ""
)
```

## Arguments

- ex:

  Expression matrix, linear units.

- scale:

  The `scale` table from
  [`linearise_channels()`](https://bhagesh-h.github.io/cyRAVEN/reference/linearise_channels.md).

- used:

  Column indices to cap. Default every column.

- warn_frac:

  Name a channel in a warning when more than this fraction of its events
  was capped. Default `0.01`.

- label:

  File name for the log line.

## Value

`ex` with the out-of-range values capped.
