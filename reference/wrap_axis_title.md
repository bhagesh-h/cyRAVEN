# Wrap a panel's axis title so a long key does not overrun its neighbour

WHY THIS IS NEEDED. The panel key becomes the y-axis title, and a y-axis
title is rotated, so the room it has is the panel's HEIGHT rather than
its width. A population name fits comfortably. The key
[`fig_functional_markers()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_functional_markers.md)
supplies does not: it is "block: population, marker", which reaches 47
characters and more, and an unwrapped title of that length runs out of
the panel and lands on the axis numbers and tag letter of the panel
beside it.

## Usage

``` r
wrap_axis_title(x, width = 30L, max_lines = 3L)
```

## Arguments

- x:

  Character scalar, the title to wrap.

- width:

  Target characters per line before wrapping starts.

- max_lines:

  Ceiling on line count. Beyond this the lines get longer rather than
  more numerous, because a title taller than the panel is wide is worse
  than a slightly long one.

## Value

`x` with newlines inserted, unchanged when it already fits.

## Details

WHY IT REBALANCES RATHER THAN WRAPPING AT A FIXED WIDTH.
[`strwrap()`](https://rdrr.io/r/base/strwrap.html) at a fixed width
leaves the last line short, so a 47-character title at width 30 gives
one full line and one stub. Choosing the line count first and then
dividing gives lines of roughly equal length, which reads better rotated
and keeps the block from looking ragged.
