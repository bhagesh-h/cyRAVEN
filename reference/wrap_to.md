# Wrap a subtitle to the width of the figure it sits on

A narrow heatmap of five markers is six inches wide, and a one-line
subtitle that explains the colour scale is longer than that; unwrapped,
ggplot2 runs it off the edge of the image.

## Usage

``` r
wrap_to(x, width_in)
```

## Arguments

- x:

  Subtitle text, or NULL.

- width_in:

  Figure width in inches.
