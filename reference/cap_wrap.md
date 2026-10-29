# Wrap a caption to the canvas it will be drawn on

WHY IT TAKES THE WIDTH. A caption is drawn at one text size with no
wrapping of its own, so a line longer than the canvas is silently
CLIPPED – the figure still writes and the end of the sentence is gone. A
fixed wrap width therefore only works for a fixed canvas, and several of
these figures size themselves from the data. Roughly seventeen
characters fit per inch at the caption's size, with the constant kept
slightly conservative because the exact figure depends on the device's
font metrics.

## Usage

``` r
cap_wrap(txt, width_in)
```

## Arguments

- txt:

  one or more strings, pasted together.

- width_in:

  the canvas width in inches, the same value passed to ggsave.
