# Pick black or white text per tile, whichever the fill can actually be read against

WHY: a single fixed text colour cannot work on a sequential fill. The
composition heatmap used one grey for every label, so numbers sitting on
the dark end of viridis were effectively invisible – grey on near-black
– while the same grey was fine on the bright end. Which cells become
unreadable depends on the data, so it is not something a fixed colour
can be chosen around.

## Usage

``` r
contrast_text_colour(values, limits, option = "D", dark = "grey10")
```

## Arguments

- values:

  numeric vector

- limits:

  the fill scale's limits; MUST match the scale, so the colour computed
  here is the colour ggplot draws

- option:

  viridis option name

- dark:

  the dark ink to use when the fill is light

## Details

HOW: compute the fill each value will actually receive, take its WCAG
relative luminance, and return whichever of white / near-black has the
higher contrast ratio against it. Comparing the two ratios rather than
thresholding on luminance means the answer stays correct if the palette
is re-themed – the viridis option is configurable, and "cividis" or
"magma" put their dark end in a different place than "D" does.
