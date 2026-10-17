# Save a ggplot at a DPI that cannot exceed the device's raster ceiling

A drop-in for
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
that lowers `dpi` – never raises it – when the requested canvas would
produce a raster larger than the graphics device can allocate (about
32,767 px per side).

## Usage

``` r
safe_ggsave(filename, plot, width, height, dpi = 300, ...)
```

## Arguments

- filename:

  Output path.

- plot:

  A ggplot or patchwork object.

- width, height:

  Canvas size in inches.

- dpi:

  Requested resolution; the effective value may be lower.

- ...:

  Passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Value

The path, invisibly.

## Details

Several figures scale their width or height with the number of samples,
populations or markers, so a large cohort or a wide panel can ask for a
canvas bigger than the device supports; the run would then die on a
figure rather than on the analysis. `ggsave()`'s own `limitsize` guards
the opposite mistake – an accidentally huge size in inches – and does
not help here.

Reducing DPI keeps the requested *physical* size, and therefore every
relative font, point and line size in the theme, exactly as specified.
At extreme scale the figure becomes less sharp; it does not become
smaller or more crowded, and the run finishes.

The ceiling is `getOption("cyRAVEN.max_raster_px")`, default 30000.
