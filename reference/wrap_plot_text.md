# Hard-wrap a title/subtitle/caption to the figure's actual width

WHY THIS IS NEEDED AT ALL: ggplot2 does not wrap plot titles, subtitles
or captions. A long one is drawn as a single line and CLIPPED at the
device edge, silently – the text is simply gone, with no warning and no
visual cue that anything is missing. On these heatmaps that lost the
last third of both subtitles, including the sentence explaining what an
even split is.

## Usage

``` r
wrap_plot_text(txt, width_in, pt = 9.5, margin_in = 0.35)
```

## Arguments

- txt:

  the string; any newlines already in it are preserved and each
  resulting line is wrapped independently (captions rely on this)

- width_in:

  the figure width being passed to ggsave

- pt:

  point size the text will be rendered at

- margin_in:

  Margin to reserve, in inches. Default `0.35`.

## Details

It is worse than an ordinary layout bug because the truncation lands
mid-sentence and reads as if the sentence ended there. A reader has no
way to tell a clipped subtitle from a badly written one.

HOW THE WIDTH IS ESTIMATED: at a given point size the mean advance width
of the default sans face is close to 0.52 \* size, so
characters-per-inch is 72 / (0.52 \* pt). `margin_in` covers the plot
margins and any inset. This is an estimate, not a measurement – strwrap
is given a deliberately conservative character count so a slightly
wider-than-average string still fits rather than spilling. Callers pair
it with plot.title.position = "plot", which starts the text at the
figure edge instead of the panel edge and so makes the full width
available.
