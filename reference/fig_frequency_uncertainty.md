# Population frequencies with the uncertainty of the gate that produced them

One row per population, one point per sample, and a bar through each
point spanning plus and minus the standard uncertainty of that sample's
frequency propagated from where its thresholds were placed.

## Usage

``` r
fig_frequency_uncertainty(
  ufreq,
  outfile,
  group_of = NULL,
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- ufreq:

  the `frequencies` element of
  [`run_gate_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_gate_uncertainty.md)

- outfile:

  Path to write the figure to.

- group_of:

  Named character vector mapping sample_id to group label.

- dpi:

  Resolution in dots per inch. Default `200`.

- colors:

  Named list of colours; defaults to the package palette. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHAT TO LOOK FOR. The question this figure answers is whether the spread
BETWEEN samples is larger than the bar WITHIN each of them. Where it is,
the variation is something the gate is measuring. Where the bars are as
wide as the scatter, the variation is the gate moving, and a group
difference in that population is not interpretable however small its
p-value.

Controls and QC-failed samples are excluded here rather than labelled,
unlike the QC figures: their percentages are fractions of an arbitrary
parent slice, so plotting them next to real samples would put two
different quantities on one axis.
