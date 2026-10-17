# How many samples each population is actually measurable in

One bar per population, split by whether that population's event count
in each sample was enough to quantify it, enough to detect it, or
neither.

## Usage

``` r
fig_detection_limits(ufreq, outfile, dpi = 200, colors = fcs_colors())
```

## Arguments

- ufreq:

  the `frequencies` element of
  [`run_gate_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_gate_uncertainty.md)

- outfile:

  Path to write the figure to.

- dpi:

  Resolution in dots per inch. Default `200`.

- colors:

  Named list of colours; defaults to the package palette. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHY A COUNT RATHER THAN A SCATTER. The limit of detection depends only
on how many parent-gate events a sample contributed, so it is a property
of the sample and identical across the populations within it. Plotting
frequencies against it would put every population on the same reference
line and say nothing extra. What varies, and what decides whether a
population is worth testing, is how many samples clear the limit at all.

HOW TO READ IT. A population quantified in every sample is measurable in
this cohort. One that is mostly below the limit of quantification is
not, and no change to the gating strategy will fix it: the numerator is
small because few cells were acquired, so the answer is a longer
acquisition or a higher `--max-events-per-file`, not a different
threshold. A population split between the two is the dangerous case,
because its group difference can be driven entirely by which samples
happened to clear the limit.
