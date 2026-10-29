# Assign colours to population labels, reserving grey for the catch-all

WHY: "Other CD45+" is not a population – it is everything the gating
strategy failed to classify. Colouring it like a lineage invites reading
it as one, and when it is large (which is itself the signal to check the
thresholds) it dominates the figure in a saturated hue. Grey recedes;
that is the point. Fixed meaning-carrying colours, consulted before any
palette is issued

## Usage

``` r
semantic_colours(levels, reference = reference_group(), colors = fcs_colors())
```

## Arguments

- levels:

  Character vector of factor levels.

- reference:

  The group every other group is compared against. Default
  [`reference_group()`](https://bhagesh-h.github.io/cyRAVEN/reference/reference_group.md).

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHY: a colour that means one group in one figure and a different group
in the next is worse than no colour at all – the reader carries the
association across the page whether or not it holds. Sex and study are
the two variables drawn repeatedly across these outputs, so both get a
fixed assignment (male/healthy green, female/study-1 red, study-2
purple) that every discrete scale checks first. Everything else still
falls through to pop_palette().

Returns NULL when the levels are not a variable we pin, which is the
signal to use the ordinary palette. The study that anchors green is
package state, not an argument: it is set once by run_cyraven() from
–reference-group (see set_reference_group()) and read here as a DEFAULT
ARGUMENT, evaluated at call time. That is what lets every discrete scale
colour the cohorts consistently without threading `reference` through a
dozen wrapper functions – the same mechanism fcs_colors() uses.
