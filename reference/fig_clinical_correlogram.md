# Spearman correlation between the clinical variables themselves

[`stats_clinical_association`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_clinical_association.md)
adjusts p-values within each variable, on the grounds that each variable
is its own question. That holds only when the variables carry different
information. A cohort where the sickest patients are also the ones who
died has one gradient and two columns describing it, and an association
found against both is one finding reported twice. This figure is where
that is visible.

## Usage

``` r
fig_clinical_correlogram(clin, outfile, dpi = 200, colors = fcs_colors())
```

## Arguments

- clin:

  named list of clinical variables, sample_id -\> value.

- outfile:

  path.

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.

## Details

Numeric variables are included, and two-level variables coded 0/1 – for
which Spearman is the rank-biserial correlation, a legitimate quantity.
Variables with three or more unordered levels are excluded and named in
the caption: there is no ordering to correlate, and coding them 1/2/3
would invent one.

The circle area is the absolute coefficient and the fill its sign. The
number is printed only where the unadjusted p-value is below 0.05; the
circle is drawn either way, because an unlabelled small circle says
"measured, not distinguishable from zero", which is a different
statement from a blank cell.
