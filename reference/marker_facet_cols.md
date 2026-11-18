# Categorical columns worth splitting a per-marker UMAP by

WHY IT AUTO-DETECTS RATHER THAN TAKING ONLY `--group-column`. A cohort
usually carries several categories that all deserve looking at –
timepoint, infection focus, phenotype, sex – and naming one of them as
the statistical group column says nothing about which are worth SEEING.
A run grouped by timepoint would otherwise draw a single pooled panel
for infection focus, and the category would be invisible in the figures
even though it is in the sheet and on the UMAP overview.

## Usage

``` r
marker_facet_cols(
  cells,
  feature_cols = character(0),
  max_levels = 8L,
  max_cols = 4L
)
```

## Arguments

- cells:

  Embedded cell table.

- feature_cols:

  Marker columns, excluded from consideration.

- max_levels:

  Most distinct values a column may have and still be split by.

- max_cols:

  Most categories to draw, so the folder cannot explode. What is dropped
  is returned in the `skipped` attribute rather than silently lost.

## Value

Character vector of column names, with a `skipped` attribute.

## Details

WHAT IS EXCLUDED, AND WHY. Numeric columns: faceting an age in years
gives one facet per distinct age. Identifiers (`sample_id`,
`patient_id`, `file`): they are not categories, and `patient_id` alone
would multiply the folder by the donor count. Anything above
`max_levels`: past that the facets are too narrow to read, and the
figure stops being a comparison.
