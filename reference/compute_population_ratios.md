# Derived abundance ratios (e.g. CD4:CD8), from the config's `ratios:` block

WHY CONFIG-DRIVEN, NOT HARDCODED: population names live in
`populations:` (section 5.4) and are study-specific – a ratio hardcoded
to "CD4 T cells"/"CD8 T cells" would silently mean nothing (or error)
under a config using different labels. Each ratio supplies its own
numerator/denominator population name, the same way each
functional_blocks entry supplies its own marker list, so this stays
generic across cohorts and panels.

## Usage

``` r
compute_population_ratios(freq, ratios)
```

## Arguments

- freq:

  population_frequencies table (sample_id, population, pct_of_cd45_pos,
  is_control)

- ratios:

  named list from cfg\$ratios; each entry has `numerator` and
  `denominator` (population names from `populations:`) and an optional
  `label`

## Value

data.frame(sample_id, ratio, population (=label), numerator,
denominator, numerator_pct, denominator_pct, value, is_control), or NULL
