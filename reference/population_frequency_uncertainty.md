# Propagate threshold uncertainty into one sample's population frequencies

For every marker a population's definition reads, the frequency is
re-scored with that marker's threshold displaced by plus and minus its
standard uncertainty, and half the resulting spread is taken as that
marker's contribution. Contributions are summed in quadrature across
markers, which is the GUM treatment of independent terms and is the same
arithmetic the operator studies apply to manual gates.

## Usage

``` r
population_frequency_uncertainty(
  tmat,
  thr,
  parent,
  spec,
  u,
  cd45_x = NULL,
  cd45_live = NULL,
  cd45_threshold = NA_real_,
  cd45_u = NA_real_,
  lod_events = 20L,
  loq_events = 50L
)
```

## Arguments

- tmat:

  transformed marker matrix for one sample

- thr:

  named threshold vector for that sample

- parent:

  logical parent-gate mask

- spec:

  population specification

- u:

  named vector of standard uncertainties, one per marker

- cd45_x:

  CD45 values for the sample, or NULL when the gate was skipped

- cd45_live:

  logical mask of live cells, the parent of the CD45 gate

- cd45_threshold:

  the CD45 cut

- cd45_u:

  standard uncertainty of the CD45 cut

- lod_events:

  events below which a population is not called detected

- loq_events:

  events below which a population is detected but not quantified

## Value

list(per_population = data.frame, budget = data.frame)

## Details

THE PARENT TERM IS INCLUDED AND MATTERS MOST. Displacing the CD45 cut
moves cells into and out of the denominator every population is
expressed against, so it perturbs every population at once. The operator
work finds the first gate dominates the budget; reporting the per-marker
terms without it would leave out the largest one.

Markers whose uncertainty is NA contribute nothing and are counted in
`n_terms_missing`, so a small total that is small only because most
terms could not be computed is distinguishable from a genuinely tight
one.

THE COUNTING TERM IS REPORTED BESIDE THE GATE TERMS, NOT MIXED INTO
THEM. `u_pct_points` keeps its meaning: gate placement only.
`u_counting_pct_points` is what the frequency carries from the number of
events behind it, and `u_total_pct_points` is their quadrature sum.
Keeping the first column fixed means every number this table published
before is still the same number.

The two are not strictly independent, since displacing a cut also
changes the count. Quadrature treats them as though they were, which is
the same approximation the GUM makes for the marker terms and is stated
here rather than left for the reader to find.
