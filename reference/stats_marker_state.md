# Sample-level differential-state test: population x marker x group

WHAT: for every (population, marker) pair, compares each non-reference
cohort against the reference using the per-sample values already in
population_marker_mfi.csv. One row per pair per comparison.

## Usage

``` r
stats_marker_state(
  mfi,
  group_of,
  reference = NULL,
  measures = c("median_asinh", "pct_positive"),
  min_cells = 20L,
  min_n = 3L
)
```

## Arguments

- mfi:

  population_marker_mfi table (sample_id, population, marker, n_cells,
  median_asinh, pct_positive, and – after the deliverable change
  documented in the README – is_control and qc_status)

- group_of:

  named character vector sample_id -\> group label

- reference:

  group every other group is compared against

- measures:

  which per-sample quantities to test

- min_cells:

  per-sample floor on the population's cell count

- min_n:

  minimum samples per group for a test to be attempted

## Value

data.frame, one row per population x marker x comparison x measure, or
NULL when nothing is testable

## Details

WHY BOTH MEASURES ARE TESTED AND REPORTED SEPARATELY: `median_asinh` the
population's median intensity for that marker. Moves when the whole
population shifts up or down – the classic "activation marker is
brighter in patients" result. `pct_positive` the fraction of that
population above the sample's own gating threshold. Moves when a SUBSET
of the population turns positive while the rest does not. A bimodal
shift (30% of cells go bright, 70% unchanged) moves pct_positive sharply
and the median barely at all; a uniform shift does the opposite.
Reporting only one measure makes the other kind of change invisible, so
both are computed and the measure is a column, not an assumption.

WHY min_cells: the per-sample median of a marker over 12 cells is noise
shaped like a measurement. Rows below the floor are dropped BEFORE
testing rather than down-weighted, because a rank test has no weighting
mechanism and silently treats a 12-cell median as equal evidence to a
4,000-cell one.

WHY THE FDR FAMILY IS THE WHOLE TABLE: this tests every marker in every
population, which at a 12-marker panel and 12 populations is ~144 pairs
x the number of comparisons. That is a screen, and its multiplicity has
to be accounted for across the screen, not within each population.
