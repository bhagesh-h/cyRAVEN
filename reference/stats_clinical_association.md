# Clinical variables against population frequencies and marker intensities

MARKER INTENSITY IS TESTED WITHIN A POPULATION. A marker's level is a
property of a defined cell population, which is how a differential-state
test asks it (diffcyt-DS, Weber et al. 2019;
[`stats_marker_state()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_marker_state.md)
here), so each population and marker pair is one test, keyed
`"<population>, <marker>"`. This used to take, per sample and marker,
the median of the declared populations' medians: populations nest and
the list is the specification's, so that value followed whichever
population type the specification listed most (for HLA-DR, the T-cell
subsets) and described no set of cells. The pipeline passes
functional_markers.csv, whose pairs are the ones the gating document
scopes as interpretable, and population_marker_mfi.csv only where no
functional block resolved.

## Usage

``` r
stats_clinical_association(
  freq,
  mfi = NULL,
  clin = list(),
  patient_of = NULL,
  value_col = NULL,
  min_n = 4L
)
```

## Arguments

- freq:

  population_frequencies table.

- mfi:

  per-sample marker medians by population, with `sample_id`,
  `population`, `marker` and `median_asinh` (functional_markers.csv or
  population_marker_mfi.csv), optional.

- clin:

  named list of clinical variables, sample_id -\> value.

- patient_of:

  named vector sample_id -\> patient_id, used to report an association
  at the donor level where a variable is a property of the donor rather
  than of the acquisition. NULL treats every row as its own subject.

- value_col:

  the frequency column to use.

- min_n:

  fewest samples for a test.

## Value

list(populations, markers)
