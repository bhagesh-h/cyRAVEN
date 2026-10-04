# Clinical variables against population frequencies and marker intensities

Clinical variables against population frequencies and marker intensities

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

  population_marker_mfi table, optional.

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
