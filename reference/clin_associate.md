# Associate per-sample values with clinical variables

Associate per-sample values with clinical variables

## Usage

``` r
clin_associate(d, key_col, value_col, clin, patient_of = NULL, min_n = 4L)
```

## Arguments

- d:

  long data frame with `sample_id`, a key column and a value column.

- key_col:

  column naming the thing measured: "population" or "marker".

- value_col:

  the per-sample quantity.

- clin:

  named list: variable -\> named vector of sample_id -\> value.

- patient_of:

  named vector sample_id -\> patient. When given, a variable that is
  CONSTANT within a patient is tested on one value per patient rather
  than one per sample. See
  [`clin_variable_unit()`](https://bhagesh-h.github.io/cyRAVEN/reference/clin_variable_unit.md).

- min_n:

  fewest observations for a test to run at all.

## Value

data frame, one row per key x variable, BH-adjusted within variable.
