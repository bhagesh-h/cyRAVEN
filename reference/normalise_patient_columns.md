# Coerce, derive and translate the subject columns of a metadata table

The half of
[`load_patient_table()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_patient_table.md)
that operates on an already-selected data.frame: numeric coercion,
placeholder rejection, age derivation from date of birth, and dictionary
translation of free-text values.

## Usage

``` r
normalise_patient_columns(
  out,
  value_map = default_value_map(),
  reference_date = Sys.Date()
)
```

## Arguments

- out:

  data.frame whose columns already carry canonical English names.

- value_map:

  The value map. Default
  [`default_value_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_value_map.md).

- reference_date:

  The reference date. Default
  [`Sys.Date()`](https://rdrr.io/r/base/Sys.time.html).

## Value

`out` with the same columns, coerced and translated, carrying an
`untranslated` attribute.

## Details

WHY IT IS A SEPARATE FUNCTION. Subject attributes now arrive by two
routes: a standalone patient table, and the subject columns of a unified
sample sheet (see
[`read_samplesheet()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_samplesheet.md)).
Both must apply the same coercions, or the same study analysed through
the two routes would produce different ages and different sex labels.
Sharing the implementation is what makes the two routes equivalent
rather than merely similar.
