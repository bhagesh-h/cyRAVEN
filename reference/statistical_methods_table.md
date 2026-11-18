# The methods catalogue

Every commonly reported test in the immunophenotyping and cytometry
literature, what this run did about it, and why. Written on every run so
the reasoning travels with the numbers.

## Usage

``` r
statistical_methods_table(
  n_groups = NA_integer_,
  paired = FALSE,
  n_tests = NA_integer_,
  ran = list()
)
```

## Arguments

- n_groups:

  Number of group levels compared, or NA.

- paired:

  Whether a paired comparison was written.

- n_tests:

  Number of tests corrected together, or NA.

- ran:

  Named list of logicals saying which analyses wrote output:
  `group_tests`, `parametric`, `compositional`, `uncertainty`,
  `normality`, `batch`, `confounding`, `rank_ancova`, `edger`. An absent
  entry takes the value a grouped run with the default options has.

## Value

data.frame.

## Details

WHAT "USED" MEANS. A method is listed as used only where this run wrote
its result, which `ran` records from the files actually written. The
parametric tests are written by default beside the rank tests
(parametric_tests.csv, posthoc_tests.csv), so they are used, as
secondary results; listing them as not used, as this table once did,
contradicted the run's own outputs.
