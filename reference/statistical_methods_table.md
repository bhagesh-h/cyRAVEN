# The methods catalogue

Every commonly reported test in the immunophenotyping and cytometry
literature, what this run did about it, and why. Written on every run so
the reasoning travels with the numbers.

## Usage

``` r
statistical_methods_table(
  n_groups = NA_integer_,
  paired = FALSE,
  n_tests = NA_integer_
)
```

## Arguments

- n_groups:

  Number of group levels compared, or NA.

- paired:

  Whether a paired design was detected.

- n_tests:

  Number of tests corrected together, or NA.

## Value

data.frame.
