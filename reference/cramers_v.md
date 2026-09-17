# Cramer's V between two categorical variables

Used here to quantify how far batch and biological group overlap. 0 is
independent, 1 is one perfectly determined by the other.

## Usage

``` r
cramers_v(a, b)
```

## Arguments

- a, b:

  character or factor vectors of equal length

## Value

numeric in 0..1, or NA when either has a single level
