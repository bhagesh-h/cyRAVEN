# The Cramer's V above which batch correction is refused

Cohen's convention puts a "large" association at w = 0.5, and Cramer's V
is w scaled by sqrt(min(levels) - 1), so the same large association is a
smaller V in a bigger table: 0.50 for two batches or two groups, 0.35
for three, 0.29 for four (Cohen 1988). A fixed V of 0.6 admitted a large
association in every table.

## Usage

``` r
cramers_v_limit(batch, group)
```

## Arguments

- batch, group:

  the two labelings.

## Value

the threshold.
