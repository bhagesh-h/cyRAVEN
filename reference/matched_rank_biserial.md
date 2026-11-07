# Matched-pairs rank-biserial correlation

(R+ - R-) / (R+ + R-), with R+ and R- the sums of the ranks of \|d\|
over the positive and negative differences, zero differences dropped as
the signed-rank test drops them (Kerby 2014, Compr Psychol 3:11.IT.3.1).

## Usage

``` r
matched_rank_biserial(d)
```

## Arguments

- d:

  paired differences.

## Value

numeric in `[-1, 1]`, or NA when no difference is non-zero.
