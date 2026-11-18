# Otsu's between-class variance threshold

Parameter-free and direction-neutral, and right where the two classes
are of comparable size. It splits a single population when the positive
fraction is small, which is why it is one candidate among several rather
than the rule.

## Usage

``` r
otsu_threshold(x, nbins = 512L)
```

## Arguments

- x:

  transformed values.

- nbins:

  histogram resolution. Default `512L`.

## Value

the cut.
