# Phenotype string for each cluster

`CD19+ HLA-DR+ CD3- CD14-`, from the positivity fractions. A marker is
called `+` above `hi`, `-` below `lo`, and omitted in between – an
omitted marker is an honest "this cluster is mixed for it", not a silent
negative.

## Usage

``` r
explore_phenotype(pos, hi = 0.65, lo = 0.2, max_markers = 6L)
```

## Arguments

- pos:

  Matrix from
  [`explore_positivity()`](https://bhagesh-h.github.io/cyRAVEN/reference/explore_positivity.md).

- hi, lo:

  Fractions bounding a positive and a negative call.

- max_markers:

  Most markers to name, highest \|evidence\| first.

## Value

Character vector, one phenotype per cluster.
