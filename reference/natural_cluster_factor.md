# Order cluster labels k1, k2, ... k30 rather than k1, k10, k11, ... k2

Cluster ids are "k" plus an integer, and sorting them as text puts k10
before k2. In a 30-cluster legend that is not a cosmetic problem: the
reader looks for k7 between k6 and k8 and finds it near the bottom, and
colour-to-cluster matching stops being possible at a glance.

## Usage

``` r
natural_cluster_factor(x)
```

## Arguments

- x:

  Character vector of cluster ids.

## Value

`x` as a factor in natural order.
