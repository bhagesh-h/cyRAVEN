# Epsilon-squared, the rank effect size for Kruskal-Wallis

WHY IT WAS ADDED. A three-level variable such as infection focus
produced a p-value and nothing else, so it appeared in the heatmap as a
grey tile: the reader could see that a test had been run and not whether
it had found anything. Epsilon-squared is `H / (n - 1)`, the standard
rank effect size for this test (Tomczak & Tomczak 2014), and is the
proportion of rank variance the grouping accounts for. It is bounded 0
to 1 and UNSIGNED, because a variable with three levels has no single
direction – which is why it is carried in its own column and never
coloured on the signed scale.

## Usage

``` r
clin_epsilon2(h, n)
```

## Arguments

- h:

  Kruskal-Wallis statistic.

- n:

  total number of observations.
