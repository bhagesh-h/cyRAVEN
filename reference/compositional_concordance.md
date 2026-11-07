# Concordance between the raw-percentage test and the CLR test

WHAT THE READER GETS: for every population x comparison, whether the
conclusion survives removing the compositional constraint, with a
one-word verdict naming which of the four cases it falls into.

## Usage

``` r
compositional_concordance(
  raw_stats,
  clr_stats,
  alpha = 0.05,
  p_col = c("p_value", "p_adj_BH")
)
```

## Arguments

- raw_stats:

  The raw stats.

- clr_stats:

  The clr stats.

- alpha:

  Significance threshold. Default `0.05`.

- p_col:

  Which p-value column to compare on. Default
  `c("p_value", "p_adj_BH")`.

## Details

WHY A TABLE AND NOT A LOG LINE: "significant on raw but not on CLR" is a
per-population fact and populations behave differently. A single summary
sentence would let a reader carry a global impression onto the one
population where it does not hold.
