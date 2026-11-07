# One population's two subsets against each other, per timepoint

WHY A DEDICATED FIGURE. A pair of subsets that divide one compartment –
Vd1 and Vd2 inside the gamma-delta compartment, CD4 and CD8 inside T
cells – is read as a balance, not as two independent numbers, and the
published comparisons report it that way: the subsets shift against each
other while the compartment stays flat. Two separate panels make that
shift hard to see, because the reader has to hold one panel in mind
while looking at the other.

## Usage

``` r
fig_subset_balance(
  freq,
  smap,
  outfile,
  pairs = list(c("Vd1 T cells", "Vd2 T cells"), c("CD4 T cells", "CD8 T cells")),
  value_col = "pct_of_cd45_pos"
)
```

## Arguments

- freq:

  population_frequencies table.

- smap:

  Sample map carrying `timepoint` and `patient_id`.

- outfile:

  Destination PNG.

- pairs:

  List of two-element character vectors naming the subsets.

- value_col:

  Column to plot.
