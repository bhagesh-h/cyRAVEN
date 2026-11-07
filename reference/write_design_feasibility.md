# Report which group comparisons can be made, and which cannot

Writes `design_feasibility.csv`: one row per group, carrying the sample
and donor counts behind it, whether the pipeline will test it, whether
an unpaired test of it is defensible, and the reason when those two
disagree.

## Usage

``` r
write_design_feasibility(
  group_of,
  gcol,
  min_n,
  reference,
  outdir,
  donor_of = NULL,
  tested = TRUE
)
```

## Arguments

- group_of:

  Named vector, sample_id to group.

- gcol:

  Name of the column the grouping came from.

- min_n:

  Minimum samples per group required to test.

- reference:

  Reference group, or NULL.

- outdir:

  Output directory.

- donor_of:

  Named vector, sample_id to donor. NULL when unavailable, in which case
  the repeated-measures check is not attempted and says so.

- tested:

  FALSE under `--no-group-tests`.

## Value

invisible data.frame, or NULL when there is no grouping to report.
