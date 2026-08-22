# Paired / repeated-measures abundance test

WHAT: when samples pair up (pre/post, longitudinal, matched donors),
compares WITHIN pair instead of between groups.

## Usage

``` r
stats_paired_comparison(
  freq,
  pair_of,
  cond_of,
  reference = NULL,
  value_col = NULL,
  min_pairs = 3L
)
```

## Arguments

- freq:

  Population frequency table, one row per sample x population.

- pair_of:

  named vector sample_id -\> pairing unit (donor)

- cond_of:

  named vector sample_id -\> condition (timepoint / arm)

- reference:

  The group every other group is compared against.

- value_col:

  Column to test or plot instead of the default measure.

- min_pairs:

  Minimum complete pairs before a paired test is attempted. Default
  `3L`.

## Details

WHY IT CANNOT BE INFERRED AND MUST BE DECLARED: nothing in an FCS file
or a filename says two tubes came from the same donor at two timepoints.
Treating paired samples as independent – which is what the baseline
necessarily does, having no notion of pairing – discards the pairing's
variance reduction AND violates the independence assumption the rank
test rests on. So this runs only when –paired-column names a column, and
it counts incomplete pairs rather than silently dropping them, because a
pairing that quietly halves n is exactly the failure this function
exists to prevent.

WHY WILCOXON SIGNED-RANK / FRIEDMAN: the paired counterparts of the rank
tests used everywhere else here. Friedman generalises to 3+ conditions
and is the conventional choice for repeated-measures frequency data at
this scale.
