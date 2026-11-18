# Pick the best available abundance measure, and say what it means

WHY THIS IS A FUNCTION AND NOT A HARDCODED COLUMN: three quantities in
the frequency table look interchangeable and are not.

## Usage

``` r
abundance_measure(freq, prefer = c("auto", "absolute", "frequency"))
```

## Arguments

- freq:

  Population frequency table, one row per sample x population.

- prefer:

  Preferred choice; auto picks the best available. Default
  `c("auto", "absolute", "frequency")`.

## Details

`count` EVENT count. Depends on how long the operator ran the tube and
on –max-events-per-file. NEVER comparable between samples; never a valid
y-axis for a group comparison. `pct_of_cd45_pos` Comparable, but
COMPOSITIONAL: populations are constrained to sum to 100%, so one
lineage expanding mathematically forces the others down. A significant
fall in "% NK" is equally consistent with NK loss and with granulocyte
gain. `cells_per_ul` Absolute. Independent per population, so a change
in one says nothing about the others. Requires wbc_per_ul.

Prefer absolute when available and fall back to frequency otherwise,
carrying the axis label and the compositional caveat with the choice so
the figure cannot silently misrepresent which of the two it is showing.
