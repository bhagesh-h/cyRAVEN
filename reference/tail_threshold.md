# Cut a robust spread above (or below) the dominant mode

The robust form of the conventional positivity rule, the negative
population's centre plus k of its standard deviations. The centre is the
dominant mode. Its spread is estimated twice, each estimate immune to a
different contamination, and the smaller is used, since contamination
only ever widens a spread:

## Usage

``` r
tail_threshold(x, k = 3.5, side = c("upper", "lower"), adjust = 2)
```

## Arguments

- x:

  transformed values.

- k:

  how many robust standard deviations from the mode. Default `3.5`.

- side:

  `"upper"` for a minority-positive marker, `"lower"` for a marker whose
  dominant mode is itself the positive population.

- adjust:

  kernel bandwidth multiplier. Default `2`.

## Value

the cut.

## Details

- the MAD of the half the positive tail cannot reach, mirrored about the
  mode, which a positive tail or shoulder on the cut side cannot
  inflate;

- the half-width at half maximum of the mode's flank on the cut side,
  converted to a standard deviation (HWHM / sqrt(2 ln 2)) after removing
  the kernel's own variance, which a second negative population on the
  far side cannot inflate. Peak width at half maximum is cytometry's
  established measure of a population's spread (the half-peak CV).

The second exists because the first fails where the negative cells are
themselves two populations with different background, as monocytes and
lymphocytes are on PBMC: the mirrored half takes in the second hump and
the cut lands past a small positive shoulder, or past every cell. On one
sepsis cohort the MAD alone put CD19 beyond every cell in 3 of 30
samples. For a single normal negative the two estimates agree.

NOT openCyto's `gate_tail`, which places the cut where the density's
first derivative, past the reference peak's steepest descent, falls
below a tolerance (cytoUtils `.cytokine_cutpoint`). On a normal negative
mode the two land near each other. Under `auto_tol` the tolerance is 1%
of the steepest slope, which a normal reaches 3.57 SD out. The default
`tol = 0.01` is an absolute slope, so where it cuts depends on the
data's units and the population's share: 3.1 SD out of a unit-variance
normal holding every event, 2.5 SD at SD 2, 3.6 SD at SD 0.5. The slope
rule also depends on the kernel bandwidth; tried on the sepsis cohort it
called 55% of one sample's cells Vd2+ at a bandwidth where it found the
B cells, and was not adopted.
