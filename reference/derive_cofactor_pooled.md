# Derive the panel cofactor from SEVERAL samples, not from whichever file was read first

WHY THIS EXISTS. The original call site derived the cofactor from
`reads[[sids[1]]]` alone – one sample sets the arcsinh transform for the
whole panel, and therefore the scale on which every threshold, every
marker median and every UMAP distance in the run is computed. If that
first file happens to be weakly stained, an unstained control, or simply
an outlier, the transform is wrong for everyone else and nothing
downstream can detect it. The choice of "first" is alphabetical, so it
is not even a considered sample – it is a filename.

## Usage

``` r
derive_cofactor_pooled(reads, sids, max_samples = 8L, warn_ratio = 2, ...)
```

## Arguments

- reads:

  named list of read objects (each carrying exprs + marker_cols)

- sids:

  sample ids belonging to this panel

- max_samples:

  Ceiling on the number of samples used. Default `8L`.

- warn_ratio:

  Ratio of per-sample values above which a warning is emitted. Default
  `2`.

- ...:

  Passed to derive_cofactor().

## Value

single numeric cofactor, with attributes recording what it came from

## Details

WHAT IT DOES INSTEAD: derives the cofactor independently on up to
`max_samples` samples and takes the median of those. The median, not the
mean, because the failure this guards against is precisely one aberrant
sample, and a mean would let it back in.

WHY IT ALSO REPORTS THE SPREAD: if per-sample cofactors vary by more
than `warn_ratio`, the panel does not have one shared background and a
single cofactor is a compromise rather than a description. That is worth
saying out loud – it usually means a gain change or an instrument
setting drifted mid-batch, which is the same thing the batch diagnostic
looks for from the other end.

BEHAVIOUR CHANGE, stated plainly: this alters derived numbers relative
to the baseline whenever samples disagree. Pass
–cofactor-from-first-sample to restore the previous single-sample
behaviour exactly.
