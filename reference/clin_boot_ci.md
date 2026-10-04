# Percentile bootstrap interval for a rank effect size

WHY AN INTERVAL AND NOT ONLY A POINT ESTIMATE. On nine samples a rho of
0.61 and a rho of 0.05 can be the same underlying quantity, and the
point estimate alone cannot say so. The interval is the part of the
answer that carries how little the cohort constrains the effect, and
reporting the effect with its interval beside the p-value is the
standard recommendation for exactly this situation.

## Usage

``` r
clin_boot_ci(stat, n, B = 2000L, min_n = 6L, seed = 42L)
```

## Arguments

- stat:

  function of an index vector returning the effect.

- n:

  number of observations to resample.

- B:

  resamples.

- min_n:

  fewest observations for an interval to be reported at all.

- seed:

  fixed so the same table always yields the same interval.

## Details

WHAT IT IS NOT. A percentile bootstrap at this sample size is itself
approximate: below about nine observations the effective resample is
smaller than the nominal one and coverage falls short of 95%. The
interval is written anyway, because a visibly wide interval is a better
description of the evidence than a bare point estimate, and it is
suppressed below `min_n` rather than being quoted at a width nobody
should trust.

The RNG stream is saved and restored, so the interval is reproducible
for a given table no matter what sampled before it, and nothing
downstream of here sees a different stream than it would have.
