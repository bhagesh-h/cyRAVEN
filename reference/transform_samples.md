# The samples a panel's transform is estimated from

The biological samples, never the declared controls. A control is in the
panel to be compared against, not to set the scale everything is
measured on: an unstained tube drags the cofactor down, and a stained
FMO or full-stain tube is one donor's distribution counted again.
Including them also changes WHICH samples the pooled estimate draws,
since it takes an evenly spaced subset of the panel. On the 30-sample
sepsis cohort adding two control tubes moved the cofactor from 261 to
326 and every cut in every sample with it, in markers the controls had
nothing to do with. When a panel holds nothing but controls they are all
it has, and are used.

## Usage

``` r
transform_samples(sids, smap = NULL)
```

## Arguments

- sids:

  sample ids in the panel.

- smap:

  the sample map, or NULL.

## Value

the sample ids to estimate from.
