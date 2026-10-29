# Report override entries that matched nothing

An override silently applying to no sample is the failure mode this
guards: the analyst believes a cut was corrected, the run says nothing,
and the uncorrected number is published.

## Usage

``` r
report_unused_overrides(overrides, applied)
```

## Arguments

- overrides:

  the `sample_overrides` list, or NULL

- applied:

  character vector of "sample\rmarker" keys that were used

## Value

invisibly, the character vector of unused keys
