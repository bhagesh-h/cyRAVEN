# Group files into panels by their exact marker set

WHY: markers define the feature space. Two files with different marker
sets have no common space, so binding them and embedding gives
meaningless coordinates – the original template's blind row-bind
produced NA columns for every non-shared channel and then aborted on its
own finiteness check.

## Usage

``` r
fingerprint_panels(reads, labels = NULL, optional = character(0))
```

## Arguments

- reads:

  Named list of objects returned by
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md).

- labels:

  Vector of labels, one per row of coords.

- optional:

  Marker names that do NOT take part in the fingerprint, so a reagent
  present in only some files does not split the cohort into two panels.
  Scored where present, UNAVAILABLE where absent. Default none.

## Value

list(assignment = named character (sample_id -\> panel), panels = list)
