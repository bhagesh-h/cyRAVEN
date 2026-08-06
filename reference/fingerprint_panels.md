# Group files into panels by their exact marker set

WHY: markers define the feature space. Two files with different marker
sets have no common space, so binding them and embedding gives
meaningless coordinates – the original template's blind row-bind
produced NA columns for every non-shared channel and then aborted on its
own finiteness check.

## Usage

``` r
fingerprint_panels(reads, labels = NULL)
```

## Arguments

- reads:

  Named list of objects returned by
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md).

- labels:

  Vector of labels, one per row of coords.

## Value

list(assignment = named character (sample_id -\> panel), panels = list)
