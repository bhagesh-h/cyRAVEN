# Analysis-scale values back to instrument units

Inverts the run's transform marker by marker. A scatter channel that
explore mode put on a log10 scale is inverted as such; a marker the
transform has no parameters for comes back NA rather than through a
wrong inverse.

## Usage

``` r
mfi_back_transform(y, marker, tr, scatter = character(0))
```

## Arguments

- y:

  Numeric vector on the analysis scale.

- marker:

  Character vector, same length, naming each value's marker.

- tr:

  A transform from
  [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md).

- scatter:

  Markers that are log10 scatter rather than fluorescence.

## Value

Numeric vector in instrument units.
