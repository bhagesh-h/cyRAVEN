# The detector a channel belongs to, with its pulse-feature suffix removed

`FL1.A` and `FL1.H` are two measurements of one detector, and so is the
pair `X..RB613...A` / `X..RB613...H` that `alter.names` makes of a
spectral channel. Stripping the suffix gives the key the two share.

## Usage

``` r
detector_of(nm)
```

## Arguments

- nm:

  `$PnN` values, as read.

## Value

character vector.
