# Canonical text form of a population specification

Used to detect that the specification itself changed between runs, which
makes a threshold comparison meaningless: a population defined
differently is a different population, however similar its frequency.

## Usage

``` r
spec_fingerprint(spec)
```

## Arguments

- spec:

  population specification

## Value

named character vector, one canonical string per population

## Details

The full text is stored rather than a hash, so a mismatch can name the
population that changed instead of only asserting that something did.
