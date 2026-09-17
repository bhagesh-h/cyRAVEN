# How far each population moves if the flagged windows are excluded

The quantity that turns a QC flag into a decision. A file with a visible
anomaly that moves no population by more than its own gate uncertainty
does not need re-acquiring.

## Usage

``` r
frequency_delta_if_cleaned(tmat, thr, parent, spec, keep)
```

## Arguments

- tmat:

  transformed marker matrix for one sample

- thr:

  named threshold vector

- parent:

  logical parent-gate mask

- spec:

  population specification

- keep:

  logical, TRUE for events to retain

## Value

data.frame with the reported and cleaned percentage per population
