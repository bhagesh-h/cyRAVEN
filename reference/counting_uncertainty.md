# Counting uncertainty and detection limits for a population frequency

The uncertainty a frequency carries from the number of events behind it,
independent of where the thresholds were placed. Returned in percentage
points so it combines directly with the gate placement uncertainty.

## Usage

``` r
counting_uncertainty(k, n, z = 1, lod_events = 20L, loq_events = 50L)
```

## Arguments

- k:

  events in the population; may be a vector

- n:

  events in the parent gate the population is expressed against; scalar
  or the same length as `k`

- z:

  coverage factor for the Wilson half-width. The default of 1 gives a
  standard uncertainty, matching the convention used for the gate terms

- lod_events:

  events below which a population is not called detected

- loq_events:

  events below which a population is detected but not quantified

## Value

a data.frame with one row per element of `k`, carrying the counts, the
standard uncertainty, both limits expressed as percentages of the parent
gate, and a verdict
