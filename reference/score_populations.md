# Evaluate the population spec against one file's thresholded markers

Evaluate the population spec against one file's thresholded markers

## Usage

``` r
score_populations(
  tmat,
  thr,
  parent,
  spec = default_population_spec(),
  hi_thr = NULL,
  quiet = FALSE
)
```

## Arguments

- tmat:

  matrix of transformed marker values (cells x markers, named cols)

- thr:

  named numeric vector of thresholds

- parent:

  logical mask (the CD45+ parent)

- spec:

  Population specification mapping population name to marker directions.
  See
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).
  Default
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).

- hi_thr:

  named numeric: upper bound for "intermediate" requirements. Derived by
  derive_intermediate_bounds(); a marker requested as "intermediate"
  without an upper bound makes its population UNAVAILABLE rather than
  silently collapsing to "above".

- quiet:

  Suppress the per-population NOTEs about `any_of` groups losing a
  member. The gate-uncertainty bootstrap re-scores every population once
  per replicate, and the note is a property of the SPECIFICATION against
  the PANEL – identical on every replicate, because neither changes.
  Emitted from there it produced 2,120 of 2,399 log lines on a 20-sample
  run, 88% of the log saying one of two things, which buries the NOTEs
  that do differ between samples. Callers that score once leave it FALSE
  and see it once per sample; callers that score in a loop pass TRUE.

## Value

list(masks, labels, unavailable)
