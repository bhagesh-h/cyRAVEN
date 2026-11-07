# Derive one set of transform parameters for a panel

Pooled across up to `max_samples` files so that every sample in a panel
is transformed on the same scale. See the note at the top of this file
for why that matters more than fitting each file well.

## Usage

``` r
derive_logicle_pooled(
  reads,
  sids,
  markers = NULL,
  m = 4.5,
  a = 0,
  max_samples = 8L,
  max_events = 50000L,
  seed = 42L
)
```

## Arguments

- reads:

  list of read_fcs_resolved() results

- sids:

  sample ids belonging to one panel

- markers:

  marker names to parameterise

- m, a:

  logicle display parameters

- max_samples:

  files pooled

- max_events:

  events drawn per file

- seed:

  RNG seed; the stream is restored on exit

## Value

named list of per-marker logicle parameters
