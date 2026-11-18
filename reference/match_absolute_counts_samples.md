# Match each –absolute-counts row to this run's sample_id

Tries patient_id first, then the acquisition filename (both
case/whitespace-insensitive, extension and a trailing " copy" stripped),
because real exports are not internally consistent about which one
labels a given row (seen in practice: most rows use the filename, a
handful use the patient_id instead). Whichever key matches wins; rows
matching neither are dropped with a NOTE naming them, rather than
silently guessed at.

## Usage

``` r
match_absolute_counts_samples(ac, smap)
```

## Arguments

- ac:

  The ac.

- smap:

  Sample map joining sample_id to patient_id.
