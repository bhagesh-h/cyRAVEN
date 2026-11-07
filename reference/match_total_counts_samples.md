# Match parsed total counts to this run's sample_id, on patient AND timepoint

The timepoint half of the key is what separates this from
[`match_absolute_counts_samples()`](https://bhagesh-h.github.io/cyRAVEN/reference/match_absolute_counts_samples.md).
A yield belongs to one draw, so HS5_018 at d0 and HS5_018 at d7 are
different measurements and must not be allowed to collapse onto each
other. Where the sheet carries no timepoint row the join falls back to
patient alone, and that is accepted only while it stays unambiguous: a
patient with several acquisitions and one unlabelled yield is reported
as unmatched rather than broadcast across their timepoints.

## Usage

``` r
match_total_counts_samples(tc, smap)
```

## Arguments

- tc:

  Parsed total counts.

- smap:

  Sample map joining sample_id to patient_id.
