# Look up a per-sample, per-marker threshold override

Reads the `sample_overrides:` block of the config, which is keyed on
sample identifier and then on marker:

## Usage

``` r
sample_override(overrides, sample_id, marker)
```

## Arguments

- overrides:

  the `sample_overrides` list, or NULL

- sample_id:

  sample identifier

- marker:

  marker name

## Value

a list carrying at least `threshold`, or NULL when none applies

## Details

    sample_overrides:
      D07:
        CCR7:
          threshold: 2.15
          reason: "valley found inside the negative mode, see gating_qc.png"
          set_by: "initials or name"

Sample and marker names are matched exactly, as they appear in
`thresholds_used.csv`. A block naming a sample or marker that does not
exist is inert rather than an error, because a config is routinely
reused across cohorts that do not all contain the same tubes; unmatched
entries are reported once by
[`report_unused_overrides()`](https://bhagesh-h.github.io/cyRAVEN/reference/report_unused_overrides.md)
rather than failing the run.
