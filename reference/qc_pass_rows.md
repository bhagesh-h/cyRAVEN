# Rows a group comparison may use: real, QC-passing samples only

Drops unstained controls and staining-QC failures from a table before it
reaches a figure or a test.

## Usage

``` r
qc_pass_rows(d)
```

## Arguments

- d:

  A table carrying `is_control` and/or `qc_status`, whichever are
  present. Tables that carry neither (for example derived population
  ratios, which are already built only from passing rows) are returned
  unchanged.

## Value

`d` with disqualified rows removed.

## Details

Both checks are needed, and `is_control` alone is not enough. Population
scoring runs for every gated sample regardless of its staining verdict –
the frequency, MFI and functional tables deliberately keep failed and
control rows so the CSV exports stay auditable. A sample that FAILED
staining (no separable CD45+ mode, or below the minimum CD45 percentage)
has `is_control = FALSE` and `qc_status = "failed"`; its percentages
carry no more evidence than a control's, and without the second check
they would enter every abundance figure and test by default rather than
only under `--include-qc-failed`. That flag works by relabelling a
forced sample's `qc_status` to `"pass"`, so filtering on `qc_status`
here is what gives the flag any meaning for figures and tests rather
than for the embedding alone.
