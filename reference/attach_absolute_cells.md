# Attach absolute cell numbers to a per-sample abundance table

Multiplies each population share by that acquisition's total yield. Rows
for samples with no total are left NA rather than dropped, so the
relative table and the absolute table always have the same shape and the
same row order and can be read side by side.

## Usage

``` r
attach_absolute_cells(ab, tc, share_col = "pct_of_gated")
```

## Arguments

- ab:

  Abundance table with sample_id and a share column.

- tc:

  Matched total counts from
  [`load_total_counts()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_total_counts.md).

- share_col:

  Column holding the share, in percent.

## Value

`ab` with total_cells, cells_absolute and count_basis added.

## Details

WHY `cells_absolute` AND NOT A REPLACEMENT FOR pct_of_gated: the two
answer different questions and have different error. The share is
measured once, by this pipeline, on events it counted itself. The
product carries the yield's error as well, from an instrument this
pipeline never saw. Overwriting the share would hide that, so both
travel together and the provenance column says so.
