# The CD45 half of [`staining_verdict()`](https://bhagesh-h.github.io/cyRAVEN/reference/staining_verdict.md)

The CD45 half of
[`staining_verdict()`](https://bhagesh-h.github.io/cyRAVEN/reference/staining_verdict.md)

## Usage

``` r
cd45_verdict(
  gate,
  declared_control = NA,
  min_cd45_pct = 5,
  force_include = FALSE
)
```

## Arguments

- gate:

  The gate.

- declared_control:

  The declared control. Default `NA`.

- min_cd45_pct:

  The min cd45 pct. Default `5`.

- force_include:

  when TRUE (–include-qc-failed), a DECLARED sample (not a control) that
  would otherwise fail is included anyway, using whatever threshold was
  already computed. It is still excluded when the sample map itself says
  this is a control tube – that case is an unstained reference, and
  embedding it puts a cloud of background-only events in the shared
  space regardless of this flag. Only meant for the "too few samples to
  afford losing any" situation; the verdict text still records the
  forcing so it stays auditable in staining_qc.csv.

## Value

as
[`staining_verdict()`](https://bhagesh-h.github.io/cyRAVEN/reference/staining_verdict.md),
without `pct_live_of_singlets`.
