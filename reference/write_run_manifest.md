# Write the run manifest

Write the run manifest

## Usage

``` r
write_run_manifest(
  path,
  opt = NULL,
  files = NULL,
  extra = NULL,
  status = c("running", "completed", "failed"),
  started = NULL
)
```

## Arguments

- path:

  output file

- opt:

  parsed options

- files:

  input FCS paths (may be NULL on the first, pre-resolution write)

- extra:

  named list of run-derived values to record (cofactors, seeds,
  embedding parameters, sample counts)

- status:

  "running" \| "completed" \| "failed"

- started:

  Run start time.
