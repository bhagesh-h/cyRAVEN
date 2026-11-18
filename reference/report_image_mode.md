# How the report carries its figures: "embed" (default) or "link"

Set per run with `--report-images`, or for a session with
`options(cyRAVEN.report_images = "link")`. Embedding makes the report
one file that survives being moved; linking makes it a few megabytes
that opens at once but shows its figures only while it sits in the
results directory.

## Usage

``` r
report_image_mode()
```
