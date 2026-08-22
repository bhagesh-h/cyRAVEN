# Write a single self-contained HTML report of everything a run produced

Presents the outputs in the order the documentation says they must be
read, because each stage can invalidate the ones after it. Sections
whose files are absent are omitted rather than shown empty.

## Usage

``` r
write_run_report(outdir, opt = NULL, verdicts = NULL, failure = NULL)
```

## Arguments

- outdir:

  the results directory, which is also where the report is written

- opt:

  parsed options, used for the header

- verdicts:

  per-sample staining verdicts, used for the summary banner

- failure:

  a condition or message; when given, the report is written as the
  record of a failed run rather than a completed one.

## Value

the path, invisibly, or NULL when there is nothing to report

## Details

The file references nothing: every figure is embedded at full resolution
as a data URI and every table as JSON, so the report can be moved,
emailed or archived on its own. Figures can be zoomed and downloaded at
full resolution; tables can be searched, sorted, paged at 10/50/100/all
rows and exported to CSV.

When `failure` is supplied the same report is written for a run that
stopped early, with a diagnosis section first and every output produced
up to the failure embedded below it. The partial output is usually where
the evidence is, so it is kept rather than discarded.
