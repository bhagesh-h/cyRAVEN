# Normalize any –absolute-counts input down to one CSV on disk

WHY funnel every format through a written CSV rather than branch the
parser on format: the source layout is irregular (blank separator rows,
a units-label spacer column) and that scan should exist exactly once,
regardless of whether the input arrived as .xlsx, .csv, or .tsv. CSV,
not TSV, because every other tabular input this pipeline reads
(–sample-map, –patient-table) is already CSV – one delimiter convention
throughout rather than a second one just for this file. Writing the
intermediate into outdir (not a tempfile) also gives the user something
to open when a match or a number looks wrong.

## Usage

``` r
absolute_counts_to_csv(path, outdir)
```

## Arguments

- path:

  File path.

- outdir:

  Directory to write outputs to.
