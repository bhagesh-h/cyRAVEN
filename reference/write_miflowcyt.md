# Write a MIFlowCyt-structured report of the run

Emits the parts of the ISAC MIFlowCyt checklist a run can establish on
its own, and marks the parts that require a human as outstanding rather
than leaving them out.

## Usage

``` r
write_miflowcyt(
  path,
  reads,
  opt = NULL,
  spec = NULL,
  transforms = NULL,
  fpr = NULL,
  thresholds = NULL
)
```

## Arguments

- path:

  output file

- reads:

  the per-sample list returned by the reading stage

- opt:

  parsed options

- spec:

  population specification, as scored

- transforms:

  named list of per-panel transform objects

- fpr:

  panel fingerprint result, for the panel-to-sample assignment

- thresholds:

  the `thresholds_used.csv` table, or NULL

## Value

the path, invisibly

## Details

The instrument and data-file sections are read from the FCS keyword
block. The data-analysis section is read from the run itself, so it
cannot disagree with what was actually computed.
