# Resolve the input file list, and REPORT what was skipped

WHY the near-miss warning: real archives contain files whose names end
in something other than a clean ".fcs" – "sample.fcs copy" from a
filesystem duplication, "sample.FCS", "sample.fcs.bak". A strict pattern
silently drops them, and silently analysing 19 of 25 files is far worse
than failing: the output looks complete and nothing in it reveals the
missing cohort. So any file whose name CONTAINS ".fcs" but did not match
the pattern is listed as a warning with the exact flag needed to include
it.

## Usage

``` r
resolve_input_files(opt)
```

## Arguments

- opt:

  Named list of options. See build_option_list() for the full set.
