# Write a sample-sheet template covering every input file

Emits the reserved columns with the filename-derived identifiers filled
in and the rest blank, so the user edits a sheet that already accounts
for every file rather than assembling one and discovering an omission at
run time.

## Usage

``` r
write_samplesheet_template(
  fcs_files,
  path,
  sample_ids = NULL,
  populations = character(0)
)
```

## Arguments

- fcs_files:

  The fcs files.

- path:

  File path to write.

- sample_ids:

  Optional identifiers; derived from the filenames otherwise.

- populations:

  Optional population names; one `count.` column each.
