# Write a sample-sheet template covering every input file

Emits EVERY reserved column, with the filename-derived identifiers
filled in and the rest blank, so the user edits a sheet that already
accounts for every file rather than assembling one and discovering an
omission at run time.

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

## Details

WHY EVERY COLUMN AND NOT A USEFUL SUBSET. The template used to carry
nine of the reserved columns. The rest are read by the pipeline exactly
the same way, but a user only learns they exist by reading the Inputs
chapter, so a run that could have reported cells per millilitre or used
an FMO control instead silently did neither. A blank column costs one
empty field and is deleted in a second; an absent column costs a re-run.
Blank columns are ignored, so a sheet returned untouched behaves exactly
as the nine-column one did.

`count.<Population>` columns are emitted only for populations the caller
names, because the column is keyed by a population that has to exist.
With none supplied one commented example header is written instead, so
the format is visible without inventing a population.
