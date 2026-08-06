# Read and validate the sample map

Schema (see README): file, well, sample_id, patient_id, timepoint,
is_control, panel. Only `file` is required. A supplied map that does not
cover every input file is a FATAL error: guessing well-\>patient
assignment from plate order would silently mislabel patients.

## Usage

``` r
load_sample_map(path, fcs_files)
```

## Arguments

- path:

  File path.

- fcs_files:

  The fcs files.
