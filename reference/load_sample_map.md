# Read and validate the sample map

Schema (see README): file, well, sample_id, patient_id, timepoint,
is_control, panel, fmo_for, control_group. Only `file` is required.

## Usage

``` r
load_sample_map(path, fcs_files)
```

## Arguments

- path:

  File path.

- fcs_files:

  The fcs files.

## Details

`fmo_for` names the markers a file is the fluorescence-minus-one control
for, comma separated; `control_group` confines a control to the samples
sharing its value. Both are optional and inert when absent. See
[`parse_fmo_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/parse_fmo_map.md).
A supplied map that does not cover every input file is a FATAL error:
guessing well-\>patient assignment from plate order would silently
mislabel patients.
