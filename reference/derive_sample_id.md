# Derive a sample label from the FCS keywords or filename WHY: wells rarely carry patient IDs; this gives a stable plot label when no sample map is supplied. \$WELLID is preferred, then \$SMNO, then the filename.

Derive a sample label from the FCS keywords or filename WHY: wells
rarely carry patient IDs; this gives a stable plot label when no sample
map is supplied. \$WELLID is preferred, then \$SMNO, then the filename.

## Usage

``` r
derive_sample_id(fname, kw = NULL)
```

## Arguments

- fname:

  The fname.

- kw:

  The kw.
