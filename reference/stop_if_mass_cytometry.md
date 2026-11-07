# Stop on mass cytometry data, which this package does not process

A CyTOF file has no scatter, carries ion counts rather than
fluorescence, and needs bead normalisation, a DNA-intercalator singlet
gate and an arcsinh cofactor of 5 before analysis (Finck et al. 2013;
Nowicka et al. 2017). None of that is implemented here, and the
fluorescence path would run on it far enough to produce numbers before
failing at the scatter gate. Recognised by FCS 3.2 `$PnTYPE` "Mass", by
`$CYT`, or by channel names written as a metal isotope with the
`Di`/`Dd` suffix the acquisition software uses.

## Usage

``` r
stop_if_mass_cytometry(nm, desc, kw, label = "")
```

## Arguments

- nm, desc:

  `$PnN` and `$PnS` values.

- kw:

  keyword list.

- label:

  file name for the message.

## Value

`invisible(FALSE)`, or an error.
