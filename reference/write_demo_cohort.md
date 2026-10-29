# Write the demonstration cohort

Writes 35 FCS files from the graft-versus-host disease cohort that ships
with flowCore, together with `samples.csv`, `panel.yaml` and
`sample_map.csv`, so a first run needs no data of your own and no
download. This is the same script the documentation has always used; the
flag only saves typing its path inside the image.

## Usage

``` r
write_demo_cohort(dir)
```

## Arguments

- dir:

  directory to write into; created if absent.

## Value

invisibly, `dir`.
