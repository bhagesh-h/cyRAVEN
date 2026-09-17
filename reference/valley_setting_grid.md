# The settings sweep that produces the method component

Three factors, each at the package default and one step either side.
Bins and smoothing set how much structure the histogram retains;
`min_rel_depth` sets how completely two modes must separate before the
gap between them is accepted as a threshold, which is the setting that
decides whether a minority positive population is found at all.

## Usage

``` r
valley_setting_grid()
```

## Value

a data.frame of setting combinations, one per row
