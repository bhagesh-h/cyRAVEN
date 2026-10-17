# Locate bead populations in one channel

Bead sets are discrete populations, so the peaks of the density are the
populations. Peaks are returned in ascending order, which is the order
assigned-value tables conventionally use.

## Usage

``` r
bead_peaks(
  x,
  n_expected,
  bins = 256L,
  smooth = 5,
  min_frac = 0.005,
  min_sep_log10 = 0.15
)
```

## Arguments

- x:

  channel values for the bead file, in linear instrument units

- n_expected:

  how many populations the bead set has

- bins, smooth:

  density resolution and smoothing

- min_frac:

  smallest share of events a population must hold

- min_sep_log10:

  minimum separation between two populations, in decades. Bead
  populations are narrow, so histogram noise readily splits one of them
  into two adjacent maxima; without a separation rule the split pair
  occupies two of the expected slots and the brightest real population
  is dropped.

## Value

numeric vector of peak locations, ascending, or NULL
