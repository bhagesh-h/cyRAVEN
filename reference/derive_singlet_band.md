# Derive the singlet band from the FSC-H/FSC-A ridge

WHY: doublet exclusion assumes height and area track together for single
cells. The RATIO at which they track is instrument- and gain-dependent –
on this cytometer the ridge sits near 0.56, not the 1.0 that a hardcoded
gate would assume, and a band centred on 1.0 discards essentially every
real cell. Centering on the measured median with a MAD-scaled width
adapts to any instrument. MAD is used rather than SD because doublets
are a heavy tail.

## Usage

``` r
derive_singlet_band(ex, sc, parent_mask, k = 3)
```

## Arguments

- ex:

  Numeric expression matrix, events x channels.

- sc:

  The sc.

- parent_mask:

  The parent mask.

- k:

  Neighbourhood size, or number of clusters, depending on the function.
  Default `3`.
