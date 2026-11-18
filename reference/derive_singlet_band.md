# Derive the singlet gate from the FSC height-to-area ridge

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

  Width of the band below the ridge, in MADs. Default `3`.

## Details

ONE-SIDED. Two cells passing together give a pulse of about twice the
area and much less than twice the height, so a doublet's height-to-area
ratio is LOW (Wersto et al. 2001). Nothing about a doublet raises it.
The gate therefore has a lower bound only; an upper bound removed single
cells from the top of the ridge for no reason.

WHY NO SIDE-SCATTER STEP. SSC-H against SSC-A after the forward-scatter
gate is common manual practice (Cossarizza et al. 2019), and was tried
here with the same rule. Scored against manual gating on 21 whole-blood
files it removed 3.8% of the low-scatter cells and 0.7% of the
granulocytes: where little light is scattered the SSC ratio is noisy,
and a cut below its median removes lymphocytes rather than doublets.
Granulocyte shares rose by about a point against the manual counts, so
it is not applied.
