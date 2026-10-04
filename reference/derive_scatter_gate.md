# Derive the leukocyte scatter gate from the FSC-A density

WHAT: lower FSC bound = the deepest valley in log10 FSC-A; upper bound =
a high percentile; SSC bounds = robust percentiles within the FSC
window.

## Usage

``` r
derive_scatter_gate(ex, sc, fsc_hi_q = 0.999, ssc_q = c(0.005, 0.995))
```

## Arguments

- ex:

  Numeric expression matrix, events x channels.

- sc:

  The sc.

- fsc_hi_q:

  The fsc hi q. Default `0.999`.

- ssc_q:

  The ssc q. Default `c(0.005, 0.995)`.

## Details

WHY THIS IS COUNTER-INTUITIVE AND IMPORTANT: on this data a large
majority of recorded events (roughly two thirds in the test batch) are
sub-cellular debris sitting BELOW the FSC valley, and that debris is
CD45-negative. A conventional "draw a box low on FSC/SSC for
lymphocytes" gate therefore selects debris, not cells. Deriving the
boundary from the valley puts the gate above the debris mode wherever
that mode happens to fall for a given instrument and sample.
