# Population cuts to take from the other samples

Two cases, decided per panel and marker by how many samples placed their
cut at a density minimum (`valley_*`, `otsu`):

## Usage

``` r
peer_cut_replacements(thr, z_max = 3.5, min_peers = 3L, min_peer_frac = 0.5)
```

## Arguments

- thr:

  Data frame with `sample_id`, `panel`, `marker`, `threshold`, `source`,
  one row per sample and marker.

- z_max:

  Robust z beyond which a tail cut is replaced when the marker has no
  majority of density minima.

- min_peers:

  Minimum number of reference cuts.

- min_peer_frac:

  Minimum share of the panel's samples the reference must cover.

## Value

A data frame of replacements (`sample_id`, `panel`, `marker`,
`old_threshold`, `old_source`, `threshold`, `n_peers`, `robust_z`,
`reason`), or `NULL` when nothing is replaced.

## Details

- At least `min_peers` and at least `min_peer_frac` of the samples did.
  The marker is then bimodal in this assay, and a sample without a
  minimum – a `quantile_fallback` or a `tailgate*` cut – takes the
  median of the minima. On the KMT2C/KMT2D cohort one sample's CD19 tail
  cut at 3.24, against minima between 2.28 and 2.50 elsewhere, counted
  0.9% B cells where manual gating counted 5.7%.

- Fewer did. The other samples' cuts, fallbacks excluded, are the
  reference: a fallback takes their median, and a tail cut more than
  `z_max` robust SDs from them (MAD floored at 0.25, as for
  `scale_outlier`) does too.

Cuts from a control, a config entry or a manual override are never
touched.
