# Each sample's composition by exclusive label

Every CD45+ cell carries exactly one label from
[`score_populations()`](https://bhagesh-h.github.io/cyRAVEN/reference/score_populations.md),
the most specific population it satisfies or the catch-all remainder.
Counted per sample, those labels partition the parent, which is what a
compositional analysis needs; the declared populations' own gates
overlap and nest, and do not.

## Usage

``` r
exclusive_label_frequencies(pops, gates, verdicts = NULL, panel_of = NULL)
```

## Arguments

- pops:

  per-sample scored populations, each with `scored$labels`.

- gates:

  per-sample gate results, each with `masks$cd45_pos`.

- verdicts:

  per-sample staining verdicts, for `is_control` and `qc_status`.

- panel_of:

  named character, sample id to panel.

## Value

data.frame with `sample_id`, `panel`, `population` (the label), `count`,
`pct_of_cd45_pos`, `is_control` and `qc_status`, every label seen in the
run present for every sample (zero where absent), or NULL.
