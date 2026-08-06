# Run the whole uncertainty analysis over a scored cohort

Driver called by
[`run_cyraven()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_cyraven.md)
after populations are scored. It touches nothing the pipeline has
already computed.

## Usage

``` r
run_gate_uncertainty(
  pops,
  gates,
  verdicts,
  panel_of,
  spec,
  B = 100L,
  seed = 42L,
  max_events = 20000L
)
```

## Arguments

- pops:

  the per-sample scoring results

- gates:

  the per-sample gate objects

- verdicts:

  the per-sample staining verdicts

- panel_of:

  named character vector mapping sample to panel

- spec:

  population specification

- B:

  bootstrap replicates

- seed:

  base seed

- max_events:

  cap on events per replicate

## Value

list(thresholds, frequencies, budget), each a data.frame or NULL
