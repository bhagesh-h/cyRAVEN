# Run explore mode

Unsupervised discovery over every eligible channel, writing into
`<outdir>/explore/`. Every argument beyond `reads`, `fpr`, `opt` and
`outdir` is optional: supplied, it makes the output better; absent, the
run degrades to what a standalone clusterer would produce and records
that it did so in `explore_provenance.csv`.

## Usage

``` r
run_explore(
  reads,
  fpr,
  opt,
  outdir,
  transforms = NULL,
  gates = NULL,
  verdicts = NULL,
  pops = NULL,
  group_of = NULL,
  confounding = NULL,
  file_paths = character(0),
  total_counts = NULL
)
```

## Arguments

- reads:

  Named list of reads.

- fpr:

  Panel resolution.

- opt:

  Option list.

- outdir:

  Run output directory.

- transforms:

  Per-panel transforms. Derived here if absent.

- gates:

  Per-sample gate objects, for thresholds and the viability marker.

- verdicts:

  Per-sample staining QC verdicts. Recorded, NOT used to exclude:
  explore does not depend on the CD45 gate, so a sample the declared
  path had to drop still contributes here. That is a capability rather
  than an oversight, which is why it becomes a column instead of a
  filter.

- pops:

  Per-sample scored populations, for the bridge tables.

- group_of:

  Named character vector, sample id -\> group.

- confounding:

  One-row data frame from the batch/group confounding check.

- file_paths:

  Full paths of the FCS files the run was given. Needed only when the
  declared pipeline has already freed the event matrices, because
  `reads[[s]]` records a basename and cannot be re-opened on its own.

## Value

Invisibly, `list(dir, files, gaps)`. `gaps` is what explore has to say
about the declared specification – populations spanning several
clusters, and clusters no population covers. The caller writes it to
`spec_gaps.csv` only under `--maybe-learn`; this function never writes
outside `explore/`.
