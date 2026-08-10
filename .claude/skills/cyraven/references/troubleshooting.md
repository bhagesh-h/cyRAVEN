# Troubleshooting

Symptom, cause, fix. Everything here is a failure mode the code raises or guards
against, not a hypothetical.

## Input discovery

### Fewer samples analysed than exist in the directory

Silently analysing 19 of 25 files is worse than failing, because the output looks
complete. The pipeline guards this: any file whose name **contains** `.fcs` but
did not match `--pattern` is listed as a warning with the flag needed to include
it. Read the log.

Real archives contain `sample.fcs copy` from filesystem duplication, `sample.FCS`,
`sample.fcs.bak`. The default pattern is `[.]fcs$`.

```bash
--pattern '[.]fcs( copy)?$'
```

### "no FCS files found, pass --files or --dir"

If `--dir` has subdirectories, the log says so. Add `--recursive`.

### Compensation controls appear as their own panel and fail QC

They are excluded by filename by default, matching
`compensation|single[ _-]?stain|comp[ _-]?control`. If yours are named otherwise
they will be admitted, and each contains beads or cells carrying one fluorophore:
no CD45 population, no viable-cell population, no marker set in common with the
panel. They form a spurious panel group and clutter every figure.

```bash
--exclude 'compensation|single_stain|my_naming_convention'
--exclude ''      # disables the default entirely
```

## Empty or nearly empty results

### The frequency table is empty, or every population is UNAVAILABLE

Almost always a nomenclature mismatch between the specification and `$PnS`, not a
biologically absent population.

Check what the panel actually resolved:

```r
rd <- read_fcs_resolved("data/sample01.fcs")
names(rd$marker_cols)
```

or run once and read `thresholds_used.csv`, which lists every resolved marker.
Then make the `populations:` keys match exactly.

The run log names each population whose markers are not all present. It is
reported UNAVAILABLE and scored as absent rather than silently dropped.

### Every sample fails staining QC with "no separable CD45+ mode"

Check in this order:

1. **Is CD45 in the panel at all?** Without it, all viable events become the
   parent and a warning is raised. Percent-of-leukocytes and percent-of-viable
   are not interchangeable denominators.
2. **Look at `gating_qc.png`.** If the CD45 density is unimodal, no valley
   exists. Declare an unstained control through `is_control` in the sample map,
   or set the threshold explicitly in the config.
3. **Check the transform.** A cofactor that is wildly wrong compresses everything
   into one mode. `--transform logicle` resolves markers that arcsinh does not.
4. **Lower the floor** with `--min-cd45-pct` if the biology genuinely has few
   leukocytes, or use `--include-qc-failed` to inspect. Forced-in samples carry
   no evidence in their percentages; `staining_qc.csv` records which were forced.

If every marker also came back `Inf` or `NaN`, the cofactor was zero. That
specific collision between `--cofactor` and `--cofactor-from-first-sample`
through R's partial matching of `$` on lists is fixed, but the symptom is worth
recognising: `asinh(x/0) = Inf`, every threshold `Inf`, every sample failing with
a message that points nowhere near the cause.

## Thresholds

### Most thresholds are `quantile_fallback`

The marker did not separate positive from negative. In order:

1. **Try `--transform logicle`.** For dim markers without clear bimodality the
   transform decides whether a stable threshold exists at all. On a seven-colour
   T cell panel, CCR7 resolved to a valley in some samples and a fallback in
   others under arcsinh, and to a valley in every sample under logicle.
2. **Supply an unstained control** (`is_control = TRUE`). Unimodal markers get
   `control_q995` instead of a guess.
3. **Accept it.** A marker falling back across most samples is not resolving in
   that panel. No gating strategy recovers it. The fix is panel design.

Frequencies from a fallback are not invalid. They carry no evidence of
separation, which is a different statement, and `source` is how the reader knows.

### A threshold is flagged in `threshold_scale_qc.csv`

Robust *z* against the panel-wide median. Thresholds are derived per sample **by
design**, so some spread is expected. Flagged rows are candidates for review, not
errors.

Do not fix this by pinning the threshold in the config. That applies to every
sample and reintroduces the fixed-coordinate bias the pipeline exists to remove.
Look at that sample's density in `gating_qc.png` first.

### `bootstrap_valley_rate` is low

Below about 0.8 the cut is unresolved, not merely imprecise. It was found by
histogram noise that happened to clear the depth rule. It looks identical to a
solid threshold in `thresholds_used.csv`. Treat any frequency reading that marker
as provisional.

## The embedding

### The UMAP is uninformative: one blob, or structure that tracks sample rather than biology

Feature selection is the usual cause. Height and width channels are redundant
duplicates of area and would double-weight every marker; scatter and Time encode
technical variation. `select_umap_features()` drops
`^FSC|^SSC|^Time|Time Stamp|^Event|Width$|Height$|-H$|-W$` and prefers a
lineage list.

Watch the log for:

```
[features] fewer than 2 preferred lineage markers present;
           falling back to all N eligible markers
```

That fallback is a warning, not a solution. Name the panel's lineage markers:

```bash
--umap-markers CD3,CD4,CD8,CD14,CD16,CD19,CD56
--umap-markers-all     # every eligible marker instead
```

If structure tracks acquisition order rather than group, quantify it before
correcting it: `--batch-column`, then read §6.2 of `interpretation.md`.

### Adding samples moved every coordinate

Expected. A fresh embedding is fitted to whatever cells it is given.

```bash
--save-umap-model m.rds    # on the first run
--umap-model m.rds         # project later batches into the same space
```

Right for adding more of the same kind of sample. Wrong for a new panel or a new
cell type; retrain.

### The UMAP changed between runs with the same seed

Any new step that consumes random draws before the embedding changes which cells
are selected. `run_cyraven()` seeds once and the UMAP cell selection draws from
that one stream. Every entry point that needs randomness must save and restore
`.Random.seed`: see `uncertainty.R` and `clustering.R` for the pattern.

If this happens after a code change, that is the cause.

## Batch correction

### "batch correction refused"

Cramér's *V* between batch and group exceeded `--batch-max-cramers-v` (0.6). At
that level of confounding, removing the batch effect and removing the biological
effect are the same operation.

This is the most informative result the diagnostic can return. `--force-batch-correction`
overrides and the override is recorded in the manifest, but the finding it
produces cannot be attributed to biology.

## External labels

### "matched fewer than 100 cells" or the join produced nothing

The join is on **sample id and event index**, never row position. cyCONDOR
subsamples at `prep_fcd(max_cell=)` and cyRAVEN takes its own cap for the
embedding, so the two tables are not row-aligned and a positional join would
relabel every cell without raising anything.

Check that the label CSV carries an event index that refers to the original file,
and that sample identifiers match the sample map. `join_external_labels()` reports
the matched fraction.

### Leave-one-donor-out is taking hours

Every fold refits the whole strategy, so cost is one fit per donor per label.

```bash
--transfer-max-donors 4
--transfer-max-cells 10000
--external-max-labels 3
```

Both caps are logged when they bind, and the summary then describes the donors
scored rather than the cohort.

## Memory and runtime

### The process was killed with no R error

The OOM reaper. Peak memory is dominated by the per-sample expression matrices
and the embedding.

```bash
--max-events-per-file 300000   # bounds per-file memory
--max-cells 100000             # bounds the embedding
--no-session                   # skip session_state.RData
```

`--keep-exprs` retains raw matrices for the whole run and should be off unless
you need them.

Note that `--max-events-per-file` lowers every parent-gate event count, which
raises the detection limit of every rare population proportionally. It samples
evenly through the acquisition rather than taking the first N, so it does not
bias toward the start of the tube.

### The run is far slower than expected

If in Docker, check for an orphaned container competing for CPU. Killing the
shell that started `docker run` does not stop the container.

```bash
docker ps
docker kill <container-id>
```

Otherwise `--threads`, or `-e OMP_NUM_THREADS=N` in the container, which defaults
to 4.

## Docker

See `docker.md` for the full list. The two that catch everyone:

**`Rscript: not a valid option`**: the image has an `ENTRYPOINT`, so your command
arrives as arguments to it. Use `--entrypoint Rscript` or `--entrypoint bash`.

**Results directory empty after a successful run**: `--outdir` pointed outside
every `-v` mount. `/results` is a declared `VOLUME`, so an unmounted path writes
into an anonymous volume that vanishes with `--rm`.

## When the numbers look wrong but nothing errored

Work through `interpretation.md` in order. The common cause is a stage that was
never read:

- Gates never inspected (`gating_qc.png`)
- A marker flagged in `threshold_drift_stats.csv`, making an abundance difference
  partly definitional
- `difference_over_gate_u` below 1, meaning the groups differ by less than the
  cut's own movement
- A `raw_only__possible_composition_artefact` verdict in
  `compositional_concordance.csv`
- `f1_median` read instead of `f1_min` in `gate_transferability_summary.csv`
