---
name: cyraven
description: Run, configure and interpret cyRAVEN, the supervised flow-cytometry immunophenotyping pipeline that turns a directory of FCS files into per-sample gates, population frequencies, a shared UMAP and between-group statistics. Execute it through Docker by default. Use when the task involves running cyRAVEN, authoring its population specification YAML or sample map, choosing a transform, reading its output tables and diagnostic figures, or diagnosing a failed run. Triggers on cyRAVEN, cyraven.R, run_cyraven, FCS batch analysis, gating strategy YAML, thresholds_used.csv, population_frequencies.csv, staining_qc.csv, threshold_uncertainty.csv, quantile_fallback, detection limits, "why is my frequency table empty", flow cytometry pipeline.
---

# cyRAVEN

A declarative gating specification is applied to a batch of FCS files. Every
marker threshold is derived inside the sample it applies to, populations are
scored as Boolean conjunctions, and the result is tested against unsupervised
structure recovered from the same cells.

The design point: manual gating is the dominant source of technical variance in
multi-sample immunophenotyping. Operators gating identical files report
population sizes differing by about 32%, and analyst subjectivity accounts for up
to 78% of technical variability once more than one person is involved (Cadwell et
al., 2021, *PDA J Pharm Sci Technol* 75:33). Transferring fixed gate coordinates
between samples does not remove that variance, it converts it into a bias that
tracks staining intensity. cyRAVEN removes the analyst from threshold placement
and quantifies what remains.

Current version: 0.4.0.

## Execute through Docker

**Default to Docker for every run.** Do not fall back to a local R installation
unless the user states that Docker is unavailable, and say so explicitly when you
do. Two dependencies determine numerical results rather than convenience: `uwot`
is a stochastic embedder whose output varies with its own version and the BLAS
beneath it, and every gate is placed at a kernel density minimum, so a different
density implementation moves thresholds and therefore the reported frequencies.
The image pins R 4.4.3, a dated CRAN snapshot and Bioconductor 3.20.

Check for an existing image before building one:

```bash
docker images | grep cyraven
```

### Build

From the package root, the directory holding `DESCRIPTION`:

```bash
docker build -f inst/scripts/Dockerfile -t cyraven:0.4.0 .
```

15 to 25 minutes on a first build. It ends by running `--help` and printing every
package version, so a broken image fails at build time.

### Run

```bash
docker run --rm \
  -v "$PWD/data:/data:ro" \
  -v "$PWD/results:/results" \
  cyraven:0.4.0 \
  --dir /data/fcs \
  --sample-map /data/sample_map.csv \
  --config /data/panel.yaml \
  --group-column cohort --reference-group "Healthy controls" \
  --outdir /results
```

Every path inside a flag is a container path. `--outdir` must fall inside a
mounted volume or the output is lost when the container exits. On Windows
PowerShell use `${PWD}`.

### Demonstration data

When the user has no data to hand, or wants to see the output shape first:

```bash
mkdir -p demo results
docker run --rm -v "$PWD/demo:/demo" \
  --entrypoint Rscript cyraven:0.4.0 \
  /opt/cyraven/src/inst/scripts/demo_data.R /demo
```

This fetches two public CytoTrol acquisitions and partitions them into eight
pseudo-samples. **The `cohort` column it writes is randomised**: every
pseudo-sample comes from the same tube, so the true group difference is zero and
any significant result is a false positive. State that whenever you present
numbers from it. It is a calibration check, not a demonstration of sensitivity.

### Iterating on the code

Mount a checkout and set `CYRAVEN_SOURCE`; the entrypoint loads that tree through
`pkgload::load_all()` instead of the installed copy, so edits take effect without
a rebuild.

```bash
docker run --rm -v "$PWD:/src:ro" -v "$PWD/results:/results" \
  -e CYRAVEN_SOURCE=/src cyraven:0.4.0 --dir /data/fcs --outdir /results
```

Full build detail, resource tuning, Windows path handling and container-specific
failures are in `references/docker.md`.

## Start here

| The task | Go to |
|---|---|
| Run it on a batch of files | The commands above, then `references/docker.md` |
| Write the population YAML or sample map | `references/gating.md` §3 |
| Decide arcsinh against logicle | `references/gating.md` §4 |
| Read the outputs | `references/interpretation.md`, in its order |
| A run failed, or the numbers look wrong | `references/troubleshooting.md` |
| Choose between cyRAVEN and cyCONDOR | `references/interpretation.md` §9 |

## Inputs

Only the FCS files are required. Everything else improves what the run can
conclude.

1. **FCS files** in a directory. `--recursive` for subdirectories. Single-stain
   compensation controls are excluded by filename by default; they are instrument
   setup files with no CD45 population and form their own spurious panel group if
   admitted. `--exclude ''` keeps everything.

2. **Sample map** (`--sample-map`). Only `file` is required.

   ```csv
   file,sample_id,patient_id,timepoint,is_control
   donor01.fcs,D01,P01,baseline,FALSE
   unstained.fcs,UNS,,,TRUE
   ```

   `is_control = TRUE` marks an unstained control tube. Those samples are
   excluded from testing without being recorded as QC failures, and their 99.5th
   percentile can anchor a threshold for a marker with no resolvable density
   minimum. `--write-sample-map` emits a template.

3. **Population spec** (`--config`). See `references/gating.md` §3. Without it the
   built-in specification is used, which is a worked example for one
   myeloid-lymphoid panel and not a sensible default for other panels.

4. **Patient table** (`--patient-table`), optional. Covariates for confounder
   screening, and `wbc_per_ul` for absolute concentrations.

`--write-config` derives everything and writes the resulting YAML without running
the analysis, which shows what the pipeline detected before time is spent.

## Options that matter

Analyses that only add output are on by default. Analyses that change an existing
number are opt-in. That rule is deliberate and worth preserving in any change.

```bash
--transform logicle        # instead of arcsinh; changes every threshold
--cluster                  # SOM clustering, then cross-check against the spec
--explain-clusters         # learn gate geometry for undescribed clusters
--external-labels x.csv    # gate a labelling from another tool, e.g. cyCONDOR
--export-gates             # write Gating-ML 2.0 for the instrument
--write-baseline b.rds     # record where this run placed its thresholds
--baseline b.rds           # test this run against that record
--fail-on-drift            # make the verdict an exit code
--batch-column run_date    # quantify batch structure, per marker and overall
--correct-batch            # correct it, subject to the confounding refusal
--lod-events 20            # events below which a population is not detected
--loq-events 50            # events below which it is detected but not quantified
--drop-unstable-events     # exclude flagged acquisition intervals; changes counts
--subsample rare           # inverse-density draw; changes every UMAP
--calibration-beads b.fcs  # convert to MESF/ERF; changes every intensity
--no-uncertainty           # skip the gate uncertainty analysis
--no-acquisition-qc        # skip the acquisition stability check
--no-spreading             # skip the spillover spreading report
--no-miflowcyt             # skip the ISAC-structured report
--no-report                # skip report.html
--max-events-per-file N    # bound memory; also raises every detection limit
--include-qc-failed        # keep samples that failed staining QC
```

Three of these change numbers the previous run reported and default to off:
`--drop-unstable-events` changes every count in an affected file,
`--subsample rare` changes which cells are embedded, and `--calibration-beads`
changes the units every intensity is expressed in. Say so when you use them.

`--cluster` is the only output that can find a population the specification does
not describe. It is off by default because it costs runtime, not because it is
peripheral to the argument the package makes.

## Analytical constraints

These are positions the tool takes, not preferences. Overriding them produces
output that looks correct and is not.

**Thresholds are not transferred between samples.** A global override in the
config applies to every sample and reintroduces the fixed-coordinate bias the
pipeline exists to remove. When one sample's cut is wrong, the remedy is in
`references/troubleshooting.md`.

**`count` is an event count, not a cell number.** It scales with acquisition
duration and with `--max-events-per-file`, so it is not comparable between
samples. Report `pct_of_cd45_pos`, or `cells_per_ul` where a haemogram supplied
`wbc_per_ul`.

**Placement precision does not imply counting sufficiency.** A cut through a wide
empty gap is well determined however few events lie beyond it. Check `detection`
before testing a rare population; a test on a population below the limit of
quantification returns a statement about acquisition depth, and re-gating does
not change it.

**A specification is a hypothesis, not a validated result.** It is declared
before the data are examined, and six outputs exist to contradict it. Plausible
frequencies are not evidence that it held.

**Inspect the gates before quoting any number.** A threshold on a distribution
shoulder rather than a density minimum is visible in `gating_qc.png` and
`recon_diagnostics.png` and in no downstream table.

## Contributing

Changes to this package are additive: every file the previous version wrote is
still written, with the same name, the same columns in the same order, and the
same values. Verify that by running both versions on one cohort and comparing
byte for byte rather than by asserting it. `NEWS.md` documents the convention.

New entry points that consume random draws must save and restore `.Random.seed`.
`run_cyraven()` seeds once and the UMAP cell selection draws from that stream, so
a step that spends draws silently redraws every embedding in the run.

## Reference files

- `references/docker.md`: build and run in detail, resource tuning, Windows path
  handling, mounting, and container-specific failures
- `references/gating.md`: the four-gate hierarchy, per-sample thresholding and the
  `source` column, writing the population YAML, three-level markers, arcsinh
  against logicle, compensation
- `references/interpretation.md`: every output file and the order to read them in,
  with the columns that decide whether a result stands
- `references/troubleshooting.md`: failure modes with their causes
