---
name: cyraven
description: Run, configure and interpret cyRAVEN — the supervised flow-cytometry immunophenotyping pipeline that turns a directory of FCS files into per-sample gates, population frequencies, a shared UMAP and between-group statistics. Use when the task involves running cyRAVEN (Docker or R), authoring its population specification YAML or sample map, choosing a transform, reading its output tables and diagnostic figures, or debugging a run. Triggers on cyRAVEN, cyraven.R, run_cyraven, FCS batch analysis, gating strategy YAML, thresholds_used.csv, population_frequencies.csv, staining_qc.csv, threshold_uncertainty.csv, quantile_fallback, "why is my frequency table empty", flow cytometry pipeline.
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

Current version: 0.2.0.

## Start here

| The task | Go to |
|---|---|
| Run it on a batch of files | Docker quickstart below, then `references/docker.md` |
| Write the population YAML or sample map | `references/gating.md` §3 |
| Decide arcsinh against logicle | `references/gating.md` §4 |
| Read the outputs | `references/interpretation.md` — follow its order, it is not arbitrary |
| A run failed, or the numbers look wrong | `references/troubleshooting.md` |
| Decide whether cyRAVEN or cyCONDOR answers the question | `references/interpretation.md` §7 |

## Docker is the recommended path

Use it unless there is a reason not to. Two of cyRAVEN's dependencies decide
numbers rather than convenience: `uwot` is a stochastic embedder whose output
changes with its own version and the BLAS beneath it, and every gate comes from a
kernel density estimate, so a different density implementation moves thresholds
and therefore the frequencies you would publish. The image pins R 4.4.3, a dated
CRAN snapshot and Bioconductor 3.20 for that reason.

### Build

From the package root, the directory holding `DESCRIPTION`:

```bash
docker build -f inst/scripts/Dockerfile -t cyraven:0.2.0 .
```

Allow 15 to 25 minutes for the first build; the dependency layer is compiled.
The build ends by running `--help` and printing every package version, so a
broken image fails at build time rather than an hour into an analysis.

### Run

Nothing is baked into the image — no FCS file, no patient table, no cohort name.
Mount data in, mount results out.

```bash
docker run --rm \
  -v "$PWD:/data:ro" \
  -v "$PWD/results:/results" \
  cyraven:0.2.0 \
  --dir /data --recursive \
  --config /data/panel.yaml \
  --sample-map /data/sample_map.csv \
  --group-column cohort --reference-group "Healthy controls" \
  --outdir /results
```

On Windows PowerShell use `${PWD}` rather than `$PWD`.

Paths inside the flags are container paths. `--config /data/panel.yaml` refers to
the mount, not the host.

### Editing code without rebuilding

Set `CYRAVEN_SOURCE` to a mounted checkout and the entrypoint loads that tree
through `pkgload::load_all()` instead of the installed copy:

```bash
docker run --rm -v "$PWD:/src:ro" -v "$PWD/results:/results" \
  -e CYRAVEN_SOURCE=/src cyraven:0.2.0 --dir /data ...
```

A rebuild is still needed once the source requires a package the image lacks.

### Without Docker

```r
install.packages("BiocManager")
BiocManager::install("flowCore")
remotes::install_github("bhagesh-h/cyRAVEN")
```

```r
library(cyRAVEN)
run_cyraven(list(
  dir = "data/", outdir = "results/", recursive = TRUE,
  sample_map = "data/sample_map.csv", config = "data/panel.yaml",
  group_column = "cohort", reference_group = "Healthy controls"
))
```

The CLI entry point is `system.file("scripts", "cyraven.R", package = "cyRAVEN")`.
`FlowSOM` is optional; without it, SOM clustering falls back to the in-package
implementation and nothing else changes.

## What you need before the first run

Only the FCS files are strictly required. Everything below improves what the run
can conclude.

1. **FCS files**, one per sample, in a directory. Use `--recursive` for
   subdirectories. Single-stain compensation controls are excluded by filename
   by default — they are instrument setup files, have no CD45 population, and
   form their own spurious panel group if admitted. `--exclude ''` keeps
   everything.

2. **Sample map** (`--sample-map`), a CSV. Only `file` is required.

   ```csv
   file,sample_id,patient_id,timepoint,is_control
   donor01.fcs,D01,P01,baseline,FALSE
   unstained.fcs,UNS,,,TRUE
   ```

   `is_control = TRUE` marks an unstained control tube. Those samples are
   excluded from testing without being recorded as QC failures, and their
   99.5th percentile can anchor a threshold for a marker that does not resolve a
   density minimum. Generate a template with `--write-sample-map`.

3. **Population spec** (`--config`), a YAML declaring each population as marker
   directions. See `references/gating.md` §3. Without it, the built-in
   specification is used, which is a worked example for one myeloid-lymphoid
   panel and not a sensible default for arbitrary panels.

4. **Patient table** (`--patient-table`), optional, keyed on patient identifier.
   Supplies covariates for confounder screening, and `wbc_per_ul` for absolute
   concentrations.

Run `--write-config` first to derive everything and write the resulting YAML
without doing the analysis. It shows what the pipeline detected before any time
is spent.

## The flags that matter

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
--batch-column run_date    # quantify batch structure
--correct-batch            # correct it, subject to the confounding refusal
--no-uncertainty           # skip the gate uncertainty analysis
--lod-events 20            # events below which a population is not detected
--loq-events 50            # events below which it is detected but not quantified
--max-events-per-file N    # bound memory; also raises every detection limit
--include-qc-failed        # keep samples that failed staining QC
--no-miflowcyt             # skip the ISAC-structured report
```

`--cluster` is the only output that can find a population the specification does
not describe. It is off by default because it costs runtime, not because it is
optional to the argument the package makes.

## Working rules

**Never transfer a threshold between samples to make numbers agree.** That
reintroduces the bias the pipeline exists to remove. If one sample's cut is
wrong, the answer is in `references/troubleshooting.md`, not a global override.

**Read `gating_qc.png` and `recon_diagnostics.png` before quoting any number
derived from them.** A threshold sitting on a shoulder rather than a minimum is
visible there and in no downstream table.

**Report `pct_of_cd45_pos`, never `count`.** `count` is an event count. It
depends on how long the tube ran and on `--max-events-per-file`, so it is not
comparable between samples and is not a cell number.

**Check `detection` before testing a rare population.** A well-placed cut says
nothing about whether enough cells were counted, and a test on a population below
the limit of quantification returns a statement about acquisition depth. Re-gating
does not fix it.

**Treat a specification as a hypothesis.** It is declared before the data are
examined, and five outputs exist to contradict it. A run where nothing
contradicts the spec is a result; a run where the contradictions were not read is
not.

**Changes to this package should be additive.** Every file the previous version
wrote should still be written, with the same name, the same columns in the same
order, and the same values. Verify it by running both versions on one cohort and
comparing byte for byte rather than by asserting it. `NEWS.md` documents this
convention and the RNG hazard that makes it easy to break: `run_cyraven()` seeds
once and the UMAP cell selection draws from that stream, so any new step that
consumes draws silently redraws every embedding in the run. Save and restore
`.Random.seed` in any new entry point.

## Reference files

- `references/docker.md` — build and run in detail, resource tuning, Windows and
  PowerShell specifics, mounting, and the Docker-specific failure modes
- `references/gating.md` — the four-gate hierarchy, per-sample thresholding and
  the `source` column, writing the population YAML, three-level markers,
  arcsinh against logicle, compensation
- `references/interpretation.md` — every output file, and the order to read them
  in, with the specific columns that decide whether a result stands
- `references/troubleshooting.md` — failure modes with their actual causes
