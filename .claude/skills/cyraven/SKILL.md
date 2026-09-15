---
name: cyraven
description: Run, configure and interpret cyRAVEN, the supervised flow-cytometry immunophenotyping pipeline that turns a directory of FCS files into per-sample gates, population frequencies, a shared UMAP and between-group statistics. Execute it through Docker by default. Use when the task involves running cyRAVEN, authoring its population specification YAML or sample map, choosing a transform, reading its output tables and diagnostic figures, or diagnosing a failed run. Triggers on cyRAVEN, cyraven.R, run_cyraven, FCS batch analysis, gating strategy YAML, thresholds_used.csv, population_frequencies.csv, staining_qc.csv, threshold_uncertainty.csv, gate_adjustments.csv, total_counts.csv, explore_cluster_identity.csv, explore_cluster_subsets.csv, spec_gaps.csv, suggested_config_next_run.yaml, quantile_fallback, --auto-fix-gates, --total-counts, --split-by-timepoint, --explore, --maybe-learn, absolute cell counts, compositional constraint, detection limits, "why is my frequency table empty", "why is every gate exactly 90%", flow cytometry pipeline.
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

Current version: 1.0.0, and the source tree is **ahead of it**. The published
image installs 1.0.0; `NEWS.md` carries a development section on top of that with
`--total-counts`, `--auto-fix-gates`, `--split-by-timepoint`, explore cluster
identity and the immune-subset annotation. Those reach a run only through a
mounted checkout (`CYRAVEN_SOURCE`, see "Iterating on the code"), so a flag
documented here can be rejected by a plain `cyraven:1.0.0` run. When a flag is
refused as `long flag "x" is invalid`, that is the cause, and it is also the
signature of the reverse case -- an old command line carrying a flag the source
has since removed. `--add-subsets` is the one that has actually been removed; see
Explore mode.

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

Prints one line per image already present, nothing if none is:

```
bhagesh/cyraven   1.0.0   8261eb7fc1e0   2.59GB
cyraven           1.0.0   8261eb7fc1e0   2.59GB
```

### Get the image

**Default to pulling.** It is a download rather than a 20-minute compile, and it
is the image the documented runs were produced with.

```bash
docker pull bhagesh/cyraven:1.0.0
```

```bash
docker tag bhagesh/cyraven:1.0.0 cyraven:1.0.0
```

The `docker tag` is only so the commands below can say `cyraven:1.0.0`.

**Build instead** when the user has edited `R/`, the Dockerfile or the config
templates, when the host is air-gapped, or when they say they would rather not
run a third-party image. A pull cannot contain local changes, so if the user is
iterating on the source, building is the only correct option. The Dockerfile
refers to paths inside the package, so the build context is the repository root,
the directory holding `DESCRIPTION`:

```bash
git clone https://github.com/bhagesh-h/cyRAVEN.git
```

```bash
cd cyRAVEN
```

```bash
docker build -f inst/scripts/Dockerfile -t cyraven:1.0.0 .
```

15 to 25 minutes on a first build. It ends by running `--help` and printing every
package version, so a broken image fails at build time.

Both routes pin R 4.4.3, the same dated CRAN snapshot and Bioconductor 3.20, so
they are numerically interchangeable. Check what is present before doing either:

```bash
docker images | grep cyraven
```

Prints one line per image already present, nothing if none is:

```
bhagesh/cyraven   1.0.0   8261eb7fc1e0   2.59GB
cyraven           1.0.0   8261eb7fc1e0   2.59GB
```

### The two inputs

A run takes the FCS directory, **one CSV** and **one YAML**.

- `--samples samples.csv`, one row per FCS file: what it is (`file`,
  `sample_id`, `is_control`, `panel`, `fmo_for`), whose it is (`patient_id`,
  `cohort`, `sex`, `age_years`, ...), any study variable under its own name, and
  externally measured counts in `count.<Population>` columns.
- `--config analysis.yaml`, the analysis: `populations:`,
  `functional_blocks:`, `ratios:`, `sample_overrides:`, `colors:`, `metadata:`,
  and a `samples:` section naming which column is the group and which the batch.

Anything varying per sample goes in the CSV; anything that is one decision for
the whole study goes in the YAML.

`--sample-map`, `--patient-table` and `--absolute-counts` still work and are not
deprecated, but **cannot be combined with `--samples`**. Prefer `--samples` for
new work; leave existing three-file invocations alone.

### Run

Never hand-author the sheet, and never run before checking. This sequence:

**1. Write a sheet with a row for every file.**

```bash
docker run --rm -v "$PWD/data:/data" cyraven:1.0.0 \
  --dir /data/fcs --recursive --write-samples /data/samples.csv
```

**2. The user fills in `patient_id`, `cohort` and any study variable.** An edit,
not a command.

**3. Validate against the FCS headers.** Seconds, and writes nothing.

```bash
docker run --rm -v "$PWD/data:/data:ro" -v "$PWD/results:/results" \
  cyraven:1.0.0 --dir /data/fcs --recursive \
  --samples /data/samples.csv --config /data/analysis.yaml \
  --outdir /results --check
```

**4. Run it.**

```bash
docker run --rm \
  -v "$PWD/data:/data:ro" \
  -v "$PWD/results:/results" \
  cyraven:1.0.0 \
  --dir /data/fcs --recursive \
  --samples /data/samples.csv \
  --config /data/analysis.yaml \
  --group-column cohort --reference-group "Healthy controls" \
  --batch-column acquisition_date --cluster \
  --no-session \
  --outdir /results
```

**Always run `--check` first and show the user its output.** It reports the
markers resolved from the files, any specification entry matching none of them,
whether the sheet covers every file, the group levels and their sizes, and the
available study variables. A marker-name mismatch caught there costs a second;
caught during a run it costs the run. It is the single highest-value habit when
driving this tool.

**Run `--list-channels` before writing the config.** It prints what each channel
resolves to, which is the name a population definition has to use, and on a
spectral panel that differs from `$PnS`: the file reads `CD45 : SparkUV-387 -
Area` and the run uses `CD45`. Writing the raw string produces a table of zeros.

It also reads every file rather than the first, so it catches a cohort splitting
into panels. Where the files disagree it names the fix rather than the symptom,
separating autofluorescence and unmixing artefacts from a real marker stained in
only some files, and printing the exact `--ignore-channels` pattern that merges
the cohort. Pass that pattern through: channels are dropped before the panel
fingerprint is computed, which is the only place it works. Anything it cannot
attribute to one of those two mechanical causes it reports as a genuine panel
difference and leaves alone.

`--no-session` is worth adding on any large cohort: `session_state.RData` can
reach several hundred MB and dominate the runtime, and no result depends on it.

Every path inside a flag is a container path. `--outdir` must fall inside a
mounted volume or the output is lost when the container exits. On Windows
PowerShell use `${PWD}`; on Git Bash prefix with `MSYS_NO_PATHCONV=1`.

Annotated templates for both files:

```bash
docker run --rm --entrypoint sh cyraven:1.0.0 -c \
  'cat /usr/local/lib/R/site-library/cyRAVEN/examples/analysis_template.yaml'
```

```bash
docker run --rm --entrypoint sh cyraven:1.0.0 -c \
  'cat /usr/local/lib/R/site-library/cyRAVEN/examples/samples_template.csv'
```

### The one hazard of the single sheet

Subject attributes repeat on every row of that subject, so two rows of one
patient can disagree about that patient's sex. cyRAVEN **stops** and names every
conflict as `subject / column: value vs value` rather than picking one. When
assembling a sheet from several sources, expect this and fix the source data; do
not work around it by deleting rows.

A blank is not a disagreement; it is filled from the rows that have a value.

### Demonstration data

When the user has no data to hand, or wants to see the output shape first:

```bash
mkdir -p demo results
```

```bash
docker run --rm -v "$PWD/demo:/demo" \
  --entrypoint Rscript cyraven:1.0.0 \
  /opt/cyraven/src/inst/scripts/demo_data.R /demo
```

No download is needed: this writes 35 samples from the graft-versus-host disease
cohort that ships with flowCore, five allogeneic transplant recipients across
seven visits stained with a four-colour myeloid panel, together with
`sample_map.csv` and `panel.yaml`. Run it with

```
--group-column cohort --reference-group "GvHD grade 1" --batch-column visit --cluster
```

Three properties constrain what may be said about the result, and you should
state them whenever you present numbers from it. Both groups are transplant
recipients, so the contrast is GvHD grade 1 against grade 3 and not disease
against health. The groups are unbalanced, 7 samples against 28, and grade is a
property of the patient, so a between-group difference is partly a between-donor
difference. And the acquisition carries pulse height with no area, so cyRAVEN
falls back to height and the singlet gate is skipped.

The group contrast on this cohort is null after correction. It demonstrates the
diagnostics on genuinely imperfect clinical data, not sensitivity to a known
effect.

### Iterating on the code

Mount a checkout and set `CYRAVEN_SOURCE`; the entrypoint loads that tree through
`pkgload::load_all()` instead of the installed copy, so edits take effect without
a rebuild.

```bash
docker run --rm -v "$PWD:/src:ro" -v "$PWD/results:/results" \
  -e CYRAVEN_SOURCE=/src cyraven:1.0.0 --dir /data/fcs --outdir /results
```

Full build detail, resource tuning, Windows path handling and container-specific
failures are in `references/docker.md`.

## Start here

| The task | Go to |
|---|---|
| Run it on a batch of files | The commands above, then `references/docker.md` |
| Build or fix the sample sheet | The Inputs section below |
| Write the population YAML | `references/gating.md` §3 |
| Decide arcsinh against logicle | `references/gating.md` §4 |
| Read the outputs | `references/interpretation.md`, in its order |
| A run failed, or the numbers look wrong | `references/troubleshooting.md` |
| Choose between cyRAVEN and cyCONDOR | `references/interpretation.md` §9 |
| A gate retains exactly 90% of its parent | `--auto-fix-gates`, below |
| Turn shares into cell numbers | Absolute counts, below |
| A repeated-measures / multi-visit study | Repeated measures, below |
| A flag was rejected as invalid | `references/troubleshooting.md` |

## Inputs

Only the FCS files are required. Everything else improves what the run can
conclude.

1. **FCS files** in a directory. `--recursive` for subdirectories. Single-stain
   compensation controls are excluded by filename by default; they are instrument
   setup files with no CD45 population and form their own spurious panel group if
   admitted. `--exclude ''` keeps everything.

   Watch the filename pattern. Files named `*.fcs copy` (a common artefact of
   copying on macOS and Windows) do not match the default `[.]fcs$` and are
   listed as skipped rather than silently dropped. The log states the fix:
   `--pattern '[.]fcs( copy)?$'`.

2. **Sample sheet** (`--samples`). Only `file` is required.

   ```csv
   file,sample_id,patient_id,is_control,cohort,sex,age_years,batch,count.Granulocytes
   donor01.fcs,D01,P01,FALSE,Patients,male,41,2025-03-11,5120
   donor01_v2.fcs,D01v2,P01,FALSE,Patients,male,41,2025-06-02,4380
   unstained.fcs,UNS,TRUE,,2025-03-11,
   ```

   Reserved acquisition columns: `file`, `sample_id`, `well`, `patient_id`,
   `panel`, `timepoint`, `is_control`, `fmo_for`, `control_group`.
   Reserved subject columns: `patient_id`, `cohort`, `sex`, `age_years`,
   `date_of_birth`, `height_cm`, `weight_kg`, `infection_focus`, `wbc_per_ul`.
   Anything else is a study variable usable as `--group-column` or
   `--batch-column`. Anything prefixed `count.` is an absolute count.

   A column that is a *measurement* rather than a grouping -- a severity score, a
   laboratory value, an outcome flag -- goes to `--clinical-columns` instead. A
   study column can only group samples, and a SOFA score cannot group anything;
   naming it there associates it with every population and every marker, with the
   test following the column's type (Spearman for numeric, Wilcoxon with Cliff's
   delta for two levels, Kruskal-Wallis with epsilon-squared for more) and each
   effect written with a bootstrap interval. Read
   `clinical_variables_correlation.png` first: p-values are adjusted within each
   variable on the assumption that the variables are separate questions.

   `is_control = TRUE` marks a control tube. Those samples are excluded from
   testing without being recorded as QC failures, and their 99.5th percentile can
   anchor a threshold for a marker with no resolvable density minimum.
   `fmo_for` names the markers a file is the fluorescence-minus-one control for.

   Headers resolve through an alias map covering English and German spellings, so
   `Patient ID` and `Geschlecht` work unchanged. `--write-samples` emits a
   template covering every file in the directory; use it rather than writing one.

3. **Population spec** (`--config`). See `references/gating.md` §3. Without it the
   built-in specification is used, which is a worked example for one
   myeloid-lymphoid panel and not a sensible default for other panels.

4. **Patient table** (`--patient-table`), optional. Covariates for confounder
   screening, and `wbc_per_ul` for absolute concentrations.

`--write-config` derives everything and writes the resulting YAML without running
the analysis, which shows what the pipeline detected before time is spent. It
needs `--outdir` as well as the output path: deriving the parameters means
reading the files, so the early stages run and write their QC output beside it.
Without `--outdir` it stops with `cannot open the connection`.

`--write-samples` and `--write-sample-map` do not need `--outdir`, because they
return before the pipeline starts.

## Options that matter

Analyses that only add output are on by default. Analyses that change an existing
number are opt-in. That rule is deliberate and worth preserving in any change.

```bash
--check                    # validate inputs from FCS headers, exit; run this first
--list-channels            # what each channel resolves to; names the panel fix
--ignore-channels 'AF*'    # drop channels BEFORE the panel fingerprint
--write-samples s.csv      # template covering every input file, then exit
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
--batch-method cluster     # fit the correction per cell type, not per file
--no-marker-group-umaps    # skip marker_umaps_by_group/: each marker pooled,
                           # and split by EVERY category the sheet carries
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
--no-session               # skip session_state.RData; often most of the runtime
--max-events-per-file N    # bound memory; also raises every detection limit
--include-qc-failed        # keep samples that failed staining QC
--read-threads 2           # parallel FCS reads; lower it when memory is tight
--auto-fix-gates           # skip a hierarchy gate the data did not support
--total-counts y.xlsx      # externally measured cell yields -> absolute counts
--split-by-timepoint       # the figure set again per visit, into by_timepoint/
--other                    # draw AND tabulate the unmatched-cell remainder
--no-group-tests           # suppress every between-group test, declared and explore
--clinical-columns sofa    # measurements to associate, not variables to group by
--paired-column patient_id # pair the test within subject
--min-group-n 3            # groups below this are not tested
--p-adjust-display raw     # show unadjusted p beside the adjusted one
--explore-k 30             # FlowSOM cluster count; explore is TOLD this number
--explore-cells-per-sample N  # bounds explore's peak memory, not just its output
--flowjo-export            # write FCS with gate columns for FlowJo
```

`--include-qc-failed` admits QC-failed samples to **every** stage, not only the
embedding; there is no embedding-only form. Under it `qc_status` reads `pass`
for all samples and only the `verdict` column still records that the sample had
no usable parent gate. Say both things whenever you use it.

Four of these change numbers the previous run reported and default to off:
`--drop-unstable-events` changes every count in an affected file,
`--subsample rare` changes which cells are embedded, `--calibration-beads`
changes the units every intensity is expressed in, and `--auto-fix-gates`
changes every count in a file whose hierarchy gate it skips. Say so when you use
them.

### `--auto-fix-gates`

A hierarchy gate whose threshold came from `quantile_fallback` has its retention
decided by the constant `fallback_q = 0.90` rather than by the data, and the
direction of the comparison decides which artefact you get:

```
live_cells   live <- parent & x <  threshold     keeps exactly 90%
cd45_pos     cd45 <- live   & x >  threshold     keeps exactly 10%
```

One constant, two opposite meanings. The signature is unmistakable and worth
checking for on any cohort: a `pct_of_parent` of **exactly 90.00 in every
sample** is not a measurement. Under the flag such a gate is skipped and its
parent carried through, on the grounds that no density minimum means one
population rather than two, and that picking any quantile fabricates a boundary
the data do not contain. `gate_adjustments.csv` names every gate changed, with
the threshold that was not applied and why.

**Population markers are deliberately untouched.** A quantile fallback on CD19 is
a poor threshold, but "no threshold" does not mean "every cell is a B cell". So
this flag repairs the hierarchy and leaves the reason a specification describes
too few cells exactly where it was -- do not present it as a fix for low
coverage.

`--cluster` is the only output that can find a population the specification does
not describe. It is off by default because it costs runtime, not because it is
peripheral to the argument the package makes.

## Repeated measures

If the sample sheet carries a `timepoint` column, `--split-by-timepoint` writes
the whole figure set once per visit into `by_timepoint/<level>/`, beside the
pooled figures.

The one thing to get right when explaining it: the embedding, the thresholds and
the statistics are **not** recomputed per visit -- only the rows drawn are
restricted. So a position on one visit's UMAP is the same position on another's,
and a threshold is the same cut. Running the pipeline separately per visit would
give you three embeddings computed from different cells, which cannot be
compared with each other at all.

Paired views that need every visit at once -- per-patient trajectories, subset
balance, marker intensity by visit -- are written to the run directory as
`timepoint_*.png` whether or not the flag is set, because splitting those by
visit is exactly what destroys them.

It is opt-in because it multiplies the figure count by the number of visits: the
view you want on a repeated-measures study, noise on a cross-sectional one that
merely records a collection date.

**Batch confounding is checked against every design variable**, not just
`--group-column`. On a repeated-measures study the timepoint is often the primary
comparison, and a cohort whose small batches each fall at one visit is confounded
exactly where it matters. `batch_group_confounding.csv` carries one row per
variable, worst first, and the log warns from Cramer's V >= 0.3.

## Absolute counts, and the constraint they lift

Every abundance this package reports is otherwise a **share**, and shares sum to
100. That is not a rounding detail: a frequency table cannot, even in principle,
distinguish one population expanding from every other population contracting,
because the composition is identical either way. The centred log-ratio fixes the
geometry of the test and not this.

Say this whenever a user reads a frequency table as though it measured cell
number -- in sepsis and other lymphopenic states the gap is the entire finding,
because the characteristic result is a global depletion rather than a
redistribution between subsets, and a share cannot represent it.

`--total-counts` takes one measured total yield per **acquisition** (a PBMC
yield, a haemocytometer count) as `.xlsx`, `.csv` or `.tsv`, and adds
`cells_absolute` beside `pct_of_cd45_pos`. Four properties decide whether the
result means anything:

- **Wide sheets blocked by timepoint are read directly**, so a sheet does not
  have to be reshaped by hand: a block-label row, then a `Patient` header beside
  a count header such as `PBMC count x 10^6`, repeated across the sheet.
- **The multiplier comes from the count header.** Where none is stated the
  values are taken at face value and a NOTE names the header it read, because
  guessing between "count" and "count x 10^6" is wrong by six orders of
  magnitude in a table that still looks internally consistent.
- **The join is on `patient_id` AND `timepoint`**, so one patient's d0 and d7
  totals cannot collapse onto each other, and an unlabelled yield facing several
  acquisitions is reported unmatched rather than broadcast across them.
- **The numbers are dual-platform** -- one measurement from this run multiplied
  by one from an instrument it never saw. Published interlaboratory CVs for that
  route are roughly **20-33%**, against 10-16% for single-platform bead counting.
  Treat any difference smaller than that as unresolved. `count_basis` records the
  route on every row carrying one.

Read `total_counts_qc.png` **first**, on its log axis: everything derived
inherits those totals' errors, and a yield in the wrong unit lands decades off
the median there while staying invisible in the derived table. Shares are never
overwritten, so a run without the flag writes exactly what it always did.

## Explore mode

`--explore` runs unsupervised discovery over **every eligible channel**,
ignoring the specification and the parent gate, and writes to
`<outdir>/explore/`. Nothing outside that directory changes, so it is always
safe to add to a run.

Reach for it when the user says any of: a population might be missing from the
spec, they have a new panel and no spec yet, they want to compare against
FlowSOM/Phenograph/cyCONDOR, or they want to know whether a declared label is
lumping two things together.

```bash
--explore                    # alongside the declared analysis
--explore --maybe-learn      # ...and let the two inform each other
--explore-only               # standalone; NO --config or --samples needed
```

Three things to get right when explaining it:

- **`--maybe-learn` is the only link between them.** With it, the declared run
  lends explore its per-sample thresholds, so clusters get a phenotype
  (`CD19+ HLA-DR+ CD3-`) measured against each sample's own cut rather than
  guessed from a heatmap; and explore lends the declared run `spec_gaps.csv`,
  naming populations that span several clusters and clusters nothing covers.
  Without it the two are computed in complete isolation, which is the right
  default: a check on the specification has to be independent of it.
- **A cluster is not a population.** Explore is a hypothesis generator. The
  intended path is: explore finds something → curate
  `explore_suggested_spec.yaml` into a real specification → run the supervised
  path, which is where it gains thresholds, uncertainty and the six checks.
  Under `--maybe-learn` that curation step is assembled for you as
  `suggested_config_next_run.yaml`: the declared specification plus one draft
  entry per cluster `spec_gaps.csv` found uncovered, in the config's own format,
  panel-qualified so `k9` in two panels cannot merge. It is a **draft** -- the
  names are placeholders and nothing in the run that writes it uses it.
- **Point at `explore/explore_report.html`**, which is separate from
  `report.html` and self-contained in the same way.

### What each cluster corresponds to

`explore_cluster_identity.png` and its CSV put clusters as rows and declared
populations as columns, grouped by the population each cluster best matches.

The match is scored by **F1, not by the largest overlap**, and the difference
matters when you explain it. A cluster of 500 cells holding 300 CD4 T cells has
CD4 as its plurality even when those are 5% of all the CD4 T cells in the run --
the label describes the cluster, but the cluster does not describe the label. F1
is high only when both hold, and precision and recall sit beside it so a
disagreement is legible (after Weber & Robinson, Cytometry A 2016;89:1084-1096).

### Immune subsets, and why there is no `--add-subsets`

`explore_cluster_subsets.csv` names each cluster from its own marker profile --
activated, exhausted and homing T cell subsets, NK CD56 bright/dim, HLA-DR low
and CD38+ monocytes, and whatever else the panel can express. Which subsets are
definable is a property of the **panel**, not of immunology, so only definitions
whose every marker is present are proposed, and the canonical groups the panel
cannot reach are named in the log with the markers they would need.

**There is no flag for this and there must not be.** It applies whenever explore
runs, annotates the cluster rather than the cell, and changes no declared output.
`--add-subsets` existed briefly and was removed. Appending the generated subsets
to the declared specification re-scored every cell. It took the frequency table
from 11 populations to 52, and `marker_state` from 275 panels to 1,275, which
then could not be rendered. It also asked a circular question: a subset named
after the marker that defines it is positive for that marker by construction.

If you meet `--add-subsets` in an old command line or an old `run_manifest.txt`,
**drop it**, and warn the user that the run it produced carries the inflated
population list.

### The two cluster heatmaps disagree on purpose

`explore_cluster_heatmap.png` is the share of each cluster above **that sample's
own threshold**. `explore_cluster_median_heatmap.png` is scaled median
expression, hierarchically ordered on both axes. Keep both. Where thresholds are
weak a cluster can read uniformly negative on the first and be cleanly separated
on the second, and that difference is a statement about the gate rather than
about the cells. The median view needs no threshold, which makes it the one to
compare against a tool that has no gating step (cyCONDOR, Phenograph).

**On cluster counts.** FlowSOM *is told* its count through `--explore-k`;
Phenograph *discovers* it from the graph. Matching them removes a gratuitous
difference between two tools; it does not make the partitions equivalent.

`--explore-only` needs nothing but `--dir`. It is the answer to "I have a folder
of FCS files and no idea what is in them".

## The statistics catalogue

Two tables are now in every run, and they answer the question a reader arriving
from an immunophenotyping paper asks first: *where is my t-test?*

- `statistical_methods.csv` names every commonly reported method, Student's and
  Welch's *t*, one-way and two-way ANOVA with Tukey, Mann-Whitney,
  Kruskal-Wallis with Dunn, chi-squared, Bonferroni, diffcyt's edgeR/voom/GLMM,
  and says whether this run computed it and why.
- `normality_tests.csv` is the evidence: Shapiro-Wilk per population per group,
  Brown-Forsythe across groups.

When quoting it, read the `interpretation` column, not `shapiro_p` alone.
Shapiro-Wilk on 4 to 10 donors has almost no power, so a non-significant result
is **not** evidence of normality, which is itself the argument for the rank
tests. Never present the catalogue as a menu to pick a test from.

## The report

`report.html` is the one file to point the user at. It is **self-contained**:
every figure embedded at full resolution, every table embedded in full, no
external reference of any kind. It can be emailed or archived on its own.

Sections are collapsible with a sidebar index. Figures share one display box
whatever their aspect ratio, zoom on click, and download at full resolution.
Tables are searchable, sortable, paged at 10/50/100/all rows, and export to CSV
exactly as filtered and sorted. Section headings state what they report rather
than posing a question, and each carries a description of how to read the
outputs beneath it.

Two things worth telling a reader facing a long report. The sidebar **resizes**:
drag the divider on its right edge, or focus it and use the arrow keys. Its
entries name every figure and table, so they are the longest strings in the
document, and the default width truncates many of them. The width is remembered
per report. A **back-to-top button** appears at the bottom right once they are a
screen down.

Its size is the sum of its figures, so tens of megabytes on a full run. That is
the price of self-containment and is expected, not a fault. A table over 8 MB is
named with its row count instead of embedded, which in practice means only the
per-cell exports.

**A failed run writes it too**, with the diagnosis first: the error, the stage
reached, the log leading up to it, what the message means for the data rather
than for R, and the next action. Everything produced before the failure is
embedded below. When a user reports a failed run, ask for `report.html`, it
carries the log, so they do not need to have kept the terminal output.

## Analytical constraints

These are positions the tool takes, not preferences. Overriding them produces
output that looks correct and is not.

**Thresholds are not transferred between samples.** A global override in the
config applies to every sample and reintroduces the fixed-coordinate bias the
pipeline exists to remove. When one sample's cut is wrong, the remedy is in
`references/troubleshooting.md`.

**`count` is an event count, not a cell number.** It scales with acquisition
duration and with `--max-events-per-file`, so it is not comparable between
samples. Report `pct_of_cd45_pos`, `cells_per_ul` where a haemogram supplied
`wbc_per_ul`, or `cells_absolute` where `--total-counts` supplied a measured
yield.

**A share cannot express a change in size.** Frequencies are constrained to sum
to 100, so no frequency table distinguishes one population expanding from every
other contracting. Where that distinction is the question, it needs
`--total-counts`; the centred log-ratio does not supply it.

**Marker symbols are matched verbatim.** The run compares the specification
against the resolved `$PnS` symbols with hyphens intact, so `TCR-Vd1` and
`HLA-DR` are what belongs in the config. Never rewrite them into R's syntactic
form (`TCR.Vd1`, `HLA.DR`) -- that is what breaks a working configuration, and
every population using the marker then reports UNAVAILABLE with no error raised.
`--check` flags the syntactic form and names the string to copy instead.

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

**Per-sample thresholds do not make samples comparable.** They place each cut
correctly *within* its own sample; they cannot repair a cohort whose staining
intensity differs several-fold between samples. When the positive fraction for a
lineage marker ranges across most of the unit interval, the frequencies are
relative between samples of the same panel and are not an absolute
immunophenotype. `scale_outlier` in `thresholds_used.csv` names every affected
threshold. The instruments for that are `--batch-column` with `--correct-batch`,
or FMO controls declared through `fmo_for` -- none of them a gating change.

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
- `references/troubleshooting.md`: failure modes with their causes, including a
  rejected flag (source newer or older than the command line) and a population
  that reports UNAVAILABLE without erroring

## Keeping this skill current

`NEWS.md` in the repository root is the authority, and its top section is
development work sitting **above** the released 1.0.0 that the image installs.
When something here disagrees with the source, the source wins. The two checks
worth running before trusting this file on an unfamiliar tree:

```bash
# every flag the current source defines
grep -o 'make_option("--[a-zA-Z0-9-]*"' R/cli-options.R | sed 's/.*--//;s/"//' | sort
```

```bash
# what has landed since the last release
awk '/^# cyRAVEN 1\.0\.0/{exit} /^## /{print}' NEWS.md
```
