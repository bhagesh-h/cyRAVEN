# cyRAVEN <img src="man/figures/logo.png" alt="cyRAVEN logo" width="200" align="right"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/pkgdown.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Documentation: <https://bhagesh-h.github.io/cyRAVEN/>**

An R package for supervised immunophenotyping of multi-sample flow cytometry
data, with per-sample gate derivation and donor-level differential abundance and
differential state testing.

A directory of FCS files becomes population frequencies, marker expression, a
shared UMAP and between-group statistics. Every threshold is placed within the
sample it applies to, every frequency carries the uncertainty of the cut and of
the event count behind it, and six diagnostics test whether the gating
specification held.

## 1. Overview

cyRAVEN applies a declarative gating specification to a batch of FCS files,
derives every marker threshold independently within each sample, quantifies
population abundance and marker expression per sample, and then tests that
specification against unsupervised structure recovered from the same cells.

The design addresses a specific failure mode. Manual gating is the dominant
source of technical variance in multi-sample immunophenotyping: operators gating
identical files report population sizes differing by approximately 32%, and
analyst subjectivity accounts for up to 78% of technical variability once more
than one person is involved (Cadwell et al., 2021, *PDA J Pharm Sci Technol*
75:33). Fixed gate coordinates transferred between samples do not remove this
variance; they convert it into a systematic bias that tracks staining intensity.

cyRAVEN removes the analyst from threshold placement and quantifies the residual
uncertainty rather than concealing it. Each cut is resampled from the events it
was derived from and re-derived over the settings that placed it, and the
resulting spread is propagated to every population that reads it, so a frequency
separated by a clean gap and one sitting on a shoulder are no longer reported to
the same apparent precision.

## 2. Quick start

Docker is the supported execution path. Two of cyRAVEN's dependencies determine
numerical results rather than convenience: `uwot` is a stochastic embedder whose
output varies with its own version and with the BLAS beneath it, and every gate
is placed at a kernel density minimum, so a different density implementation
moves thresholds and therefore the frequencies that would be published. The image
pins R 4.4.3, a dated CRAN snapshot and Bioconductor 3.20 for that reason. Local
installation is documented in [section 7](#7-installation-without-docker).

The commands below reproduce a complete run on a public dataset in four steps and
require only Docker.

### 2.1 Build the image

```bash
git clone https://github.com/bhagesh-h/cyRAVEN.git
cd cyRAVEN
docker build -f inst/scripts/Dockerfile -t cyraven:0.3.0 .
```

The first build compiles the dependency stack and takes 15 to 25 minutes. It
finishes by running `--help` and printing every package version, so a broken
image fails at build time rather than during an analysis.

### 2.2 Fetch the demonstration cohort

```bash
mkdir -p demo results
docker run --rm -v "$PWD/demo:/demo" \
  --entrypoint Rscript cyraven:0.3.0 \
  /opt/cyraven/src/inst/scripts/demo_data.R /demo
```

This downloads two CytoTrol control acquisitions from RGLab's
[flowWorkspaceData](https://github.com/RGLab/flowWorkspaceData) repository
(Artistic-2.0, no account required), partitions each at random into four
pseudo-samples, and writes `sample_map.csv` and `panel.yaml` beside them.

The `cohort` column it writes is randomised. Every pseudo-sample originates from
the same tube, so the true between-group difference is zero by construction. The
example is therefore a calibration check: the correct outcome is that no
population differs significantly, and any result reported as significant is a
false positive. It cannot demonstrate sensitivity, because no real effect is
present.

### 2.3 Run the pipeline

```bash
docker run --rm \
  -v "$PWD/demo:/data:ro" \
  -v "$PWD/results:/results" \
  cyraven:0.3.0 \
  --dir /data/fcs \
  --sample-map /data/sample_map.csv \
  --config /data/panel.yaml \
  --group-column cohort --reference-group GroupA \
  --outdir /results
```

Runtime is about three minutes for these eight samples on four cores. Every path
inside a flag is a path inside the container, and `--outdir` must fall within a
mounted volume or the output is discarded when the container exits.

On Windows PowerShell, substitute `${PWD}` for `$PWD`.

### 2.4 Read the result

```bash
ls results/
```

Inspect `recon_diagnostics.png` and `gating_qc.png` before any quantity derived
from them. A threshold placed on a distribution shoulder rather than a density
minimum is visible there and in no downstream table.

The reading order for the remaining outputs is given in the
[Diagnostics article](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.html);
it is not arbitrary, because each stage can invalidate the ones after it.
[Section 5](#5-output) shows what this run produces.

## 3. Input

| Argument | Content |
|---|---|
| `--dir` | Directory of FCS files, one per sample |
| `--sample-map` | CSV linking each filename to a sample identifier and study group |
| `--config` | YAML declaring each population as a set of marker directions |
| `--patient-table` | Optional clinical covariates, keyed on patient identifier |

Only the FCS files are required. A population is declared as a conjunction of
marker directions evaluated within the CD45<sup>+</sup> parent gate:

```yaml
populations:
  CD4 T cells:
    CD3: above
    CD4: above
    CD8: below
```

`--write-sample-map` and `--write-config` emit templates and exit. The
specification syntax, including three-level markers and the `any_of` form, is
documented in the
[Gating article](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html).

## 4. Options

Analyses that only add output are enabled by default. Analyses that alter an
existing quantity are opt-in.

```bash
  --transform logicle       # auto-logicle instead of arcsinh
  --cluster                 # SOM clustering and gate concordance
  --explain-clusters        # learn gate geometry for undescribed clusters
  --external-labels x.csv   # learn gates for a labelling from another tool
  --export-gates            # write learned gates as Gating-ML 2.0
  --write-baseline b.rds    # record where this run placed its thresholds
  --baseline b.rds          # test this run against that record
  --batch-column run_date   # quantify batch structure
  --correct-batch           # correct it, subject to the confounding guard
  --lod-events 20           # events below which a population is not detected
  --auto-subcluster-k       # silhouette-selected subcluster count
  --save-umap-model         # persist the embedding for later batches
  --no-uncertainty          # skip the gate uncertainty analysis
```

`docker run --rm cyraven:0.3.0 --help` lists every option.

## 5. Output

The run in [section 2](#2-quick-start) writes 16 tables and 18 figures, together
with `run_manifest.txt`, `miflowcyt.md` and `session_state.RData`. The count
varies with the panel and the flags in force: outputs whose inputs are absent are
skipped and the omission is logged.

The figures below are the unmodified output of that run. Two properties of the
demonstration panel affect how they read. It carries no CD45 and no viability
dye, so the parent gate is scatter and singlets rather than live leukocytes,
while the axis label and the `pct_of_cd45_pos` column retain the CD45 naming. And
the group labels are randomised, so the correct result is that nothing differs.

### 5.1 Gate inspection

Per sample: the scatter gate boundary, the singlet band, and every marker
threshold drawn on the density it was derived from. Read before any quantity
derived from it, because a threshold placed on a shoulder rather than a minimum
appears here and in no table.

<img src="man/figures/demo_gating_qc.png" alt="Per-sample marker densities with the derived threshold marked on each" width="100%"/>

### 5.2 Abundance and its uncertainty

Each per-sample frequency carries the standard uncertainty propagated from the
thresholds behind it.

In this run 7 of 10 populations have an uncertainty as wide as the spread between
samples, which is the expected result here: the pseudo-samples are drawn from one
tube, so the between-sample spread is close to zero and gate placement is the
larger term. On a real cohort the same reading identifies populations whose
apparent variation is the cut moving rather than the biology.

<img src="man/figures/demo_frequency_uncertainty.png" alt="Population frequencies per sample with gate placement uncertainty as horizontal bars" width="100%"/>

### 5.3 Measurability at this acquisition depth

Placement precision and counting sufficiency are separate guarantees. A cut
through a wide empty gap is well determined however few events lie beyond it, so
each population is also classified against the limits of detection and
quantification set by its own parent gate.

All ten populations clear both limits in all eight samples here, at roughly
30,000 events per sample. The figure is diagnostic when they do not: a population
mostly below the limit of quantification cannot be recovered by re-gating,
because the numerator is small for want of acquired cells.

<img src="man/figures/demo_detection_limits.png" alt="Populations counted by whether their event count clears the limits of detection and quantification" width="100%"/>

### 5.4 Phenotype concordance

Marker intensity against the declared populations, z-scored across the run.
Identity is declared rather than inferred, so this figure serves the inverse
function of its equivalent in clustering-first analysis: a column labelled CD4 T
cells that does not show elevated CD4 and depressed CD8 falsifies the gate that
produced it.

<img src="man/figures/demo_population_marker_heatmap.png" alt="Heatmap of marker intensity against declared populations" width="100%"/>

### 5.5 Shared embedding

One UMAP per marker panel, computed across all samples so that cross-sample
comparison is meaningful.

<img src="man/figures/demo_umap_overview.png" alt="Shared UMAP embedding coloured by population, sample and group" width="100%"/>

### 5.6 Calibration of this run

Because the group labels are randomised, the number of significant results is a
measurement of the false-positive rate rather than a finding. Of 10 populations
tested, 1 reached raw *p* < 0.05 and none survived Benjamini-Hochberg correction,
which is what a correctly calibrated test returns at this number of comparisons.
`compositional_concordance.csv` classifies that one as significant on raw
percentages only, the configuration produced by the compositional constraint
rather than by an independent change.

### 5.7 Principal tables

| File | Content |
|---|---|
| `population_frequencies.csv` | Population percentage of parent, per sample, with gate uncertainty, counting uncertainty, and detection verdict |
| `population_marker_mfi.csv` | Median intensity and percent positive, per sample, population and marker |
| `group_comparison_stats.csv` | Per-population test with Cliff's delta, adjusted p, and the difference expressed in units of gate and of total uncertainty |
| `thresholds_used.csv` | Every threshold, its derivation and its outlier status |
| `threshold_uncertainty.csv` | Sampling and method components of each threshold, and how often resampling recovered it |
| `uncertainty_budget.csv` | Which threshold each population's uncertainty comes from |
| `specification_conformance.csv` | This run against an accepted baseline, written by `--baseline` |
| `marker_batch_drift.csv` | Per-marker distributional distance between batches, written by `--batch-column` |
| `gate_transferability.csv` | Held-out-donor performance of a learned gate, written by `--external-labels` |
| `cluster_gate_agreement_*.csv` | Declared labels against unsupervised clusters |
| `run_manifest.txt` | Package versions, git commit, invocation, options |
| `miflowcyt.md` | ISAC-structured report of the instrument configuration and the analysis |

Every file the pipeline writes, with its columns and the flag that produces it,
is documented in the
[Output article](https://bhagesh-h.github.io/cyRAVEN/articles/outputs.html).

## 6. Method summary

Each item below is treated in depth in the linked article.

**Preprocessing and gating.** Channels are resolved to marker symbols from `$PnS`
with fallback to `$PnN`, and the acquisition spillover matrix is applied when
present. Events pass a four-level hierarchy: scatter, singlet band, viability,
and CD45 positivity. Every marker threshold is then placed at the density minimum
separating negative from positive modes, computed within each sample.
[Gating](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html)

**Uncertainty.** Each threshold is resampled from the events it was derived from
and re-derived across the settings that placed it. The two spreads combine in
quadrature and propagate to every population reading the marker, together with
the parent gate. A counting term derived from the Wilson interval is reported
beside them. The reported frequency is never altered: perturbation runs on
copies.
[Gating §2.4](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html)

**Inference.** Populations are aggregated to one value per sample before testing,
following Weber et al. (2019, *Commun Biol* 2:183), so donors rather than cells
are the replicates. Group comparisons use rank tests with Benjamini-Hochberg
correction within test families, and abundance is additionally tested on centred
log-ratios to separate compositional artefacts from independent change.
[Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.html)

**Falsification.** Six outputs exist to contradict the supplied specification:
threshold drift between groups, gate and counting uncertainty, phenotype
concordance, gate against unsupervised cluster concordance, learned gate
geometry, and conformance against a baseline from an accepted earlier run.
[Diagnostics](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.html)

## 7. Installation without Docker

Reproducibility of numerical output is not guaranteed outside the pinned
environment, for the reasons given in [section 2](#2-quick-start).

```r
install.packages("BiocManager")
BiocManager::install("flowCore")
remotes::install_github("bhagesh-h/cyRAVEN")
```

```r
library(cyRAVEN)

run_cyraven(list(
  dir             = "demo/fcs",
  outdir          = "results",
  sample_map      = "demo/sample_map.csv",
  config          = "demo/panel.yaml",
  group_column    = "cohort",
  reference_group = "GroupA"
))
```

`FlowSOM` is an optional dependency. When absent, self-organising map clustering
falls back to the implementation within the package and all other behaviour is
unchanged.

## 8. Documentation

| Article | Content |
|---|---|
| [Workflow](https://bhagesh-h.github.io/cyRAVEN/articles/cyRAVEN.html) | The ten pipeline stages with the function implementing each, and executable examples of cofactor estimation, density minimum detection and group comparison |
| [Gating](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html) | Gate hierarchy and its behaviour when CD45 or viability markers are absent; per-sample thresholding and the `source` column; threshold precision and counting sufficiency; stability against a baseline; specification syntax; arcsinh against logicle |
| [Diagnostics](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.html) | The checks in reading order: gate inspection, staining QC, phenotype concordance, threshold drift, gate uncertainty, detection limits, the four gate-cluster concordance patterns, held-out-donor transferability, learned gate geometry, covariate screening, batch structure, conformance, and provenance |
| [Output](https://bhagesh-h.github.io/cyRAVEN/articles/outputs.html) | Every file the pipeline writes, its columns, the flag producing it, and why event counts are not cell counts |
| [Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.html) | Sample-level aggregation; rank tests against moderated *t* and the sample size at which the trade reverses; the compositional constraint; covariate diagnosis against adjustment; multiplicity within test families; differences expressed in units of uncertainty |
| [Interoperability](https://bhagesh-h.github.io/cyRAVEN/articles/with-cycondor.html) | Method selection by question type; handing a cyCONDOR clustering to cyRAVEN to obtain an executable gate with held-out-donor performance; the four sources of divergence between their embeddings |
| [Claude skill](https://bhagesh-h.github.io/cyRAVEN/articles/claude-skill.html) | Installing and using the bundled Claude Code skill, which executes through Docker by default |
| [Scope](https://bhagesh-h.github.io/cyRAVEN/articles/scope.html) | Nine excluded methods with the reasoning and the condition under which each becomes appropriate |

## 9. Citation

```r
citation("cyRAVEN")
```

Cite FlowSOM (Van Gassen et al., 2015, *Cytometry A* 87:636) when using the
unsupervised clustering, and diffcyt (Weber et al., 2019, *Commun Biol* 2:183)
for the aggregation strategy underlying the differential state tests.

## 10. Licence

GPL-3. See [LICENSE](LICENSE).
