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

Turns a directory of FCS files into population frequencies, marker expression, a
shared UMAP and between-group statistics. Every threshold is placed within the
sample it applies to, and four diagnostics test whether the gating strategy held.

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

cyRAVEN removes the analyst from threshold placement and reports the residual
uncertainty rather than concealing it.

## 2. Installation

```r
install.packages("BiocManager")
BiocManager::install("flowCore")
remotes::install_github("bhagesh-h/cyRAVEN")
```

`FlowSOM` is an optional dependency. When absent, self-organising map clustering
falls back to the implementation within the package; all other behaviour is
unchanged.

## 3. Usage

```r
library(cyRAVEN)

run_cyraven(list(
  dir             = "data/",
  outdir          = "results/",
  recursive       = TRUE,
  sample_map      = "data/sample_map.csv",
  config          = "data/panel.yaml",
  group_column    = "cohort",
  reference_group = "Healthy controls"
))
```

Command-line equivalent:

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts","cyraven.R",package="cyRAVEN"))')" \
  --dir data/ --outdir results/ --recursive \
  --sample-map data/sample_map.csv \
  --config data/panel.yaml \
  --group-column cohort --reference-group "Healthy controls"
```

### 3.1 Input

| Argument | Content |
|---|---|
| `--dir` | Directory of FCS files, one per sample |
| `--sample-map` | CSV linking each filename to a sample identifier and study group |
| `--config` | YAML declaring each population as a set of marker directions |
| `--patient-table` | Optional clinical covariates, keyed on patient identifier |

A population is declared as a conjunction of marker directions evaluated within
the CD45<sup>+</sup> parent gate:

```yaml
populations:
  CD4 T cells:
    CD3: above
    CD4: above
    CD8: below
```

### 3.2 Options

Analyses that only add output are enabled by default. Analyses that alter an
existing quantity are opt-in.

```bash
  --transform logicle       # auto-logicle instead of arcsinh
  --cluster                 # SOM clustering and gate concordance
  --explain-clusters        # learn gate geometry for undescribed clusters
  --batch-column run_date   # quantify batch structure
  --correct-batch           # correct it, subject to the confounding guard
  --auto-subcluster-k       # silhouette-selected subcluster count
  --save-umap-model         # persist the embedding for later batches
```

## 4. Methods

### 4.1 Preprocessing

Channels are resolved to marker symbols from `$PnS` with fallback to `$PnN`. The
acquisition spillover matrix is applied when present and its absence is reported
rather than assumed. Events are filtered through a four-level hierarchy: scatter
gate on the log<sub>10</sub> FSC-A density minimum, singlet band on the FSC-H
against FSC-A ratio at median ± *k*·MAD, viability gate at the dye density
minimum, and CD45 positivity.

### 4.2 Transformation

Arcsinh with a cofactor estimated from the data by bisection against a target
background interquartile range, or the automatic logicle rule with the
linearisation width taken from the fifth percentile of the negative population.
Logicle parameters are pooled across the panel rather than fitted per file, which
preserves the comparability of cross-sample medians.

### 4.3 Thresholding

Each marker threshold is placed at the density minimum separating negative from
positive modes, computed within each sample from cells passing the parent gate.
Where no minimum exists the threshold falls back to a quantile and is labelled as
such in the `source` column of `thresholds_used.csv`. Thresholds deviating from
the panel-wide median are flagged in `threshold_scale_qc.csv` by robust *z*.

### 4.4 Inference

Populations are scored as Boolean conjunctions and aggregated to one value per
sample before testing, following the strategy of Weber et al. (2019,
*Commun Biol* 2:183). Group comparisons use Wilcoxon rank-sum and
Kruskal-Wallis with Benjamini-Hochberg correction within test families. Effect
sizes accompany every p-value. Quantities that cannot be aggregated to the sample
level are reported descriptively without inferential statistics.

Abundance is additionally tested on centred log-ratios to separate compositional
artefacts from independent change, with concordance between the two
parameterisations reported per population.

### 4.5 Falsification

Four outputs exist to contradict the specification supplied by the user:

**Threshold drift** tests whether per-sample thresholds differ systematically
between study groups. A flagged marker indicates that part of any abundance
difference is definitional rather than biological.

**Phenotype concordance** displays measured marker intensity against the
populations declared to express them, reading the same specification used for
scoring.

**Gate against cluster concordance** cross-tabulates the declared labels against
an unsupervised SOM clustering computed without reference to them. A cluster
dominated by the unassigned label identifies a population absent from the
specification; a label containing far fewer cells than the cluster it dominates
identifies a misplaced threshold.

**Learned gate geometry** (`--explain-clusters`) derives a sequence of
two-marker convex polygon gates selecting an undescribed cluster, with
performance measured on held-out cells and reported alongside the resubstitution
value.

### 4.6 Batch correction

Batch structure is quantified as iLISI against a permutation null and batch
against group association as Cramér's *V*. Correction by monotone quantile
alignment proceeds only when *V* falls below a configurable threshold. Above it
the operation is refused, on the grounds that removing the batch and removing the
effect under study are not separable at that level of confounding. Correction is
applied to the shared embedding and clustering, not to per-sample thresholds,
which are already batch-local by construction.

## 5. Output

A default run produces 23 tables and 20 figures.

| File | Content |
|---|---|
| `population_frequencies.csv` | Population percentage of parent, per sample |
| `population_marker_mfi.csv` | Median intensity and percent positive, per sample, population and marker |
| `group_comparison_stats.csv` | Per-population test with Cliff's delta and adjusted p |
| `thresholds_used.csv` | Every threshold, its derivation and its outlier status |
| `cluster_gate_agreement_*.csv` | Declared labels against unsupervised clusters |
| `batch_group_confounding.csv` | Cramér's *V* and correction verdict |
| `umap_overview.png` | Shared embedding by population and by sample |
| `population_marker_heatmap.png` | Marker intensity against declared populations |
| `run_manifest.txt` | Package versions, git commit, invocation, options |

`recon_diagnostics.png` and `gating_qc.png` should be inspected before any
quantity derived from them.

## 6. Documentation

| Article | Content |
|---|---|
| [Workflow](https://bhagesh-h.github.io/cyRAVEN/articles/cyRAVEN.html) | The nine pipeline stages with the function implementing each, and executable examples of cofactor estimation, density minimum detection and group comparison |
| [Gating](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html) | Gate hierarchy and its behaviour when CD45 or viability markers are absent; per-sample thresholding and interpretation of the `source` column; population specification syntax including three-level markers; arcsinh against logicle with the CCR7 memory-subset case |
| [Diagnostics](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.html) | Nine checks in reading order: gate inspection, staining QC, phenotype concordance, threshold drift, the four gate-cluster concordance patterns, learned gate geometry, covariate screening, iLISI batch quantification with the Cramér's *V* refusal rule, and run provenance |
| [Output](https://bhagesh-h.github.io/cyRAVEN/articles/outputs.html) | Every file the pipeline writes, its columns, the flag producing it, and why event counts are not cell counts |
| [Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.html) | Sample-level aggregation with a worked example; rank tests against moderated *t* and the sample size at which the trade reverses; compositional constraint and the limits of the log-ratio; covariate diagnosis against adjustment; multiplicity within test families |
| [Interoperability](https://bhagesh-h.github.io/cyRAVEN/articles/with-cycondor.html) | Method selection by question type; three worked analyses; running cyRAVEN and cyCONDOR from one sample map; the four sources of divergence between their embeddings; cluster count as an analytical choice |
| [Scope](https://bhagesh-h.github.io/cyRAVEN/articles/scope.html) | Nine excluded methods with the reasoning and the condition under which each becomes appropriate |

## 7. Citation

```r
citation("cyRAVEN")
```

Cite FlowSOM (Van Gassen et al., 2015, *Cytometry A* 87:636) when using the
unsupervised clustering, and diffcyt (Weber et al., 2019, *Commun Biol* 2:183)
for the aggregation strategy underlying the differential state tests.

## 8. Licence

GPL-3. See [LICENSE](LICENSE).
