# cyRAVEN

**Documentation: <https://bhagesh-h.github.io/cyRAVEN/>**

An R package for supervised immunophenotyping of multi-sample flow
cytometry data, with per-sample gate derivation and donor-level
differential abundance and differential state testing.

Turns a directory of FCS files into population frequencies, marker
expression, a shared UMAP and between-group statistics. Every threshold
is placed within the sample it applies to, every frequency carries the
uncertainty of the cut that produced it, and five diagnostics test
whether the gating strategy held.

## 1. Overview

cyRAVEN applies a declarative gating specification to a batch of FCS
files, derives every marker threshold independently within each sample,
quantifies population abundance and marker expression per sample, and
then tests that specification against unsupervised structure recovered
from the same cells.

The design addresses a specific failure mode. Manual gating is the
dominant source of technical variance in multi-sample immunophenotyping:
operators gating identical files report population sizes differing by
approximately 32%, and analyst subjectivity accounts for up to 78% of
technical variability once more than one person is involved (Cadwell et
al., 2021, *PDA J Pharm Sci Technol* 75:33). Fixed gate coordinates
transferred between samples do not remove this variance; they convert it
into a systematic bias that tracks staining intensity.

cyRAVEN removes the analyst from threshold placement and quantifies the
residual uncertainty rather than concealing it. Each cut is resampled
from the events it was derived from and re-derived over the settings
that placed it, and the resulting spread is propagated to every
population that reads it, so a frequency separated by a clean gap and
one sitting on a shoulder are no longer reported to the same apparent
precision.

## 2. Installation

``` r

install.packages("BiocManager")
BiocManager::install("flowCore")
remotes::install_github("bhagesh-h/cyRAVEN")
```

`FlowSOM` is an optional dependency. When absent, self-organising map
clustering falls back to the implementation within the package; all
other behaviour is unchanged.

## 3. Usage

``` r

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

``` bash
Rscript "$(Rscript -e 'cat(system.file("scripts","cyraven.R",package="cyRAVEN"))')" \
  --dir data/ --outdir results/ --recursive \
  --sample-map data/sample_map.csv \
  --config data/panel.yaml \
  --group-column cohort --reference-group "Healthy controls"
```

### 3.1 Input

| Argument | Content |
|----|----|
| `--dir` | Directory of FCS files, one per sample |
| `--sample-map` | CSV linking each filename to a sample identifier and study group |
| `--config` | YAML declaring each population as a set of marker directions |
| `--patient-table` | Optional clinical covariates, keyed on patient identifier |

A population is declared as a conjunction of marker directions evaluated
within the CD45⁺ parent gate:

``` yaml
populations:
  CD4 T cells:
    CD3: above
    CD4: above
    CD8: below
```

### 3.2 Options

Analyses that only add output are enabled by default. Analyses that
alter an existing quantity are opt-in.

``` bash
  --transform logicle       # auto-logicle instead of arcsinh
  --cluster                 # SOM clustering and gate concordance
  --explain-clusters        # learn gate geometry for undescribed clusters
  --external-labels x.csv   # learn gates for a labelling from another tool
  --export-gates            # write learned gates as Gating-ML 2.0
  --write-baseline b.rds    # record where this run placed its thresholds
  --baseline b.rds          # test this run against that record
  --batch-column run_date   # quantify batch structure
  --correct-batch           # correct it, subject to the confounding guard
  --auto-subcluster-k       # silhouette-selected subcluster count
  --save-umap-model         # persist the embedding for later batches
  --no-uncertainty          # skip the gate uncertainty analysis
```

## 4. Methods

### 4.1 Preprocessing

Channels are resolved to marker symbols from `$PnS` with fallback to
`$PnN`. The acquisition spillover matrix is applied when present and its
absence is reported rather than assumed. Events are filtered through a
four-level hierarchy: scatter gate on the log₁₀ FSC-A density minimum,
singlet band on the FSC-H against FSC-A ratio at median ± *k*·MAD,
viability gate at the dye density minimum, and CD45 positivity.

### 4.2 Transformation

Arcsinh with a cofactor estimated from the data by bisection against a
target background interquartile range, or the automatic logicle rule
with the linearisation width taken from the fifth percentile of the
negative population. Logicle parameters are pooled across the panel
rather than fitted per file, which preserves the comparability of
cross-sample medians.

### 4.3 Thresholding

Each marker threshold is placed at the density minimum separating
negative from positive modes, computed within each sample from cells
passing the parent gate. Where no minimum exists the threshold falls
back to a quantile and is labelled as such in the `source` column of
`thresholds_used.csv`. Thresholds deviating from the panel-wide median
are flagged in `threshold_scale_qc.csv` by robust *z*.

### 4.4 Inference

Populations are scored as Boolean conjunctions and aggregated to one
value per sample before testing, following the strategy of Weber et
al. (2019, *Commun Biol* 2:183). Group comparisons use Wilcoxon rank-sum
and Kruskal-Wallis with Benjamini-Hochberg correction within test
families. Effect sizes accompany every p-value. Quantities that cannot
be aggregated to the sample level are reported descriptively without
inferential statistics.

Abundance is additionally tested on centred log-ratios to separate
compositional artefacts from independent change, with concordance
between the two parameterisations reported per population.

### 4.5 Uncertainty

Each threshold is resampled from the parent-gate events it was derived
from, and re-derived across the settings that placed it: histogram
resolution, smoothing width, and how completely two modes must separate
before the gap between them is accepted. The two spreads combine in
quadrature and propagate to every population that reads the marker,
together with the parent gate, which enters all of them because it sets
the denominator.

The reported value is unchanged. Perturbation runs on copies, so a
frequency is identical with the analysis on or off, and what is added is
a second number beside it. `difference_over_gate_u` in the group
comparison states how many multiples of that uncertainty a between-group
difference amounts to; below one, the groups differ by less than the
distance the cut itself moves.

### 4.6 Falsification

Five outputs exist to contradict the specification supplied by the user:

**Threshold drift** tests whether per-sample thresholds differ
systematically between study groups. A flagged marker indicates that
part of any abundance difference is definitional rather than biological.

**Gate uncertainty** states how far each cut moves under resampling of
the events that determined it, and what that costs every population
reading it. Where a between-group difference is smaller than the
movement of its own threshold, the comparison is not resolvable by this
strategy at this staining quality.

**Phenotype concordance** displays measured marker intensity against the
populations declared to express them, reading the same specification
used for scoring.

**Gate against cluster concordance** cross-tabulates the declared labels
against an unsupervised SOM clustering computed without reference to
them. A cluster dominated by the unassigned label identifies a
population absent from the specification; a label containing far fewer
cells than the cluster it dominates identifies a misplaced threshold.

**Learned gate geometry** (`--explain-clusters`) derives a sequence of
two-marker convex polygon gates selecting an undescribed cluster, with
performance measured on held-out cells and reported alongside the
resubstitution value.

### 4.7 Conformance

A baseline written from an accepted run records where each marker’s
threshold sat, how variable it was, how often it required the quantile
fallback, and what the populations came out at. Later runs are measured
against it.

This is a second reference, not a replacement for the first. The
within-run check compares each sample against its peers and identifies
one deviant tube. It cannot identify a cohort that moved as a whole,
since the leave-one-out peer median moves with it, which is the
situation after a laser service or a reagent lot change. A changed
transform is reported once as not comparable rather than as every marker
having drifted, and a redefined population is reported as a different
measurement rather than a moving one.

### 4.8 Gate learning from external labels

`--external-labels` accepts a cell labelling produced elsewhere, for
example a cyCONDOR clustering, and derives an executable gating strategy
for each label. The join is on sample and event index rather than row
position, since the two tools subsample independently.

Transferability is measured by refitting the strategy with one donor
withheld and scoring it on that donor. Cells held out of a fit come from
the same donors acquired in the same tubes, so they share every source
of between-donor variation the gate will meet in use; the per-donor
minimum is the quantity that predicts the next patient. `--export-gates`
writes the result as ISAC Gating-ML 2.0 in the linear units the FCS file
stores, with polygon edges subdivided before inversion because an edge
that is straight on the analysis scale is a curve in linear units.

### 4.9 Batch correction

Batch structure is quantified as iLISI against a permutation null and
batch against group association as Cramér’s *V*. Correction by monotone
quantile alignment proceeds only when *V* falls below a configurable
threshold. Above it the operation is refused, on the grounds that
removing the batch and removing the effect under study are not separable
at that level of confounding. Correction is applied to the shared
embedding and clustering, not to per-sample thresholds, which are
already batch-local by construction.

## 5. Output

A default run produces 25 tables and 22 figures.

| File | Content |
|----|----|
| `population_frequencies.csv` | Population percentage of parent, per sample, with its gate uncertainty |
| `population_marker_mfi.csv` | Median intensity and percent positive, per sample, population and marker |
| `group_comparison_stats.csv` | Per-population test with Cliff’s delta, adjusted p, and the difference expressed in units of gate uncertainty |
| `thresholds_used.csv` | Every threshold, its derivation and its outlier status |
| `threshold_uncertainty.csv` | Sampling and method components of each threshold, and how often resampling found it |
| `uncertainty_budget.csv` | Which threshold each population’s uncertainty comes from |
| `specification_conformance.csv` | This run against an accepted baseline, written by `--baseline` |
| `gate_transferability.csv` | Held-out-donor performance of a learned gate, written by `--external-labels` |
| `cluster_gate_agreement_*.csv` | Declared labels against unsupervised clusters |
| `batch_group_confounding.csv` | Cramér’s *V* and correction verdict |
| `umap_overview.png` | Shared embedding by population and by sample |
| `frequency_uncertainty.png` | Between-sample spread against within-sample gate uncertainty |
| `population_marker_heatmap.png` | Marker intensity against declared populations |
| `run_manifest.txt` | Package versions, git commit, invocation, options |

`recon_diagnostics.png` and `gating_qc.png` should be inspected before
any quantity derived from them.

## 6. Documentation

| Article | Content |
|----|----|
| [Workflow](https://bhagesh-h.github.io/cyRAVEN/articles/cyRAVEN.html) | The ten pipeline stages with the function implementing each, and executable examples of cofactor estimation, density minimum detection and group comparison |
| [Gating](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html) | Gate hierarchy and its behaviour when CD45 or viability markers are absent; per-sample thresholding and interpretation of the `source` column; how far a threshold moves under resampling and what that costs the population; stability against a baseline from an earlier run; population specification syntax including three-level markers; arcsinh against logicle with the CCR7 memory-subset case |
| [Diagnostics](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.html) | Twelve checks in reading order: gate inspection, staining QC, phenotype concordance, threshold drift, gate uncertainty and its budget, the four gate-cluster concordance patterns, held-out-donor gate transferability, learned gate geometry, covariate screening, iLISI batch quantification with the Cramér’s *V* refusal rule, conformance against a baseline, and run provenance |
| [Output](https://bhagesh-h.github.io/cyRAVEN/articles/outputs.html) | Every file the pipeline writes, its columns, the flag producing it, and why event counts are not cell counts |
| [Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.html) | Sample-level aggregation with a worked example; rank tests against moderated *t* and the sample size at which the trade reverses; compositional constraint and the limits of the log-ratio; covariate diagnosis against adjustment; multiplicity within test families; the difference expressed in units of gate uncertainty |
| [Interoperability](https://bhagesh-h.github.io/cyRAVEN/articles/with-cycondor.html) | Method selection by question type; three worked analyses; handing a cyCONDOR clustering to cyRAVEN to get an executable gate with held-out-donor performance; running both from one sample map; the four sources of divergence between their embeddings; cluster count as an analytical choice |
| [Scope](https://bhagesh-h.github.io/cyRAVEN/articles/scope.html) | Nine excluded methods with the reasoning and the condition under which each becomes appropriate |

## 7. Citation

``` r

citation("cyRAVEN")
```

Cite FlowSOM (Van Gassen et al., 2015, *Cytometry A* 87:636) when using
the unsupervised clustering, and diffcyt (Weber et al., 2019, *Commun
Biol* 2:183) for the aggregation strategy underlying the differential
state tests.

## 8. Licence

GPL-3. See [LICENSE](https://bhagesh-h.github.io/cyRAVEN/LICENSE).
