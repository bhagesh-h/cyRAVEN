# cyRAVEN

<img src="man/figures/logo.png" alt="cyRAVEN" width="200" align="right"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**cyRAVEN** takes a directory of FCS files and returns a gated, embedded,
statistically tested result set: per-sample gate thresholds derived from each
sample's own marker density, populations scored from a declarative YAML
specification, one shared UMAP per marker panel, and differential abundance and
differential state tested on **per-sample** summaries so that donors, not cells,
are the replicates.

Every stage is an exported function, so the pipeline can be run whole, run in
part, inspected between steps, or have any one rule replaced without forking it.

## 1. What this is, in one paragraph

Most cytometry toolboxes cluster first and name the clusters afterwards. cyRAVEN
gates first, against a threshold specification you can read, and then spends its
effort on the question that ordering raises: **can the gates be trusted?** So the
package ships the checks that a gating-first workflow needs and a clustering-first
one has no place for — a phenotype heatmap that can contradict a gate label, a
threshold-drift test that asks whether the *gate* moved with cohort, and an
unsupervised clustering run purely as a cross-check that is allowed to disagree.
Alongside that sit the things any small-cohort study needs and rarely gets:
sample-level differential state, compositional (CLR) frequency testing,
confounding diagnostics, and a batch effect that is measured rather than silently
corrected.

## 2. Installing and running

### Installation

```r
# flowCore comes from Bioconductor, so put it on the search path first
install.packages("BiocManager")
BiocManager::install("flowCore")

# then the package
remotes::install_github("bhagesh-h/cyRAVEN")
```

`FlowSOM` is suggested rather than required: without it the SOM clustering falls
back to the package's own implementation, and everything else is unaffected.

### Running

One call does everything:

```r
library(cyRAVEN)

run_cyraven(list(
  dir             = "data/",
  outdir          = "results/",
  recursive       = TRUE,
  sample_map      = "data/sample_map.csv",
  patient_table   = "data/patient_table.csv",
  config          = system.file("config", "config_cohorts.yaml", package = "cyRAVEN"),
  group_column    = "cohort",
  reference_group = "Healthy controls"
))
```

The same thing from a shell, via the command-line front end the package installs:

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts","cyraven.R",package="cyRAVEN"))')" \
  --dir data/ --outdir results/ --recursive \
  --sample-map data/sample_map.csv \
  --patient-table data/patient_table.csv \
  --group-column cohort --reference-group "Healthy controls"
```

Every analysis that only *adds* output is **on by default**; anything that
changes an existing number is **off by default** or has an escape hatch:

```bash
  --cluster                 # unsupervised clustering + gate cross-check
  --explain-clusters        # learn a gating strategy for undescribed clusters
  --batch-column run_date   # batch-mixing diagnostic
  --rank-ancova             # exploratory covariate-adjusted tests
  --auto-subcluster-k       # data-driven subcluster count
  --save-umap-model results/umap.model
```

### Verbosity

Progress goes through `message()`, so `suppressMessages()` works, and the level
is an option rather than an argument threaded through every call:

```r
options(cyRAVEN.verbose = "none")     # or "inform" (default), or "debug"
```

### Memory

The main operational constraint: `--memory` must be sized against
`--max-events-per-file`, because gating holds a transformed matrix per sample. At
25 files with **no** events cap, 10 GB is killed at STEP 4 (exit 137, no message
— the kernel does not ask). A production run of that size wants
`--max-events-per-file 300000` with `--memory 14g`.

`--cluster` adds to the peak, but far less than gating does: it clusters the
already-subsampled embedding matrix, not the full event set.

### Docker

```bash
docker build -f inst/scripts/Dockerfile -t cyraven:0.1.0 .
docker run --rm -v "$PWD:/data:ro" -v "$PWD/results:/results" cyraven:0.1.0 \
  --dir /data --recursive --outdir /results
```

The image pins R, CRAN and Bioconductor, and installs exactly the dependencies
`DESCRIPTION` declares — there is no second list to drift out of step. Set
`-e CYRAVEN_SOURCE=/src` with your checkout mounted at `/src` to run a live source
tree without rebuilding.

### Package layout

| Path | Role |
|---|---|
| `R/` | the package: 24 thematic source files, every stage an exported function |
| `inst/scripts/cyraven.R` | command-line front end — a thin wrapper over `run_cyraven()` |
| `inst/scripts/Dockerfile` | pinned, reproducible environment |
| `inst/config/` | the worked-example population specification |
| `inst/examples/` | template `--sample-map` / `--patient-table` / `--absolute-counts` files |
| `tests/testthat/` | planted-truth fixtures plus guards on each optimisation |
| `vignettes/` | getting started, and the statistical rationale |

## 3. Implemented

### 3.1 Differential state — testing marker intensity, correctly

**The gap.** Abundance testing tells you a population changed *size*
(`group_comparison.png`, Wilcoxon on per-sample abundance — proper sample-level
statistics). It could *show* you a marker shifted inside a population
(`umap_multigraph_overlay.png`) and *rank* those shifts
(`subcluster_marker_shifts.csv`). It could not **test** one, because the ranking
is computed over pooled cells: n is the number of cells, so one donor
contributing 4,000 cells can produce a large Cliff's delta alone. That is
pseudoreplication. It is common enough that the cell-level Wilcoxon helpers
in this field ship with a warning telling you not to interpret their p-values.

**What was done.** The established fix, from diffcyt's differential-state
pathway (Weber et al. 2019): aggregate each marker to one value per sample per
population *first*, then test those. This pipeline already computed that table —
`population_marker_mfi.csv` is exactly per-sample × population × marker medians —
so only the testing was missing.

Two measures are tested, because they detect different things:

| measure | moves when |
|---|---|
| `median_asinh` | the whole population shifts brighter or dimmer |
| `pct_positive` | a *subset* of the population turns positive |

A bimodal shift (30% of cells go bright, 70% unchanged) moves `pct_positive`
sharply and the median barely at all. Reporting one measure makes the other kind
of change invisible, so the measure is a column, not an assumption.

**Why Wilcoxon and not limma/LMM as diffcyt uses.** limma's moderated *t*
borrows variance across markers via an empirical-Bayes prior — a good trade at
mass-cytometry marker counts (30–40) and designed-experiment sample counts. At a
dozen markers and single-digit donors per cohort, the prior is estimated from too
few markers to help and the normality assumption is doing real work rather than
being a formality. The abundance side of the package settled this question
already and chose Wilcoxon + Kruskal-Wallis. Using the same tests here means an abundance
result and a state result are directly comparable, and a reader learns one set of
caveats rather than two.

**Outputs:** `marker_state_stats.csv`, `marker_state.png`
(drawn by the package's own `fig_group_comparison()`, so the layout is
identical to the abundance figure by construction).

### 3.2 Phenotype heatmap — does the gate label match the phenotype?

**The gap.** A marker-by-population heatmap — markers in rows, populations in
columns, z-scored — is the standard annotation figure in this field. The gating
stage already produced the data; what was missing was the figure.

**Why it matters more here than in a clustering-first workflow.** There,
populations come from
unsupervised clustering and the heatmap's job is to *name* them. Here they come
from Boolean threshold gates written down in advance, so the heatmap's job is to
*check* them: if the column labelled "CD4 T cells" does not show high CD4 and low
CD8, the gate is wrong, and that is visible in one glance instead of being
inferred from a frequency that looks surprising. A pipeline that assigns labels
from thresholds needs a figure that can contradict the labels.

The heatmap **outlines** the cells each population's spec requires to be
positive, reading the same `populations:` block the scoring uses — so the audit
cannot drift from the definitions it audits.

Alongside it, a cohort-composition heatmap
normalises every cohort to the same notional cell count before taking each
population's split, so unequal sample numbers and acquisition depth cannot drive
the pattern.

**Outputs:** `population_marker_heatmap.png`, `cohort_composition_heatmap.png`

### 3.3 Compositional frequencies — percentages that must sum to 100

**The gap.** `abundance_measure()` states the problem exactly
right — *"frequencies are compositional: a fall in one population may reflect
expansion of another, not its own loss"* — and then nothing acts on it. Every
population's percentage was tested as if it could move independently, and it
cannot: a granulocyte expansion mechanically depresses every lymphocyte
percentage, and those depressions test significant with no lymphocyte having
changed.

**What was done.** The centred log-ratio (Aitchison 1986), with multiplicative
zero replacement (Martín-Fernández et al. 2003) so genuinely-absent rare
populations are not dropped — those being exactly the populations most likely to
differ.

**Both tests are reported, and the comparison is the finding.** A population
significant on raw percentage but not on CLR is a candidate artefact of the
constraint; one significant on both is robust to it. Swapping the test silently
would hide that distinction.

**What CLR does not fix**, stated plainly because a transform that looks like a
solution invites over-reading: it does not recover absolute cell numbers. If
every population doubles, the composition is unchanged and CLR sees nothing. Only
cells/µL separates "expanded" from "expanded relative to the rest" — which is why
`abundance_measure()` already prefers `cells_per_ul` when the patient table
supplies `wbc_per_ul`. **CLR is the right test for frequency data; absolute
counts remain the better data.**

**Outputs:** `compositional_clr_stats.csv`, `compositional_concordance.csv`
(verdict per population: `robust_to_composition`,
`raw_only__possible_composition_artefact`, `clr_only__was_masked_by_composition`,
`not_significant_either`)

### 3.4 Confounding — is a cohort difference really an age difference?

**The gap.** The patient table carries age and sex; nothing used them for
colouring and never for inference.

**Why the deliverable is a diagnostic rather than an adjustment.** At single-digit
n per cohort, covariate *adjustment* is close to unusable: a model with
`group + age + sex` spends most of its residual degrees of freedom on nuisance
terms, and if age is strongly associated with cohort — which in a
syndrome-vs-adult-control design it very often is — the two are not separable at
*any* n. The adjusted estimate is extrapolation, and it arrives wearing a
confidence interval.

So the order is deliberate:

1. **Always** — report whether each covariate differs between cohorts (Fisher /
   Kruskal-Wallis) *and* whether it tracks the outcome (Spearman /
   Kruskal-Wallis). Those two facts **together** are what makes something a
   confounder, and both are estimable at small n. Either alone is harmless, and
   flagging either alone would cry wolf on every study with unequal ages.
2. **On request** (`--rank-ancova`), and only where residual df permit — a
   rank-based ANCOVA, labelled `EXPLORATORY` in the output rows themselves. Where
   df are insufficient it writes `NOT FITTED` with the reason rather than a
   number that looks comparable to the others in the folder.

**Outputs:** `confounding_diagnostics.csv` (with a `confounder_risk` verdict),
`covariate_adjusted_stats.csv` (opt-in)

### 3.5 Paired and repeated-measures designs

**Why it must be declared and cannot be inferred.** Nothing in an FCS file or a
filename says two tubes came from the same donor at two timepoints. Treating
paired samples as independent — which an unpaired analysis necessarily does, having no
notion of pairing — discards the pairing's variance reduction *and* violates the
independence assumption the rank test rests on.

Runs only under `--paired-column` + `--condition-column`. Wilcoxon signed-rank,
Friedman for 3+ conditions, the conventional choice at this scale.
Replicate tubes of the same donor-timepoint are averaged before reshaping, and
incomplete pairs are **counted in the output**, not silently dropped — a pairing
that quietly halves n is exactly the failure this exists to prevent.

**Output:** `paired_comparison_stats.csv`

### 3.6 Batch effect — diagnosed, deliberately not corrected

**The gap.** The comparable toolboxes ship batch *correctors* — Harmony and
CytoNorm. cyRAVEN ships neither and, more to the point, had no way to tell
whether it needed one: an
island on the UMAP that is an acquisition batch is drawn identically to an island
that is a phenotype.

**Why this implements the diagnostic and not the correction.** Correction is
neither free nor neutral. Harmony moves cells so batch labels mix; when batch is
confounded with the biological group — samples of one cohort acquired on one set
of days, which is the normal way a clinical cohort is collected — *"removing the
batch"* and *"removing the effect you are looking for"* are the same operation,
and the algorithm cannot tell them apart. Running a corrector under that
confounding produces a clean-looking UMAP with the finding deleted.

So: measure the batch effect, measure the **confounding between batch and
cohort**, report both. If they are separable and mixing is poor, correction is
worth considering and the numbers say so. If they are confounded, no correction
is safe at any setting — and that is the single most important thing this
diagnostic can tell you.

**The metric** is iLISI (Korsunsky et al. 2019 — the metric Harmony itself is
evaluated with), against a **permutation null**: the achievable score depends on
how many batches there are and how unevenly sized they are, so "1.7" means
nothing on its own. Shuffling batch labels and recomputing gives the score this
dataset would produce with no batch structure at all. Confounding is reported as
Cramér's V with a plain-language verdict.

Falls back to the FCS `$DATE` keyword when `--batch-column` is not supplied —
a free batch variable, used only when it actually varies.

**Outputs:** `batch_mixing_stats.csv`, `batch_group_confounding.csv`,
`batch_diagnostic.png`

### 3.7 Threshold drift — does the *gate* move with cohort?

**Specific to this package; a clustering-first workflow has no equivalent**
because it never gates. One model is fitted to pooled data there, so every cell
is judged by the same rule. cyRAVEN derives a threshold **per sample** from that
sample's own density valley — which is what a human gater does, adapts correctly
to staining variation, and is the right default.

It also means "CD4-positive" is not literally the same predicate in every sample.
If the derived CD4 threshold sits systematically higher in one cohort, that
cohort shows fewer CD4 T cells **for that reason alone**, and the abundance test
reports it as biology with a small p-value. Nothing downstream can detect it,
because by then the thresholds have been applied and discarded.

The check is cheap and the data was already being written. A marker is flagged
only when the between-cohort gap is both statistically detectable *and* larger
than the within-cohort scatter — a gap of 0.3 asinh units is nothing if samples
within a cohort already vary by 0.5, and decisive if they vary by 0.02. Each
flagged marker names **which populations use it**.

**Outputs:** `threshold_drift_stats.csv`, `threshold_drift.png`

### 3.8 Unsupervised clustering — finding what the config does not describe

**The largest gap.** Every gated population is a Boolean conjunction of
thresholds written down in advance. That has real advantages — named,
reproducible, meaning exactly what the gating document says — and one decisive
weakness: **it cannot find anything the spec does not describe.** Cells matching
no definition become "Other CD45+", excluded from the UMAPs by default, so
unexpected biology is not merely unlabelled — it is not drawn.

It also has no way to be wrong out loud. A mis-set CD4 cut and a genuine absence
of CD4 T cells produce identical output. Unsupervised clustering breaks the tie:
a CD4 T-cell cluster appears at its true size whether or not the CD4 threshold
was set well, and the **disagreement** between cluster and gate is the diagnosis.

**Algorithm:** FlowSOM's two-stage design (Van Gassen et al. 2015 — winner of the
Weber & Robinson 2016 benchmark on accuracy and runtime): a self-organising map
quantises the marker space into ~100 micro-clusters, then those 100 prototypes are
metaclustered down to the requested count. The expensive step touches each cell
once per epoch; the clustering step operates on 100 vectors, not 10⁶.

Uses the **FlowSOM package when installed**; otherwise runs an equivalent
built-in batch SOM + Ward.D2 metaclustering, so no mandatory dependency is added.
The two are algorithmically equivalent and *not* numerically identical — which
one ran is recorded in the output, because a cluster number from one is not a
cluster number from the other.

**Clusters the marker matrix, not the UMAP coordinates.** UMAP is a
visualisation: its distances are not metric, it does not preserve density, and
the gaps between islands are partly an artefact of `min_dist`. Clustering it
would give clusters of the *picture*. FlowSOM and Phenograph both cluster the
high-dimensional space and use the embedding only to display.

**How to read the cross-check:**

| pattern | meaning |
|---|---|
| cluster >80% one gate label | gate and data agree |
| cluster mostly "Other CD45+" | **a real population the spec does not describe** |
| gate label split across clusters | the label is a union of distinct phenotypes |
| gate label holds few cells of a cluster it dominates | **the threshold is wrong** — the cells are there and cluster together; the Boolean rule is rejecting them |

That last row is the concrete motivation: if CD4 T cells score 0.29% of cells
while a ~20% cluster is CD4-bright and labelled "Other CD45+", the population is
not absent — the cut is misplaced, and the two numbers side by side say so
immediately.

**Outputs:** `unsupervised_clusters.csv`, `cluster_gate_agreement_clusters.csv`,
`cluster_gate_agreement_populations.csv`, `unsupervised_clusters.png`
(cluster panel beside gate-label panel — the comparison is the whole value)

#### 3.8.1 Learning a gate back out of a cluster — `--explain-clusters`

**The gap §3.8 opens and cannot close.** The agreement table can tell you that
cluster 7 matches no described population. It cannot tell you **what cluster 7
is**. The finding arrives as a cluster number and a cell count — nothing a
cytometrist can act on, because there is no way to go back to the instrument, or
to a FlowJo workspace, and select those cells.

`--explain-clusters` closes the loop in the other direction: for every cluster
the spec fails to describe, it **learns** a short sequence of two-marker gates
that selects it, and reports how well each one does.

**Why a convex polygon and not more thresholds.** Every gated population is a
conjunction of one-dimensional cuts, which makes each gate an axis-aligned
rectangle. Real boundaries are frequently diagonal — CD4/CD8, CD14/CD16,
FSC-A/FSC-H — and a rectangle placed over a diagonal boundary must either admit
contaminants or reject real cells. There is no third option. The package already
concedes this once: `derive_singlet_band()` is a hand-written diagonal band,
bolted on because a rectangle could not express it. A convex polygon is the
general form of that special case, and it is still a gate a human can draw.

**Algorithm.** A gate is the intersection of *k* half-planes; membership is
relaxed to `∏ σ(s·(wⱼ·x + bⱼ))` so it can be differentiated, and fitted by
class-weighted cross-entropy plus a tightness penalty whose strength is searched
against held-out F1 rather than fixed. Four of the eight normals are pinned to
±PC1/±PC2 of the target cloud — without that the feasible region need not be
**bounded**, and an optimiser will happily return a wedge running off to
infinity that scores well and is not a gate. Their offsets stay free, so the box
can still move and resize; the other four planes cut corners off it. Then the
next marker pair is chosen among the surviving cells and the whole thing repeats.

**No new dependencies, and no automatic differentiation.** The objective's
gradient is derived in closed form — four lines of matrix algebra — and the
problem is 24 parameters, which `stats::optim`'s L-BFGS-B solves directly. A
minibatch first-order method is the right tool when the parameter count is large
and the data does not fit in memory; here neither is true. `stats::prcomp`
initialises, `grDevices::chull` tightens, both already imported.

**Every metric is computed on held-out cells.** Eight free half-planes fitted to
a few thousand cells can memorise them, and a convex hull drawn around the
survivors memorises them completely. Fitting and scoring on the same cells
reports a number that is optimistic by construction, so cells are split before
the fit and the in-sample value is printed **beside** the held-out one — the gap
between them is information, and hiding it would be the whole problem.

> **This is descriptive.** A learned gate says where a group of cells sits in
> marker space. It is not evidence that the group is real, and it carries no
> p-value. Nothing it produces re-enters the analysis: the scored populations,
> the frequencies and every tested result are byte-identical whether this runs
> or not. It proposes; a human disposes.

Recall is measured against the **original** population size at every level, so a
strategy that discards half the cells at each step cannot report 100% recall
three times. The search stops as soon as another gate fails to earn its place.

**Outputs:** `cluster_gate_proposals.csv` (marker pair, held-out
precision/recall/F1, in-sample F1, cumulative metrics per level),
`cluster_gate_polygons.csv` (vertices in arcsinh units, ready to draw),
`cluster_gate_strategy_<k>.png` (one panel per gate, polygon over the cells it
was fitted to)

### 3.9 Data-driven subcluster count

`subcluster_by_reference()` fixed k = 3 for every population. Nothing checked
whether three compartments is what the data contains, and every subcluster label
in the multigraph overlay and every row of `subcluster_marker_shifts.csv` rests on
that number.

`--auto-subcluster-k` chooses k per population by mean silhouette on the
reference group's cells. **Honest limit**, because a "chosen" k invites more
trust than a fixed one: silhouette favours small, compact, roughly spherical
clusters — and so does k-means, so the criterion and the algorithm share a bias.
It is a better default than a hardcoded 3, not an oracle. The full score curve is
written out so the choice is inspectable.

**Opt-in**, because it changes the subcluster lettering, and published overlays
are indexed by that lettering.

**Output:** `subcluster_k_selection.csv`

### 3.10 UMAP model persistence

`run_umap()` discards uwot's model by default, so adding one sample re-embeds everything
and every coordinate moved. Cluster 4 in this run was not cluster 4 in the next,
and two results folders from the same study could not be laid side by side. The
README documents adding batches over time, which makes this a live problem.

`--save-umap-model` / `--umap-model` implement the established pattern: train
once with uwot's `ret_model`, then project later data into the same space.

**What projection is and is not:** projected cells are placed by the existing
model and do not move the manifold. That is the point — coordinates stay
comparable — and also the limitation: a population present only in the new samples
has no region of its own and will be placed among whatever it is nearest.
**Project when the cohort grows; retrain when the biology or the panel changes.**

The scaling constants travel *with* the model. Re-deriving median/MAD from a new
batch would scale it against itself and land it in a subtly different space — a
silent arithmetic error that looks exactly like a batch effect.

Uses `uwot::save_uwot()` rather than `saveRDS()`: a uwot model holds an external
pointer to its nearest-neighbour index, which `saveRDS` serialises as a null
pointer that fails on load, usually at transform time, well after the file looked
fine.

> **`--save-umap-model` changes the embedding — measured, not assumed.**
> Asking uwot to retain the model puts it on a different internal path that
> consumes the RNG stream differently. On 3,000 × 6 synthetic cells at a fixed
> seed, `ret_model = TRUE` and `ret_model = FALSE` are each perfectly
> reproducible with themselves but differ from each other: max coordinate
> difference 3.7, per-axis correlation 0.96. The embedding is equally valid — not
> degraded — but it is **not the same picture**, so clusters may be renumbered
> and islands may move. **Decide to persist the model at the start of a study,
> not after the figures have been circulated.** Without the flag, the default
> path is unchanged.



### 3.11 Pooled cofactor derivation

**Behaviour change — the one place a default number moves.**

The cofactor used to be derived from `reads[[sids[1]]]` — one sample
sets the transform for the whole panel, and therefore the scale on which every
threshold, marker median and UMAP distance is computed. If that file happens to
be weakly stained, an unstained control, or simply an outlier, the transform is
wrong for everyone else and nothing downstream can detect it. The choice of
"first" is alphabetical, so it is not even a considered sample — it is a filename.

Now: the **median** of up to 8 per-sample derivations, spaced evenly across the
sample order so a batch acquired in cohort order does not take its transform from
one cohort. Median, not mean, because the failure being guarded against is
precisely one aberrant sample.

If per-sample cofactors span more than 2×, that is logged: the panel does not
have one shared background, and a single cofactor is a compromise rather than a
description — usually a gain change mid-batch.

`--cofactor-from-first-sample` restores the previous behaviour exactly.

### 3.12 Run manifest

`sessionInfo()` used to be captured only inside `save_session()`, so a default
run left no record of the R version, package versions, command line, or input
files. A methods section needs all four.

`run_manifest.txt` is written **before** the expensive steps and rewritten on the
way out via `on.exit`, so it fires on error and interrupt too. A run that crashes
leaves `status: failed`; a manifest still reading `status: running` is itself the
diagnosis. Records R version, platform, every installed package version, the git
commit **and whether the checkout was dirty**, the full command line, every
option, every input file with size and mtime, and the derived cofactors.

## 4. Deliverables

### New files in `results/`

| File | From |
|---|---|
| `marker_state_stats.csv` | §3.1 |
| `marker_state.png` | §3.1 |
| `population_marker_heatmap.png` | §3.2 |
| `cohort_composition_heatmap.png` | §3.2 |
| `compositional_clr_stats.csv` | §3.3 |
| `compositional_concordance.csv` | §3.3 |
| `confounding_diagnostics.csv` | §3.4 |
| `covariate_adjusted_stats.csv` | §3.4 (`--rank-ancova`) |
| `paired_comparison_stats.csv` | §3.5 (`--paired-column`) |
| `batch_mixing_stats.csv` | §3.6 |
| `batch_group_confounding.csv` | §3.6 |
| `batch_diagnostic.png` | §3.6 |
| `threshold_drift_stats.csv` | §3.7 |
| `threshold_drift.png` | §3.7 |
| `unsupervised_clusters.csv` | §3.8 (`--cluster`) |
| `cluster_gate_agreement_clusters.csv` | §3.8 (`--cluster`) |
| `cluster_gate_agreement_populations.csv` | §3.8 (`--cluster`) |
| `unsupervised_clusters.png` | §3.8 (`--cluster`) |
| `cluster_gate_proposals.csv` | §3.8.1 (`--explain-clusters`) |
| `cluster_gate_polygons.csv` | §3.8.1 (`--explain-clusters`) |
| `cluster_gate_strategy_<k>.png` | §3.8.1 (`--explain-clusters`) |
| `subcluster_k_selection.csv` | §3.9 (`--auto-subcluster-k`) |
| `umap.model` + `umap.model.meta.rds` | §3.10 (`--save-umap-model`) |
| `run_manifest.txt` | §3.12 |

Every one of these is produced under `tryCatch`, after all primary outputs are
already on disk. A failure in any of them logs a warning naming the analysis and
leaves everything else untouched — the same contract the FlowJo export runs
under, and for the same reason: an addition that can destroy the thing it was
added to is not an addition.

### Changed deliverables

| File | Change | Breaking? |
|---|---|---|
| `population_marker_mfi.csv` | **+2 columns** `is_control`, `qc_status` | No — appended. Without them `qc_pass_rows()` passes the table through unchanged, so an unstained control's marker medians would enter the DS tests and the phenotype heatmap silently. Brings it in line with `functional_markers.csv`, which already carried both. |
| every derived number | cofactor now pooled across samples (§3.11) | Values may shift where samples disagree. `--cofactor-from-first-sample` reverts. |
| `config_cohorts.yaml` | **+3 keys** `heatmap_low/mid/high` | No — optional, defaults apply if absent. |

### Unchanged

Every pre-existing figure keeps its layout, colours, panel order, facet arrangement
and filename. Semantic colours still hold: healthy = green, subject 1 = red,
subject 2 = purple; male = green, female = red. The new diverging heatmap scale
is deliberately **not** drawn from `study_palette` — a colour that means "subject
1" in one figure must not mean "above average" in another.

## 5. Verification

Four levels.

**1. `R CMD check`.** Clean on the package's own container image.

**2. Unit — planted-truth fixtures.** `testthat`, run the usual way:

```r
devtools::test()          # or: R CMD check
```

The synthetic data carries known structure (a +1.5 asinh CD4 shift in one cohort
only, a genuinely depressed CD4 frequency, a CD4 threshold that drifts by cohort,
batch perfectly confounded with cohort) and the assertions check that each
function *finds that specific thing* — not merely that it returns a data frame.
Also checked: n is donors and not cells; controls and QC failures are excluded;
CLR sums to zero within sample; rank-biserial is bounded; RNG streams are
restored by every function that consumes them; figures actually write. Three
further tests guard the optimisations in §3.11 and the clustering and LISI inner
loops, asserting bit-identical agreement with the straightforward code they
replaced.

**3. End-to-end** on 25 real FCS files with `--cluster --rank-ancova`: completes,
producing all 18 new outputs plus every pre-existing output.

**4. Equivalence against the pre-refactor implementation — the important one.**
Run with `--cofactor-from-first-sample` against the same inputs:

> **25 of 26 outputs byte-for-byte identical.** The single difference is
> `population_marker_mfi.csv`, and on its original 7 columns it is also
> byte-identical — the diff is exactly the two documented added columns. Same
> cofactor (75.3), same QC verdicts, same figures.

So every behavioural difference in normal operation is attributable to the
pooled cofactor (§3.11) and nothing else.

### One bug this verification caught

The end-to-end run initially failed with all 25 samples reporting *"no separable
CD45+ mode"* and every threshold `Inf`. Cause: `$` on an R list partial-matches,
and optparse **omits options whose default is `NULL`** from the list it returns.
The existing `opt$cofactor` (default `NULL`, therefore absent) silently
partial-matched the newly added `cofactor_from_first_sample` and returned
`FALSE`. `FALSE` is not `NULL`, so the pipeline accepted it as a user-supplied
cofactor of 0, making every marker `asinh(x/0) = Inf`.

The failure was total, silent, and blamed the data. Fixed three ways: the dest no
longer begins with `cofactor`; the call site uses `opt[["cofactor", exact =
TRUE]]`; and a non-positive cofactor now stops the run naming the value. The
code had met this trap once before — see `opt[["umap_markers", exact =
TRUE]]` comment — and an audit of every option name against every `opt$` access
found no other instance.

Unit tests could not have caught this: it lives entirely in the CLI wiring.

## 6. Not implemented, and why

### 6.1 Batch *correction* (Harmony / CytoNorm)
Deliberate, not deferred — see §3.6. Under batch–cohort confounding, correction
and deleting the finding are the same operation. The diagnostic tells you which
situation you are in; that is the honest deliverable at this cohort design.

### 6.2 Phenograph clustering
FlowSOM covers the same need and is faster at these cell counts. Phenograph
(Louvain on a k-NN graph) picks its own cluster count, which is genuinely useful,
but requires `Rphenograph`/`Rphenoannoy` — both GitHub-only, neither on CRAN or
Bioconductor, which would break the container's reproducibility guarantee.
**Add if:** the cluster count itself becomes a question rather than a setting.

### 6.3 diffcyt (limma / LMM differential testing)
The rank tests in §3.1 and §3.3 answer the same questions without the
distributional assumptions, which is the right trade at single-digit n — see
§3.1. **Add if:** cohorts reach ~15+ samples each, where the empirical-Bayes
prior starts paying for itself and mixed models can absorb repeated measures
properly.

### 6.4 Per-marker cofactors
`derive_cofactor()` now *returns* per-marker candidates (`attr(cf,
"per_marker")`), but applying them is not wired up. Doing so changes the arcsinh
scale per channel, which shifts every threshold, every marker median and every
UMAP distance — a large methodological change that cannot be validated without
re-running and re-reading every gating QC figure. The data to make the decision
is now collected. **Add if:** the gating QC shows channels whose backgrounds
genuinely differ in scale.

### 6.5 Pseudotime / trajectory analysis
The comparable toolboxes wrap diffusion maps and Slingshot. Trajectories
assume a continuous
developmental process; this panel measures terminally-differentiated peripheral
blood populations, where a fitted trajectory would be an artefact of the fitting.
**Add if:** the panel gains a developmental axis.

### 6.6 Supervised cell-type classifiers (Astir, CytoDx, random forest)
Label-transfer models are a standard offering elsewhere. Here the labels come
from a written gate spec, so a classifier trained on them can only reproduce the
spec — including its errors. §3.8 is the better cross-check because it is
independent of the labels. **Add if:** an external reference-annotated dataset
becomes available to train against.

### 6.7 FlowJo `.wsp` workspace ingestion
Reading gates back out of a FlowJo workspace is the inverse of what cyRAVEN
does (`--flowjo-export`). Parsing `.wsp` needs `CytoML` +
`flowWorkspace`, a heavy dependency pair. **Add if:** manual FlowJo gates need to
be the input rather than the output.

### 6.8 Cell-level differential expression
Deliberately absent, for exactly the pseudoreplication reason in §3.1 — which
is why the tools that do provide a cell-level Wilcoxon ship it with a warning
against using it. The cell-level view
is already available descriptively in `subcluster_marker_shifts.csv`, which
carries effect sizes and **no p-values** — that is the correct treatment.

## 7. Statistical conventions

Applied consistently across every test the package runs:

- **Replicates are samples, never cells.** Anything computed over pooled cells is
  labelled descriptive and carries no p-value.
- **Rank-based tests throughout** — Wilcoxon rank-sum, Kruskal-Wallis, Spearman,
  and their paired counterparts. No normality assumption at these n.
- **Benjamini-Hochberg** within each test family. DS testing adjusts *within*
  each measure, not across both: `median_asinh` and `pct_positive` are two views
  of the same cells and are strongly correlated, so pooling them would inflate the
  family with non-independent tests and over-penalise every hit.
- **Effect sizes reported alongside every p-value** — Cliff's delta over samples,
  rank-biserial for paired.
- **Controls and QC failures excluded before testing**, via `qc_pass_rows()`.
- **Differences, not ratios**, for arcsinh values: the scale runs negative, where
  a ratio is uninterpretable and changes sign for no biological reason.

## 8. References

- Van Gassen S. *et al.* (2015) FlowSOM. *Cytometry A* 87:636–645.
- Weber L.M. *et al.* (2019) diffcyt. *Commun Biol* 2:183.
- Weber L.M. & Robinson M.D. (2016) Comparison of clustering methods for
  high-dimensional cytometry data. *Cytometry A* 89:1084–1096.
- Korsunsky I. *et al.* (2019) Harmony. *Nat Methods* 16:1289–1296.
- Van Gassen S. *et al.* (2019) CytoNorm. *Cytometry A* 97:268–278.
- Aitchison J. (1986) *The Statistical Analysis of Compositional Data.*
- Martín-Fernández J.A. *et al.* (2003) Dealing with zeros. *Math Geol*
  35:253–278.
- McInnes L. *et al.* (2018) UMAP. arXiv:1802.03426.

## 9. Citing cyRAVEN

```r
citation("cyRAVEN")
```

The methods the package implements belong to their original authors; please cite
them as well where relevant. `citation()` lists them.

## 10. Licence

GPL-3.0-or-later. The full text is in [LICENSE.md](LICENSE.md).

This is a deliberate choice rather than a default. cyRAVEN links against `uwot`,
which is GPL-3; a permissive licence on top of that would misdescribe what a user
of the *installed* package may actually do. Copyleft also matches what the code
is for — a result you cannot reproduce is not a result, and requiring that
modifications stay open keeps a published figure traceable to source that is
still available to whoever tries to check it.
