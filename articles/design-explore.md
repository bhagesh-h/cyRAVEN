# Why explore mode is built the way it is

The design record for explore mode: what the supervised and unsupervised
halves can lend each other, which of those exchanges are implemented,
which are deliberately refused, and where that leaves cyRAVEN against
the tools it sits beside.

Written because “add unsupervised clustering” is the easy half of the
problem. The hard half is deciding what the two halves are allowed to
tell each other without destroying the guarantee that made either worth
having.

## 1. The two failure modes

Supervised and unsupervised cytometry analysis fail in opposite
directions, and neither can detect its own failure.

**Supervised analysis cannot find what it was not told to look for.**
You declare twelve populations, it scores twelve populations. A
thirteenth cell type sits in the data unremarked. Worse, a declared
label may cover two distinct phenotypes, and their total can stay flat
while one rises and the other falls. That is a real pattern that no
amount of testing on the declared table can surface, because the table
has one column where the biology has two.

**Unsupervised analysis cannot tell you what it found.** FlowSOM returns
cluster 17. Whether cluster 17 is a B cell subset, a doublet artefact,
or the same cells as cluster 4 under a different sampling seed is left
to a human squinting at a heatmap. And because the clusters are defined
on the data, the number of them that reach significance moves with
parameters nobody reports. Measured elsewhere in this project, one
cyCONDOR run went from **1 significant cluster to 8** by changing only
how many cells were drawn.

Running both and reading them side by side does not fix either. What
fixes them is deciding, deliberately, which specific facts may cross
between them.

## 2. The exchange, enumerated

Every transfer that is technically possible between the two halves, with
the verdict on each.

### 2.1 Declared -\> explore

| What could be lent | Verdict | Why |
|----|----|----|
| **The transform and its estimated cofactor** | **Always** | A property of the panel derived from the data, not a product of the declaration. An unsupervised tool assumes a transform; cyRAVEN estimates one. Lending it costs nothing and removes an arbitrary constant |
| **Per-sample thresholds** | **Under `--maybe-learn`** | The single largest gain. See section 3 |
| **Staining QC verdicts** | **Under `--maybe-learn`, as a column** | Recorded, never used to exclude. See section 4 |
| **Batch/group confounding verdict** | **Under `--maybe-learn`** | Without it explore reintroduces the exact failure the declared path refuses. See section 5 |
| **Group labels** | **Always** | They come from the sample sheet, not from the analysis. Withholding them would mean explore could not test anything |
| **The parent gate** | **Never** | Explore exists to look *outside* it. Applying it would make explore a re-clustering of the declared result |
| **The population definitions** | **Never** | That is the thing being checked |

### 2.2 Explore -\> declared

| What could be lent | Verdict | Why |
|----|----|----|
| **Populations spanning several clusters** | **Under `--maybe-learn`**, as `spec_gaps.csv` | The finding the declared path structurally cannot make about itself |
| **Clusters no population covers** | **Under `--maybe-learn`**, same file | Names what the specification missed |
| **Cluster-derived population definitions, adopted automatically** | **Never** | Circular. See section 6 |
| **Re-gating, re-thresholding, or altering any declared frequency** | **Never** | A run whose numbers were changed by a clustering is no longer the run its manifest describes |

The asymmetry is deliberate. Explore may **describe** the declared
analysis. It may never **alter** it.

## 3. Per-sample thresholds: the exchange that matters most

A standalone clusterer names a cluster by pooling every cell, computing
a median per marker, and putting a colour scale on it. Three things are
wrong with that on multi-sample data, and all three are things cyRAVEN
has already solved for the declared path:

1.  **The pooled median is not a threshold.** It is the middle of the
    distribution, which for a marker where 80% of cells are positive
    sits inside the positive population.
2.  **It ignores staining variation between samples.** A bright sample
    and a dim one are compared against one number, which is precisely
    the systematic bias that per-sample gating exists to remove.
3.  **It produces a colour, not a call.** “Reddish for CD19” is not a
    measurement.

Under `--maybe-learn`, explore reads `thresholds_used.csv` and calls
each cell against **its own sample’s cut**. A cluster’s phenotype
becomes:

    CD19+ HLA-DR+ CD3- CD14-

where every `+` means *this fraction of the cluster’s cells were above
the threshold derived in their own sample*, recorded in
`frac_pos.<marker>` columns. A marker the cluster is genuinely mixed for
is **omitted** rather than called either way, so the phenotype string
admits what it does not know.

That is a claim no standalone clusterer can make, and it is not a clever
algorithm. It is the reuse of work the declared path already did.

## 4. Two things explore can do that the declared path cannot

**It keeps the samples staining QC excluded.** A sample with no
resolvable CD45 mode has no usable parent gate, so every declared
percentage from it would be a fraction of an arbitrary slice of the
scatter plot. The declared path must drop it. Explore does not depend on
that gate, so those samples still contribute, and `staining_qc_verdict`
is a *column* rather than a filter. On the KMT2 cohort analysed in this
project that is three donors the declared path had to discard and
explore can still describe.

**It looks outside the parent gate.** Anything the CD45 gate excluded is
invisible to the declared analysis by construction. Explore embeds it,
which is the only way a population that fails the parent gate can ever
be seen.

## 5. Carrying the confounding verdict

The most likely way explore mode could damage this package is by
becoming a p-value generator that ignores everything the declared path
established.

`explore_cluster_stats.csv` therefore carries `batch_group_cramers_v`
and `batch_group_verdict` as columns of the results table itself, not as
a note elsewhere. On a cohort where batch and group are confounded, a
*q* \< 0.05 on a cluster means exactly as little as a *q* \< 0.05 on a
declared population, and the table says so on the same row.

Explore’s statistics are also donor-level, using the same
[`stats_group_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_group_comparison.md)
the declared path uses: one value per sample, Kruskal-Wallis omnibus,
Wilcoxon pairwise, Cliff’s delta, Benjamini-Hochberg. Not a
reimplementation: the same function, so the two cannot drift apart.

## 6. Why the draft specification is never adopted automatically

`explore_suggested_spec.yaml` turns clusters into declarable population
definitions. The obvious next step, feeding it straight back into the
declared pipeline, is refused, and this is the most important design
decision in the feature.

cyRAVEN’s argument is that populations are declared **before** the data
are examined, and six diagnostics then try to falsify that declaration.
Derive the specification from the same data and the diagnostics become
circular:

- **Phenotype concordance** asks whether a population expresses the
  markers its definition demands. Trivially yes, if the definition was
  read off those cells.
- **Gate-cluster concordance** compares declared gates against
  unsupervised clusters. If the gates came from those clusters, it
  compares a thing to itself.

There is a sharper version. Explore’s clustering is **group-blind**. It
never sees cohort labels, so deriving a definition from cluster
*phenotypes* is a mild form of double-dipping. But
`explore_cluster_stats.csv` **does** test by group. A user who reads
that table, sees which clusters differ, and declares those as
populations has selected on the outcome, and every subsequent p-value is
invalid.

The defensible uses, ranked:

| Workflow | Verdict |
|----|----|
| Explore on cohort A -\> curate -\> declared run on **cohort B** | Clean. The intended loop |
| Explore -\> curate -\> declared run on the **same** cohort, reporting frequencies | Acceptable; state that the spec was data-derived |
| Explore -\> read the group statistics -\> declare the clusters that differed | **Invalid.** Selection on the outcome |

Hence: a YAML file whose header says it is a suggestion, requiring a
human edit, with no flag anywhere that adopts it.

## 7. Choices inherited rather than invented

Three decisions were taken from the unsupervised literature because they
are better than what a supervised tool would reach for.

**Cluster everything, including scatter and viability.** No lineage
preference. This is what lets the QC gate find debris (leukocyte marker
low), dead cells (viability high) and granulocytes, which no antibody in
most panels identifies. The cost is real: some clusters split on scatter
rather than lineage. `--explore-markers` overrides it.

**Gate on clusters, not events.** Cluster coarsely first, then judge
*whole clusters* by their median profile, so no per-event cutoff is
invented. The debris and dead calls come from 2-means on the per-cluster
medians when nothing better is available, and from per-sample thresholds
when it is. Both are recorded in the `basis` column of
`explore_qc_clusters.csv`.

**Refuse a gate that keeps almost nothing.** Below
`--explore-min-retained` (default 5%), the gate is ignored and the run
flagged: at that point the markers are likelier mis-named than the data
bad.

One decision goes the other way. **Cells are re-equalised per sample
after the gate.** cyCONDOR does not do this, and its faceted panels then
differ in density for a reason that is group size rather than biology.
That trap was documented while running it during this project.
Equalising costs cells and buys panels that can honestly be read side by
side.

## 8. Where this leaves cyRAVEN

| Capability | cyRAVEN | FlowSOM / Phenograph | cyCONDOR | diffcyt | CATALYST |
|----|:--:|:--:|:--:|:--:|:--:|
| Per-sample gate derivation | yes |  |  |  |  |
| Uncertainty propagated to every frequency | yes |  |  |  |  |
| Detection limits per population | yes |  |  |  |  |
| Donor-level inference by construction | yes |  | part | yes | part |
| Refuses batch correction when confounded | yes |  |  |  |  |
| Unsupervised discovery | yes | yes | yes | yes | yes |
| Clusters named from per-sample thresholds | yes |  |  |  |  |
| Cluster -\> executable gate geometry | yes |  |  |  |  |
| Gating-ML 2.0 export | yes |  |  |  |  |
| Empirical-Bayes moderation across clusters |  |  |  | yes |  |
| Trajectory inference |  |  | yes |  |  |
| Data-determined cluster count |  | yes | yes |  | yes |

`yes` present. `part` possible but not enforced by the tool.

The gap this feature closes is the “unsupervised discovery” row. The row
that remains distinctive is **clusters named from per-sample
thresholds**, which requires both halves in one tool and is why bolting
FlowSOM onto a supervised pipeline is not the same thing.

### What is deliberately not adopted

**diffcyt’s empirical-Bayes moderation** (edgeR, limma-voom, GLMM). It
borrows strength across clusters via a shared prior, which pays off at
mass-cytometry marker counts with hundreds of clusters. On a declared
panel of a dozen populations there is little to borrow, and the
moderation is harder to explain than it is worth. Reachable through
`--external-labels` for anyone who wants it.

**A data-determined cluster count.** Phenograph and consensus methods
choose *k* from the data. `--explore-k` is a setting, and
`explore_provenance.csv` says so, because a *k* chosen from the data is
one more thing the significance count depends on.

**Trajectory inference.** Diffusion maps and pseudotime need a
differentiation continuum the design supports. On cross-sectional human
cohorts of 10 to 40 samples it produces a picture, not a measurement.

## 9. The claim worth testing

If there is a methods paper here, the claim is not “cyRAVEN clusters
cells”. It is:

> Threshold placement uncertainty, derived per sample and propagated to
> reported frequencies, predicts the between-operator variance that
> manual gating produces, and that populations whose propagated
> uncertainty is large are the ones where published cohorts disagree.

That is falsifiable and, as far as this project’s reading goes,
untested. It would need the same files gated by several operators, with
the propagated uncertainty compared against the observed spread of their
results. The inter-laboratory data to do it exists in the literature
this package already cites.

The second claim, weaker but easier to support: **a refusal rule is a
feature.** A tool that declines to run batch correction above a measured
Cramér’s *V* between batch and group, and says why, produces fewer
results and more defensible ones. That one needs a demonstration rather
than a benchmark: a cohort where correction would have produced a clean,
wrong answer.

## See also

- [Explore
  mode](https://bhagesh-h.github.io/cyRAVEN/articles/explore.md): how to
  run it
- [Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.md):
  why the donor is the unit of replication
- [Interoperability](https://bhagesh-h.github.io/cyRAVEN/articles/with-cycondor.md):
  handing a cyCONDOR clustering back
- [Scope](https://bhagesh-h.github.io/cyRAVEN/articles/scope.md): nine
  excluded methods and the conditions for each
