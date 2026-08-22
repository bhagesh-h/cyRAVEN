# Known limitations, and what was deliberately left out

Every caveat below exists somewhere else in this repository already, in
a code comment, a NEWS entry or a function’s documentation. Scattered
like that, a new user has no way to see the aggregate risk before
deciding whether to trust a result. This page is the single place that
collects them.

The same page also records what was considered and **not** added, and
why. That half exists so the same arguments do not get relitigated every
six months.

## Part 1. Limitations

### A control used to be believed too easily

[`staining_verdict()`](https://bhagesh-h.github.io/cyRAVEN/reference/staining_verdict.md)
read an **absent** `is_control` column as “might be a control”, because
only `is_control = FALSE` positively asserts that a file is a biological
sample. On a sheet without that column, any sample whose CD45 gate found
no density minimum was promoted to the unstained reference for its whole
panel, and every threshold there became the 99.5th percentile of that
one sample. Where it was really a stained sample, every population
collapsed, with no error raised and no warning printed.

**How to check whether a result you already have was affected.** Open
`thresholds_used.csv`. If `source` reads `control_q995` on most markers,
those thresholds came from a reference. Find out which sample supplied
it before reading a single frequency. A cohort that declares no control
should show `valley` and `quantile_fallback` only.

Fixed. `--discover-controls` restores the old behaviour for anyone who
needs it.

**The version number will not tell you whether you are affected, so
check the file.** A `1.0.0` image was published before this fix and
replaced afterwards under the same tag, so two pulls of
`bhagesh/cyraven:1.0.0` can differ. The `source` column in
`thresholds_used.csv` is the reliable test, and it costs a few seconds.

### Automated gates have not been benchmarked against independent expert gating

[`learn_convex_gate()`](https://bhagesh-h.github.io/cyRAVEN/reference/learn_convex_gate.md)
and
[`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)
are validated on held-out cells, and, where donor count allows,
leave-one-donor-out. Both of those are internal: the only cohort they
have been measured on is the demonstration one. Nothing here establishes
that the defaults generalise to a different panel, a different
instrument, or a different laboratory’s staining.

Two consequences worth stating plainly. The hardcoded constants in
[`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md),
`min_rel_depth` among them, may have been tuned to one cohort without
anyone intending to tune them. And a learned gate is descriptive: it
carries no p-value and is not evidence that the population is real.

**Read the minimum, not the median.** On a 21-donor cohort, three
clusters gave median F1 of 0.987, 0.731 and 0.189 with worst-donor F1 of
0.953, 0.111 and 0.008. Cells held out from a fit come from the same
donors, acquired in the same tubes on the same day, so they flatter the
model. Only the per-donor minimum tells you whether a gate survives
contact with the next patient.

### No mixed-effects differential abundance

The tests are rank-based: Wilcoxon rank-sum, Kruskal-Wallis, Wilcoxon
signed-rank and Friedman when paired, with Benjamini-Hochberg within
families. Parametric equivalents sit beside them with their assumptions
attached.

None of those can model donor and batch as **crossed random effects**.
They can stratify or pair, which is not the same thing. A GLMM with a
count offset is the standard remedy, and `statistical_methods.csv` names
it in the catalogue without computing it. That gap is deliberate for now
rather than overlooked, and it is recorded in the backlog.

### A single convex gate can misrepresent a multimodal population

[`learn_convex_gate()`](https://bhagesh-h.github.io/cyRAVEN/reference/learn_convex_gate.md)
fits one convex polygon per population. Nothing currently checks whether
the label being fitted is a single coherent region in the chosen marker
pair before fitting it. A hull drawn around two disjoint clusters of the
same nominal population reports a large, low-precision polygon that
“explains” the label without describing the geometry.

If a population is known to be phenotypically heterogeneous, check
`spec_gaps.csv` from explore mode, which reports how many effective
clusters a declared label spans.

### FMO disagreement is computed but not enforced

`fmo_agreement.csv` reports the distance between a derived cut and its
FMO-anchored equivalent in units of that threshold’s own uncertainty.
Beyond about three units, one of them is wrong. Nothing acts on that:
the file is written and the run proceeds. An unattended pipeline will
not stop on it, unlike `--fail-on-drift` for baseline conformance.

### Spectral panels can fragment a cohort silently

The panel fingerprint is the set of resolved marker names, so a channel
present in some files and absent in others makes them separate panels,
each with its own cofactor, thresholds and embedding. Spectral unmixing
writes extracted autofluorescence back as extra channels and how many
appear varies per acquisition, so twelve comparable files can become
seven panels of one to three files each.

Nothing fails when this happens. Every file loads and every table is
written. `--list-channels` now names the offending channels and prints
the exact `--ignore-channels` pattern that merges the cohort, so run it
before a long run.

There is no unmixing-specific QC. `spreading.R` diagnoses conventional
spillover spreading well, but poor reference spectra and
autofluorescence misattribution are distinct failure modes it does not
capture.

### Everything is documented on one cohort

The README, the worked example and the container smoke test all run on
the same demonstration cohort: one panel, one biological question.
Behaviour on a larger marker panel, on spectral channels, or on a design
with three or more groups is exercised by tests but not by a documented
worked example.

### Batch correction changes intensities, not conclusions, and can be wrong

`--correct-batch` aligns marker distributions and refuses above a
Cramér’s *V* threshold between batch and group, because at that overlap
removing the batch and removing the finding are the same operation.
`--batch-method cluster` fits one map per cell type rather than one per
file, which is the better choice when a detector shift moves bright and
dim populations by different amounts.

Two things to know. The clustering it fits is derived from the data, so
a cohort with genuinely batch-specific biology will have that biology
partly absorbed into the correction. And no correction is applied to the
gate hierarchy, which is derived per sample and is already batch-local;
correction is descriptive of intensity only.

## Part 2. Deliberately not included

Nothing below is missing by accident. Each was considered against what
the package already does and rejected for the stated reason. If
circumstances change, the reason is what has to change first.

### Packages

| Package | Overlaps | Why not |
|----|----|----|
| **AutoSpill** | compensation | Installs only from GitHub, which breaks the dated-snapshot guarantee the Docker image rests on. For spectral data, unmixing has already happened on the instrument |
| **CytoNorm** | batch normalisation | Same GitHub-only problem, and its API is file-based: `QuantileNorm.train()` takes FCS paths and `QuantileNorm.normalize()` writes new FCS to disk, while correction here happens on an in-memory matrix already read, transformed and gated. Its *method*, one map per cell type, is implemented natively as `--batch-method cluster` |
| **fdaNorm**, **gaussNorm** | normalisation | No longer supported, and `fdaNorm` is gone from recent flowStats versions |
| **flowTrans** | transformation | The cofactor is derived per panel from the data, so there is nothing to fit. The package also carries a multi-year unfixed biexponential bug |
| **ggcyto**, **flowViz** | plotting | One figure system and one theme already. A second plotting idiom costs consistency and adds nothing a reader sees |
| **Rtsne** | dimensionality reduction | UMAP only, deliberately. Two embeddings of the same cells invite comparison between two things that are not comparable |
| **diffcyt**, **CATALYST** | differential abundance | They test clusters *defined on the subsample*, so the subsample changes what is tested. That is the failure the donor-level design exists to avoid. Note this rejects the paradigm, not the mathematics: a GLMM on our own per-population counts remains open |
| **FastPG**, **Rphenograph** | clustering | FlowSOM is a cross-check here, not the primary method, so clustering speed is not the binding constraint |
| **openCyto** | gating | CSV templates applied uniformly across samples. The premise of this package is that one gate copied across samples is the problem |
| **flowAI**, **flowClean** | acquisition QC | They remove anomalous events. Acquisition QC here reports what removing them *would* change and deletes nothing unless asked. Possible optional second opinion later |
| **PICAFlow**, **CytoPipeline**, **CytoExploreR** | whole workflow | Alternatives to this package, not components of it |
| **CytoML**, **flowWorkspace** | FlowJo import, Gating-ML | **Not rejected, deferred.** Reading a `.wsp` into `--external-labels` would let a hand-drawn strategy be scored leave-one-donor-out. Held back on dependency weight: it pulls `cytolib`, `RProtoBufLib` and `ncdfFlow`, which compile C++ |

### Capabilities

| Capability | Why not |
|----|----|
| **Trajectory inference / pseudotime** | Ordering cells along a trajectory assumes the population is a continuum that was sampled densely enough to reconstruct. A gated cytometry panel with a dozen markers usually is not, and a trajectory will still be returned |
| **Batch correction without a confounding check** | Available in several other tools. Where batch and group overlap far enough, correcting removes the finding, so the check comes first and the refusal is the feature |
| **Automatic removal of anomalous events** | The impact of removal is reported first. `--drop-unstable-events` performs it only when asked |
| **A p-value on anything that cannot be collapsed to one value per sample** | Testing over events treats one deeply acquired donor as many independent observations. Those quantities get an effect size and no p-value |
| **Trained classifiers for label transfer to a new cohort** | A classifier transfers a labelling; it does not transfer the reason for it. What moves between cohorts here is a gating strategy, which is inspectable, exportable to the instrument, and can be scored per withheld donor. A model that cannot be read is a worse artefact even when it is more accurate |
| **Cross-panel benchmarking against public expert gates** | Wanted, and blocked on a rule that is not worth breaking: third-party FCS does not go into this repository. It needs an external harness, not a vignette |
| **A GUI** | The point of the package is that a run is one command with a recorded invocation. A GUI makes the choices unrecorded again, which is the problem it was built to solve |

## Where to look when something seems wrong

| Symptom | File to open first |
|----|----|
| Frequencies near zero across the board | `thresholds_used.csv`, `source` column |
| A population is missing entirely | `staining_qc.csv`, then `gating_qc.png` |
| More panels than expected | `--list-channels`, which now names the fix |
| A comparison you expected is absent | `design_feasibility.csv` |
| A difference that looks too good | `threshold_drift_stats.csv`, then `difference_over_gate_u` |
| Correction ran or did not | `batch_correction.csv` |
