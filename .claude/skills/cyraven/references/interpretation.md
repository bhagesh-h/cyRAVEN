# Reading a cyRAVEN result

A default run writes 25 tables and 22 figures. The order below is not arbitrary —
each stage can invalidate everything after it, so working backwards from the
p-values wastes time on results that were already dead.

Sections 1 to 3 reject most compromised runs. Read them first, every time.

---

## 1. Did the gates land in the right place?

**`recon_diagnostics.png`, `gating_qc.png`** — written unconditionally. Per
sample: the scatter gate boundary, the singlet band, and every marker threshold
drawn on the density it came from.

A threshold sitting on a distribution shoulder rather than a minimum is visible
here and **in no downstream table**. Nothing else in the run can tell you.

Stop here if the gates are wrong. Every number below inherits them.

## 2. Is the staining usable?

**`staining_qc.csv`** — one verdict per sample. A sample with no resolvable CD45⁺
mode has no usable parent gate. It is excluded and the exclusion reported, rather
than contributing a spurious frequency to a group mean.

Three states that are easy to confuse:

| Column | Meaning |
|---|---|
| `is_control = TRUE` | An actual control tube. Excluded from testing, not a failure |
| `qc_status = "failed"` | A declared sample that failed staining QC. Its percentages carry no evidence |
| neither | An ordinary sample |

`--include-qc-failed` forces failures back in for inspection. Their
`pct_of_cd45_pos` still carries no evidence; the verdict column says which were
forced.

`--min-cd45-pct` sets the floor, default 5%.

## 3. Do the populations express what they were declared to express?

**`population_marker_heatmap.png`** — marker intensity against declared
populations, z-scored across the run.

In clustering-first analysis this figure *assigns* identity. Here identity is
already declared, so it does the inverse: a column labelled CD4 T cells that does
not show elevated CD4 and depressed CD8 **falsifies the gate that produced it**.

The figure outlines the markers each population's spec requires to be positive,
reading the same `populations:` block used for scoring, so the audit cannot
diverge from the definitions being audited.

**`cohort_composition_heatmap.png`** normalises each group to a common notional
cell count before partitioning by population, so unequal sample numbers or
acquisition depth cannot generate the pattern.

---

## 4. The four falsification checks

### 4.1 Threshold drift

**`threshold_drift_stats.csv`**, **`threshold_drift.png`**

Per-sample thresholding carries a specific hazard: if thresholds for a marker
differ systematically between study groups, the populations defined by that
marker are **not identically specified across the comparison**, and any abundance
difference is partly definitional rather than biological.

A flagged marker invalidates naive interpretation of every population that reads
it. No clustering method has an equivalent output, because none carries a prior
definition capable of drifting.

### 4.2 Gate uncertainty

**`threshold_uncertainty.csv`**, **`uncertainty_budget.csv`**,
**`frequency_uncertainty.png`**, **`uncertainty_budget.png`**

Threshold drift asks whether a cut sits in a different place for one group than
another. This asks the prior question: how well is it determined at all.

Each cut is re-derived from resamples of the events it was computed on
(`u_sampling`), and again across the settings that placed it — histogram
resolution, smoothing, how completely two modes must separate (`u_method`). The
two combine in quadrature per the GUM convention.

**Read `bootstrap_valley_rate` before `u_combined`.** A cut recovered in every
resample marks a real boundary. One recovered in half of them was found by
histogram noise that happened to clear the depth rule — and both are written to
`thresholds_used.csv` in exactly the same way. Below about 0.8, treat the
threshold as unresolved rather than imprecise.

`uncertainty_budget.csv` says which threshold each population's uncertainty comes
from. The CD45 parent gate is a term in every population, because displacing it
moves cells into and out of the denominator all of them are expressed against.
Expect it to be the largest contributor; that is what the operator-variability
literature reports for manually gated data too.

`frequency_uncertainty.png` puts the two quantities side by side: the spread
between samples, and the bar within each. **Where the bars are as wide as the
scatter, the variation on display is the cut moving, not the biology.**

None of this alters a threshold. The value in `thresholds_used.csv` is
`density_valley()` at its defaults on the real events either way; the
perturbation runs on copies. `--no-uncertainty` reproduces the previous output
byte for byte.

The comparison worth making: operator studies following the same convention
report expanded uncertainty rising from about 12% on a three-gate strategy to
about 16% on a five-gate one, with the first gate contributing most (Whitmore et
al., 2021, *Methods Protoc* 4:24).

### 4.3 Cluster concordance

**`cluster_gate_agreement_clusters.csv`**,
**`cluster_gate_agreement_populations.csv`** — requires `--cluster`.

A self-organising map clustering of the same cells, computed without reference to
the specification, then cross-tabulated. Four configurations are diagnostic:

| Pattern | Reading |
|---|---|
| Cluster mostly one declared label | Concordant. Spec and data agree |
| Cluster dominated by the unassigned label | **A population the spec does not describe.** Supervised analysis cannot produce this finding on its own |
| One label spread across several clusters | Definition coarser than the phenotypic structure. If the subsets move in opposite directions, the aggregate reports no change |
| Label holds far fewer cells than the cluster it dominates | **Misplaced threshold.** If CD4 T cells score 0.3% while a 20% cluster is CD4-bright and unassigned, the population is there and the Boolean rule is rejecting it |

The last one is the diagnosis a frequency table structurally cannot deliver.

### 4.4 Gate geometry and transferability

**`cluster_gate_proposals.csv`**, **`cluster_gate_polygons.csv`**,
**`cluster_gate_strategy_<k>.png`** — requires `--explain-clusters`.

Converts "cluster 7 matches nothing" into gate geometry a cytometrist can draw:
the two most discriminating markers, a convex polygon in that plane, then
recurse. Convex polygons rather than rectangles because CD4/CD8, CD14/CD16 and
FSC-A/FSC-H all separate along oblique boundaries, and a rectangle on an oblique
boundary must either admit contaminants or reject real cells.

Every metric is computed on events held out of the fit, with the resubstitution
value printed beside it. Eight free half-planes can memorise a few thousand
events, so the gap between the two quantifies overfitting.

Descriptive only. No p-value; scored frequencies and test results are identical
whether it runs or not.

**`gate_transferability.csv`**, **`gate_transferability_summary.csv`** — requires
`--external-labels`.

Refits the whole strategy with one donor withheld and scores it on that donor,
once per donor.

**Read `f1_min`, not `f1_median`.** Cells held out of a fit come from the same
donors, acquired in the same tubes on the same day, so they share every source of
between-donor variation the gate will meet in use. A strategy can score an
excellent held-out F1 on *cells* and fail on the next patient. A gate scoring 0.9
on every donor and a gate scoring 1.0 on nine donors and 0.1 on one share a
median and are not the same gate.

`--transfer-max-donors` (8) and `--transfer-max-cells` (20000) bound the cost —
every fold refits, so it is one fit per donor per label. The donor subset is
drawn at random rather than taken as the largest donors, because a worst-donor
statistic computed on the best-represented donors is optimistic in exactly the
direction that matters. Both caps are logged when they bind.

---

## 5. The results

### 5.1 Abundance

**`population_frequencies.csv`** — `pct_of_cd45_pos` per sample per population,
with `u_pct_points`, `pct_lo`, `pct_hi` from the gate uncertainty.

**Report `pct_of_cd45_pos`. Never report `count`.** `count` is an event count. It
depends on acquisition duration and `--max-events-per-file`, so it is not
comparable between samples and is not a cell number. `cells_per_ul` appears only
when the patient table supplies `wbc_per_ul`.

**`group_comparison_stats.csv`** — Wilcoxon rank-sum or Kruskal-Wallis,
aggregated to one value per sample before testing so donors are the replicates,
not cells.

Read the columns in this order:

1. **`p_adj`** (Benjamini-Hochberg, within test family). Not `p_value`.
2. **`cliffs_delta`** — how large the difference is relative to between-donor
   spread.
3. **`difference_over_gate_u`** — how many multiples of the gate's own
   uncertainty the difference amounts to. **Below 1, the groups differ by less
   than the distance the cut itself travels under resampling.**
4. **`difference_over_total_u`** — the same against gate placement and counting
   together. Always the smaller of the two, and the one to act on where they
   disagree.

A result surviving all four is one where the groups differ, the difference is
large relative to donor variation, and it is larger than both the distance the
gate can move and the noise in the count.

### 5.1a Was there enough to count?

Placement and sufficiency are separate guarantees. A cut through a wide empty gap
is well determined however few cells lie beyond it, so `u_pct_points` can be tiny
for a population built on nine events.

| Column | Reading |
|---|---|
| `u_pct_points` | Gate placement only. What it always meant |
| `u_counting_pct_points` | What the frequency carries from its event count (Wilson half-width at 1σ) |
| `u_total_pct_points` | The two in quadrature |
| `lod_pct`, `loq_pct` | Twenty and fifty events as a percentage of that sample's parent gate |
| `detection` | `quantified`, `detected, below LOQ`, or `below LOD` |

**`detection_limits.png`** counts, per population, how many samples clear each
limit. A population mostly below the limit of quantification cannot be rescued by
re-gating — the numerator is small because few cells were acquired, so the fix is
a longer acquisition or a higher `--max-events-per-file`. A population split
across the limit is the dangerous one: its group difference can be produced
entirely by which samples happened to clear it.

Both limits are computed against the parent-gate events **this run** saw, so
subsampling raises them proportionally. That reports the resolution of the
analysis actually performed, not of the files on disk.

`--lod-events` (20) and `--loq-events` (50) change the conventions.

`difference_over_gate_u` is a screen, not a test. The uncertainty is partly shared
across the run — one panel, one transform, one placement rule — so it cancels in
part when a difference is taken and the ratio is conservative. Below 1 is a
reason to open `threshold_uncertainty.csv` before interpreting, not grounds to
discard.

### 5.2 Compositionality

**`compositional_concordance.csv`**

Population percentages within a sample sum to 100 and cannot vary independently.
A genuine granulocyte expansion mechanically depresses every lymphocyte
percentage, and those depressions test significant while absolute lymphocyte
counts per microlitre are unchanged.

Both parameterisations are tested — raw percentages and centred log-ratio — and
each result is classified:

| Class | Reading |
|---|---|
| `robust_to_composition` | Survives both. Supports interpretation as an independent change |
| `raw_only__possible_composition_artefact` | Significant on percentages, not on log-ratios. The configuration above |
| `clr_only__was_masked_by_composition` | The converse: an independent change hidden by the constraint |

The log-ratio does not recover absolute abundance. Proportional expansion of
every population leaves the composition invariant and is undetectable in
frequency data by construction. Distinguishing expansion from relative expansion
needs cells per microlitre.

### 5.3 Marker expression

**`marker_state_stats.csv`**, **`population_marker_mfi.csv`**

Median intensity and percent positive, per sample per population per marker.
Multiplicity is adjusted **within each measure**, not across both, because the
two are strongly correlated summaries of the same events and pooling them would
over-penalise every result.

Differences rather than ratios are reported for transformed intensities: arcsinh
and logicle admit negative values, on which a ratio is undefined in
interpretation and changes sign without a change in biology.

**`subcluster_marker_shifts.csv`** is the cell-level view. Effect sizes, no
p-values, hypothesis-generating. Event counts are set by acquisition duration, so
a test over pooled events derives its degrees of freedom from instrument time.

---

## 6. Confounding

### 6.1 Covariates

**`confounding_diagnostics.csv`**

A variable confounds only when it **both** differs between groups **and**
associates with the outcome. Either alone is inert, and flagging either alone
would raise an alarm on any study with unequal age distributions. Both conditions
are reported, with a `confounder_risk` verdict.

Adjustment is opt-in and deliberately separate. At single-digit *n* per group a
model carrying group with age and sex spends most residual degrees of freedom on
nuisance terms. Where age is strongly associated with group, the parameters are
not separable at any sample size — there are no young controls from which the age
effect could be estimated independently, so the adjusted estimate extrapolates
beyond the observed range and reports an interval that does not say so.

`--rank-ancova` fits it anyway where degrees of freedom permit, labels every row
`EXPLORATORY`, and records `NOT FITTED` with a reason where they do not.

### 6.2 Batch

**`batch_mixing_stats.csv`**, **`batch_group_confounding.csv`** — requires
`--batch-column`.

Two quantities, both required:

**Magnitude.** iLISI against a permutation null. The null is necessary because
the attainable score depends on the number and relative sizes of the batches, so
an absolute iLISI means nothing on its own.

**Separability.** Cramér's *V* between batch and study group. Where patients and
controls were acquired in distinct periods, *V* approaches 1 and batch is not
distinguishable from the comparison of interest.

**Which channel.** Both quantities above describe the shared embedding and name
no marker, so neither corresponds to an action at the bench.

**`marker_batch_drift.csv`** compares each marker's distribution between batches
by Earth Mover's distance. Read `emd_over_mad` — the distance divided by the
marker's own pooled MAD, which is the form comparable across markers. At or above
0.5 the marker is flagged; `worst_pair` names the two batches.

**`threshold_batch_drift.csv`** is the §4.1 threshold test grouped by batch
instead of by study group.

Both are needed. A threshold is one number per sample, so it registers only drift
that moves the cut. A marker can change its spread, grow a tail, or lose the
separation between its modes while the density minimum between them stays exactly
where it was — the distributional test sees that, the threshold test cannot.

A flag in either is a statement about the assay, not the donors — **unless**
Cramér's *V* is high, in which case neither separates a reagent lot from biology.

`--correct-batch` aligns each marker across batches by monotone quantile mapping;
monotonicity means within-batch cell ordering cannot change, so location and
scale move without structure being invented. **Above `--batch-max-cramers-v`
(default 0.6) correction is refused**, because at that level of confounding
removing the batch effect and removing the biological effect are the same
operation. `--force-batch-correction` overrides and the override is recorded in
the manifest.

Correction applies to the shared embedding, the clustering and the gate-cluster
concordance. Not to per-sample frequencies, marker medians or the differential
tests — those come from per-sample thresholds and are batch-local by
construction. The quantity a batch effect distorts is the single embedding
computed across all samples.

## 7. Conformance across runs

**`specification_conformance.csv`** — requires `--baseline`.

Every check above is internal to one run. `threshold_scale_qc.csv` compares each
threshold against the other samples of the same panel, which finds a single
deviant tube. It cannot find a cohort that moved as a whole, because the
leave-one-out peer median moves with it. After a laser service, a reagent lot
change or six months of drift, **every sample can agree with its peers and
disagree with the assay as it was validated**.

`--write-baseline` records where an accepted run placed each threshold, how
variable it was, how often it needed the fallback, and what the populations came
out at. `--baseline` measures a later run against it.

A failure is not a statement that the run is bad. It says the two runs no longer
place their cuts in the same place, so their frequencies are not the same
measurement and should not be pooled until someone has looked.

Two verdicts are withdrawals rather than measurements: a transform differing from
the baseline's puts thresholds on a different scale, and a redefined population is
a different population. Both read `not comparable`.

The baseline holds summaries and the specification text — no event-level or
patient data — so it belongs in version control beside the config it describes.

`--fail-on-drift` turns the verdict into an exit code, raised after every output
has been written.

## 8. Provenance and reporting

**`run_manifest.txt`** — R version, platform, the version of every package loaded
at run time, git commit and working-tree state where the code was a checkout, the
full invocation, and every option in force.

Written before the expensive stages and rewritten at completion, so an
interrupted run leaves `status: failed` rather than an absent or misleading
record.

**`miflowcyt.md`** — the same run restated against the ISAC MIFlowCyt checklist,
which *Cytometry A*, *Nature* and PLOS check at submission.

| Checklist section | Where it comes from |
|---|---|
| 1. Experiment overview | Marked `TO BE COMPLETED` — a run cannot know intent |
| 2. Specimens and reagents | Marked `TO BE COMPLETED` — clone, vendor, lot |
| 3. Instrumentation | FCS keywords: `$CYT`, `$CYTSN`, `$SYS`, `$DATE`, `$OP`, `$TOT`, and a per-panel detector table from `$PnN`/`$PnS`/`$PnV`/`$PnR`/`$PnE` |
| 4. Data analysis | The run itself: compensation state, transform and derived parameters, gate hierarchy, population spec, threshold source tally |

Two conventions worth knowing when reading it. A section that needs a human says
`TO BE COMPLETED` rather than being omitted, because an absent section reads as
one that did not apply. A keyword the file does not carry says `not recorded in
the FCS file` rather than being dropped, because a missing `$PnV` or `$SPILLOVER`
is often the explanation for something further down.

Compensation is reported as applied, absent, or **mixed**. Mixed means some files
carry a `$SPILLOVER` matrix and others do not, which is worth resolving before
the results are used.

`--no-miflowcyt` skips it.

---

## 9. cyRAVEN or cyCONDOR

They answer different questions on the same files.

| | cyCONDOR | cyRAVEN |
|---|---|---|
| Question | What is in the data? | Did the declaration hold? |
| Direction | Cluster, then name the clusters | Declare, then try to falsify |
| Finds a population nobody expected | Yes, that is its purpose | Only through `--cluster` |
| Says a gate is mis-specified | No prior definition to contradict | Five outputs exist for it |
| Per-sample thresholds | No | Yes, every one |
| Uncertainty on a frequency | No | Yes |
| Reads a FlowJo workspace | Yes | No |
| Trajectory / pseudotime | Yes | No, deliberately |

Use both. A cyCONDOR clustering handed to cyRAVEN through `--external-labels`
comes back as an executable gating strategy with held-out-donor performance and,
with `--export-gates`, a Gating-ML document the instrument can run. The join is on
sample and event index rather than row position, because the two tools subsample
independently and a positional join relabels every cell without raising anything.
