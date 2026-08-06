# Output

Reference for every file a run writes. A default run on eight samples
produces 25 tables and 22 figures.

## 1. Quality control

| File | Content |
|----|----|
| `recon_diagnostics.png` | Gate boundaries per sample |
| `gating_qc.png` | Thresholds superimposed on the densities they derive from |
| `staining_qc.csv` | Per-sample verdict and exclusion reason |
| `thresholds_used.csv` | Threshold, derivation, cofactor and outlier status per sample and marker |
| `threshold_scale_qc.csv` | Panel median and robust *z* per marker |

These constrain the interpretation of every subsequent file. A run whose
gates are misplaced produces internally consistent statistics that are
unusable.

## 2. Uncertainty

| File | Content |
|----|----|
| `threshold_uncertainty.csv` | Sampling and method components per sample and marker, their quadrature sum, the valley’s relative depth, and the fraction of resamples that found it |
| `uncertainty_budget.csv` | Contribution of each threshold to each population’s uncertainty, per sample |
| `frequency_uncertainty.png` | Per-sample frequencies with their uncertainty, one row per population |
| `uncertainty_budget.png` | The same budget as a stacked median contribution |

`bootstrap_valley_rate` is the column to read first. A threshold found
in every resample is a population boundary; one found in half of them is
histogram noise that happened to clear the depth rule, and the two are
indistinguishable in `thresholds_used.csv`.

`u_combined` is `NA` rather than zero where it could not be computed,
which happens for a threshold taken from a separate control tube.
`n_terms_missing` on the frequency table counts how many of a
population’s markers contributed nothing, so a small total that is small
only because most terms are absent is distinguishable from a genuinely
tight one.

Written by default. `--no-uncertainty` skips it and restores the
previous output exactly.

## 3. Conformance

Written by `--baseline`.

| File | Content |
|----|----|
| `specification_conformance.csv` | Per marker: this run’s threshold against the baseline’s, scaled by the baseline’s own spread, with a verdict |
| `specification_conformance_populations.csv` | The same for population frequencies |
| `specification_changes.csv` | Populations added, removed or redefined since the baseline |

Three verdicts. `pass` is within tolerance, `qualify` is a marker worth
looking at before pooling the two runs, `fail` says the cuts are no
longer in the same place and the frequencies are not the same
measurement.

`not comparable` is a fourth value and means the comparison was
withdrawn rather than made: the transform differs from the baseline’s,
so thresholds are on different scales, or the population was redefined,
so it is a different population and a drift statistic would be a
category error.

`--write-baseline` writes the reference itself. It contains summaries
and the specification text, no event-level or patient data, so it can be
version-controlled beside the config it describes.

## 4. Abundance

| File | Content |
|----|----|
| `population_frequencies.csv` | Percentage of parent per population per sample |
| `population_marker_mfi.csv` | Median transformed intensity and percent positive per sample, population and marker |
| `functional_markers.csv` | The same restricted to functional marker blocks declared in the config |
| `population_ratios.csv` | Derived ratios where declared |
| `gate_counts.csv` | Event counts at each level of the gate hierarchy |
| `population_frequencies.png` | Pooled composition across samples |
| `population_marker_heatmap.png` | Marker intensity against declared populations, z-scored |

`count` is an event count and scales with acquisition duration. It is
not an abundance. Report `pct_of_cd45_pos` unless absolute counts were
supplied through `--absolute-counts`.

## 5. Embedding

| File | Content |
|----|----|
| `cells_umap.csv` | Per-cell coordinates with sample, population and panel |
| `umap_overview.png` | Embedding by population and by sample |
| `umap_markers.png` | One panel per marker |
| `umap_density.png` | Cell density in embedding space |
| `umap_density_by_group.png` | The same, faceted by study group |
| `umap_overview_by_group.png` | Combined embedding with one column per group |
| `umap_multigraph_overlay.png` | Per-compartment marker distributions against the reference group |
| `umap.model`, `umap.model.meta.rds` | Persisted model, written by `--save-umap-model` |

A persisted model permits projection of subsequent batches into the same
coordinate space, holding cluster positions fixed between runs.
Projection is refused where the new data lack markers the model used,
since zero-filling would produce coordinates that are geometrically
valid and biologically meaningless.

## 6. Inference

| File | Content |
|----|----|
| `group_comparison_stats.csv` | Abundance between groups, per population |
| `group_comparison.png` | The same as per-sample distributions |
| `marker_state_stats.csv` | Differential state per population and marker |
| `marker_state.png` | The same as per-sample distributions |
| `functional_markers_stats.csv` | Functional blocks between groups |
| `population_ratios_stats.csv` | Declared ratios between groups |
| `compositional_clr_stats.csv` | Abundance tests on centred log-ratios |
| `compositional_concordance.csv` | Classification of each result against both parameterisations |
| `paired_comparison_stats.csv` | Paired designs, written by `--paired-column` |
| `covariate_adjusted_stats.csv` | Rank ANCOVA, written by `--rank-ancova` |
| `subcluster_marker_shifts.csv` | Pooled-event marker shifts per compartment |

Every file above except the last uses one value per sample.
`subcluster_marker_shifts.csv` operates on pooled events, carries effect
sizes without p-values, and answers a different question from the tests
above it.

## 7. Diagnostics

| File | Content |
|----|----|
| `threshold_drift_stats.csv`, `threshold_drift.png` | Whether thresholds separate by group |
| `confounding_diagnostics.csv` | Covariate imbalance, outcome association, and verdict |
| `batch_mixing_stats.csv` | iLISI against a permutation null |
| `batch_group_confounding.csv` | Cramér’s *V* and correction verdict |
| `batch_diagnostic.png` | Observed against null mixing |
| `run_manifest.txt` | Versions, commit, invocation, options |

## 8. Clustering

Written under `--cluster` and `--explain-clusters`.

| File | Content |
|----|----|
| `unsupervised_clusters.csv`, `.png` | Cluster assignment per cell |
| `cluster_gate_agreement_clusters.csv` | Dominant label and purity per cluster |
| `cluster_gate_agreement_populations.csv` | Recovery of each declared population |
| `subcluster_k_selection.csv` | Silhouette-selected k, from `--auto-subcluster-k` |
| `cluster_gate_proposals.csv` | Learned gate geometry with held-out metrics |
| `cluster_gate_polygons.csv` | Polygon vertices in transformed units |
| `cluster_gate_strategy_<k>.png` | One figure per explained cluster |

## 9. External labels

Written by `--external-labels`, and by `--export-gates` for the last
two.

| File | Content |
|----|----|
| `external_label_gates.csv` | Learned strategy per supplied label, with held-out metrics at each depth |
| `external_label_polygons.csv` | Polygon vertices on the analysis scale |
| `gate_transferability.csv` | Precision, recall and F1 on each donor, from a strategy refitted without that donor |
| `gate_transferability_summary.csv` | Minimum, median, maximum and IQR of F1 across donors |
| `external_label_strategy_<label>.png` | One figure per label |
| `*.gatingml.xml` | ISAC Gating-ML 2.0, linear units, levels chained parent to child |
| `*_polygons_linear.csv` | The same vertices as a table, for redrawing by hand |

Read `f1_min` from the summary, not `f1_median`. A gate scoring 0.9 on
every donor and a gate scoring 1.0 on nine donors and 0.1 on one have
the same median and are not the same gate.

The Gating-ML vertices are subdivided along each edge. An edge that is
straight on the analysis scale is a curve in the linear units the FCS
file stores, so a polygon built from the corners alone would describe a
different region from the one that was fitted and validated. Dimensions
are named by marker symbol, which is what cyRAVEN resolves from `$PnS`.

## 10. Auxiliary

| File | Content |
|----|----|
| `patient_metadata_english.csv` | Patient table after column mapping and value translation |
| `absolute_counts*.csv`, `.png` | Written by `--absolute-counts` |
| `flowjo/` | UMAP-annotated FCS, written by `--flowjo-export` |
| `config_derived.yaml` | Derived parameters, written by `--write-config` |
| `sample_map_template.csv` | Written by `--write-sample-map` |

`--write-config` and `--write-sample-map` terminate after writing and
cannot be combined with a full run.

The FlowJo export writes one pooled file, one per sample and one per
group. Open `_ALL_SAMPLES.fcs`, plot UMAP-1 against UMAP-2, and decode
the Population, SampleID and CohortID parameters using
`population_codes.csv`.

## 11. Failure handling

Analyses beyond the core pipeline execute after the primary outputs are
written and are individually wrapped. A failure logs a warning naming
the analysis and leaves all other output intact. An addition capable of
destroying the result it augments is not an addition.
