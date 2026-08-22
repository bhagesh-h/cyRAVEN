# Package index

## Run the pipeline

One call does everything. The stage functions below are the same code,
split so a step can be run, inspected or replaced on its own.

- [`run_cyraven()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_cyraven.md)
  : run_cyraven

## Reading and transformation

Getting FCS files into a matrix on a sensible scale. The transform is
selectable, and its parameters are derived once per panel so every
sample stays on one ruler.

- [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)
  : Read one FCS file and resolve marker symbols
- [`fingerprint_panels()`](https://bhagesh-h.github.io/cyRAVEN/reference/fingerprint_panels.md)
  : Group files into panels by their exact marker set
- [`maybe_compensate()`](https://bhagesh-h.github.io/cyRAVEN/reference/maybe_compensate.md)
  : Apply the spillover matrix if one is present WHY: spectral
  instruments (e.g. Sony ID7000) write already-unmixed data and carry no
  \$SPILLOVER keyword. Compensating twice, or erroring because the
  matrix is absent, are both wrong – detect and report.
- [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md)
  : Build the intensity transform the whole run uses
- [`derive_cofactor()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_cofactor.md)
  : Derive the arcsinh cofactor from the data
- [`derive_cofactor_pooled()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_cofactor_pooled.md)
  : Derive the panel cofactor from SEVERAL samples, not from whichever
  file was read first
- [`derive_logicle_pooled()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_logicle_pooled.md)
  : Derive one set of transform parameters for a panel

## Gating

Thresholds derived from each sample’s own density rather than
transferred from fixed coordinates. The reference for a marker that
resolves no minimum can be an unstained tube or a fluorescence-minus-one
control, and a single sample’s cut can be corrected by hand with the
correction recorded.

- [`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md)
  : Find the deepest density valley between two modes

- [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md)
  : Resolve a threshold for one marker, recording where the value came
  from

- [`apply_gate_hierarchy()`](https://bhagesh-h.github.io/cyRAVEN/reference/apply_gate_hierarchy.md)
  : Apply the full gate hierarchy to one file

- [`staining_verdict()`](https://bhagesh-h.github.io/cyRAVEN/reference/staining_verdict.md)
  : Staining QC verdict per file

- [`sample_override()`](https://bhagesh-h.github.io/cyRAVEN/reference/sample_override.md)
  : Look up a per-sample, per-marker threshold override

- [`parse_fmo_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/parse_fmo_map.md)
  :

  Parse the `fmo_for` column into a marker-to-file map

- [`fmo_for_sample()`](https://bhagesh-h.github.io/cyRAVEN/reference/fmo_for_sample.md)
  : Which FMO file, if any, controls a given marker for a given sample

- [`fmo_agreement()`](https://bhagesh-h.github.io/cyRAVEN/reference/fmo_agreement.md)
  : Distance between a derived cut and its FMO-anchored equivalent

## Acquisition and panel quality

Two questions the gating stage cannot answer on its own: whether the
instrument was doing the same thing throughout the acquisition, and
whether a marker that resolves no density minimum is failing because of
the panel.

- [`find_time_column()`](https://bhagesh-h.github.io/cyRAVEN/reference/find_time_column.md)
  : Locate the Time channel in a read FCS file
- [`detect_time_anomalies()`](https://bhagesh-h.github.io/cyRAVEN/reference/detect_time_anomalies.md)
  : Flag acquisition intervals whose rate or signal departs from the
  file
- [`run_acquisition_qc()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_acquisition_qc.md)
  : Run acquisition-time QC across a cohort
- [`frequency_delta_if_cleaned()`](https://bhagesh-h.github.io/cyRAVEN/reference/frequency_delta_if_cleaned.md)
  : How far each population moves if the flagged windows are excluded
- [`fig_acquisition_qc()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_acquisition_qc.md)
  : Event rate across the acquisition, with the flagged intervals marked
- [`spreading_pairs()`](https://bhagesh-h.github.io/cyRAVEN/reference/spreading_pairs.md)
  : Spreading received by each marker from each other marker
- [`run_spreading_report()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_spreading_report.md)
  : Spreading across a cohort, and which markers it is costing

## Calibration

Channel units are comparable within an instrument and never between two.
Beads of assigned fluorophore equivalence convert them to units that
are.

- [`bead_peaks()`](https://bhagesh-h.github.io/cyRAVEN/reference/bead_peaks.md)
  : Locate bead populations in one channel
- [`fit_bead_calibration()`](https://bhagesh-h.github.io/cyRAVEN/reference/fit_bead_calibration.md)
  : Fit a channel-units to assigned-units calibration from a bead
  acquisition
- [`apply_bead_calibration()`](https://bhagesh-h.github.io/cyRAVEN/reference/apply_bead_calibration.md)
  : Apply a fitted calibration to a sample's linear channel values

## Uncertainty

How far a threshold can move, how many cells were behind the frequency
it produced, and what each of those costs. Every cut is resampled and
re-derived over the settings that placed it, and the resulting spread is
propagated to the populations that read it; the counting term is
reported beside it, because a well-placed cut and a sufficient number of
events are different guarantees. None of this changes a reported number;
it puts a second number next to it.

- [`threshold_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/threshold_uncertainty.md)
  : Uncertainty in one per-sample threshold
- [`counting_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/counting_uncertainty.md)
  : Counting uncertainty and detection limits for a population frequency
- [`population_frequency_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/population_frequency_uncertainty.md)
  : Propagate threshold uncertainty into one sample's population
  frequencies
- [`run_gate_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_gate_uncertainty.md)
  : Run the whole uncertainty analysis over a scored cohort
- [`annotate_gate_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/annotate_gate_uncertainty.md)
  : Compare a between-group difference against the gate uncertainty
  behind it
- [`valley_setting_grid()`](https://bhagesh-h.github.io/cyRAVEN/reference/valley_setting_grid.md)
  : The settings sweep that produces the method component

## Populations

Scoring a declarative population specification.

- [`score_populations()`](https://bhagesh-h.github.io/cyRAVEN/reference/score_populations.md)
  : Evaluate the population spec against one file's thresholded markers

- [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md)
  : Built-in population spec, derived from the supplied gating strategy

- [`default_functional_blocks()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_functional_blocks.md)
  : Functional marker blocks, transcribed from the document's four
  scoping rules

- [`expected_positive_markers()`](https://bhagesh-h.github.io/cyRAVEN/reference/expected_positive_markers.md)
  : Markers each population's gate definition requires to be POSITIVE

- [`compute_population_ratios()`](https://bhagesh-h.github.io/cyRAVEN/reference/compute_population_ratios.md)
  :

  Derived abundance ratios (e.g. CD4:CD8), from the config's `ratios:`
  block

## Embedding

One shared UMAP per marker panel, plus model persistence for later
batches.

- [`select_umap_features()`](https://bhagesh-h.github.io/cyRAVEN/reference/select_umap_features.md)
  : Choose the markers that define the embedding space
- [`plan_subsample()`](https://bhagesh-h.github.io/cyRAVEN/reference/plan_subsample.md)
  : Draw a size-balanced subsample of gated cells across samples
- [`density_weights()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_weights.md)
  : Inverse-density sampling weights in marker space
- [`draw_subsample_rare()`](https://bhagesh-h.github.io/cyRAVEN/reference/draw_subsample_rare.md)
  : Draw a subsample that preserves sparse regions of marker space
- [`run_umap()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_umap.md)
  : Build the UMAP input matrix and run the embedding
- [`save_umap_model()`](https://bhagesh-h.github.io/cyRAVEN/reference/save_umap_model.md)
  : Save a trained UMAP model plus the scaling needed to reproduce its
  input
- [`load_umap_model()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_umap_model.md)
  : Load a saved UMAP model
- [`project_umap()`](https://bhagesh-h.github.io/cyRAVEN/reference/project_umap.md)
  : Project new cells into a saved embedding

## Clustering

Unsupervised clustering as a cross-check on the gate specification, and
reference-defined subclustering.

- [`run_unsupervised_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_unsupervised_clusters.md)
  : Unsupervised clustering of the embedded cells
- [`cluster_gate_agreement()`](https://bhagesh-h.github.io/cyRAVEN/reference/cluster_gate_agreement.md)
  : Cross-tabulate unsupervised clusters against gate labels
- [`subcluster_by_reference()`](https://bhagesh-h.github.io/cyRAVEN/reference/subcluster_by_reference.md)
  : Split each cluster into subclusters DEFINED ON THE REFERENCE GROUP
- [`choose_subcluster_k()`](https://bhagesh-h.github.io/cyRAVEN/reference/choose_subcluster_k.md)
  : Choose k per population from the reference group's cells

## Learned gates

The inverse of the gating stage: given a set of labelled cells, learn a
two-marker gating strategy that selects them. This is what turns “this
cluster matches no described population” into a gate someone can draw.
Descriptive throughout – it proposes gates and never changes a scored
one.

- [`learn_convex_gate()`](https://bhagesh-h.github.io/cyRAVEN/reference/learn_convex_gate.md)
  : Learn a convex two-marker gate for a labelled set of cells
- [`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)
  : Learn a hierarchical gating strategy for one labelled population
- [`explain_unmatched_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_unmatched_clusters.md)
  : Propose gating strategies for clusters the spec does not describe
- [`rank_gate_markers()`](https://bhagesh-h.github.io/cyRAVEN/reference/rank_gate_markers.md)
  : Rank markers by how differently they are distributed in targets and
  the rest
- [`fig_gate_strategy()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_gate_strategy.md)
  : Draw a proposed gating strategy, one panel per gate

## Labels from another tool

Take a cell labelling this package did not produce, a cyCONDOR
clustering for instance, learn a gate that selects it, and find out
whether that gate holds on a donor it was not fitted to. Held-out cells
come from the same donors as the training cells and overstate
transferability; refitting with one donor withheld does not.

- [`read_external_labels()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_external_labels.md)
  : Read cell labels produced by another tool
- [`join_external_labels()`](https://bhagesh-h.github.io/cyRAVEN/reference/join_external_labels.md)
  : Attach external labels to the embedding cell table
- [`explain_external_labels()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_external_labels.md)
  : Learn and validate a gate for each externally supplied label
- [`gate_transferability()`](https://bhagesh-h.github.io/cyRAVEN/reference/gate_transferability.md)
  : Does the gate work on a donor it was not fitted to
- [`apply_gate_strategy()`](https://bhagesh-h.github.io/cyRAVEN/reference/apply_gate_strategy.md)
  : Apply a learned gating strategy to new cells

## Gate export

A learned gate that stays in R cannot be sorted on. These write it out
in the units the FCS file stores, as ISAC Gating-ML 2.0 and as a plain
table of vertices.

- [`write_gating_ml()`](https://bhagesh-h.github.io/cyRAVEN/reference/write_gating_ml.md)
  : Write learned gating strategies as a Gating-ML 2.0 document
- [`polygons_linear_table()`](https://bhagesh-h.github.io/cyRAVEN/reference/polygons_linear_table.md)
  : The same polygons as a plain table in linear units
- [`polygon_to_linear()`](https://bhagesh-h.github.io/cyRAVEN/reference/polygon_to_linear.md)
  : Subdivide polygon edges and return them in linear units

## Statistics

Every test below uses per-sample values, so the replicates are donors
rather than cells.

- [`abundance_measure()`](https://bhagesh-h.github.io/cyRAVEN/reference/abundance_measure.md)
  : Pick the best available abundance measure, and say what it means
- [`stats_group_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_group_comparison.md)
  : Test each population's abundance between groups
- [`fig_group_volcano()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_group_volcano.md)
  : Between-group differences for every population in one figure
- [`stats_marker_state()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_marker_state.md)
  : Sample-level differential-state test: population x marker x group
- [`stats_subcluster_shifts()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_subcluster_shifts.md)
  : Per (cluster, subcluster, marker, study): how far is this study from
  reference
- [`clr_frequencies()`](https://bhagesh-h.github.io/cyRAVEN/reference/clr_frequencies.md)
  : Zero-safe centred log-ratio of each sample's population composition
- [`compositional_concordance()`](https://bhagesh-h.github.io/cyRAVEN/reference/compositional_concordance.md)
  : Concordance between the raw-percentage test and the CLR test
- [`stats_confounding()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_confounding.md)
  : Is a covariate confounded with cohort, and does it track the
  outcome?
- [`stats_rank_ancova()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_rank_ancova.md)
  : Rank ANCOVA: cohort effect after adjusting for covariates
  (EXPLORATORY)
- [`stats_paired_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_paired_comparison.md)
  : Paired / repeated-measures abundance test

## Clinical variables

A severity score, a laboratory value or an outcome flag is neither the
study group nor a confounder: a confounder is screened to decide whether
a group difference can be believed, a clinical variable is the question.
The test follows the variable’s type – Spearman for numeric, Wilcoxon
with Cliff’s delta for two levels, Kruskal-Wallis for more – and
multiplicity is corrected within each variable rather than across all of
them. Association only: a 28-day flag is a two-group comparison, and no
time-to-event model is fitted because the sheet carries no follow-up
time. The figures are in reading order: the variables against each
other, then the association across all of them, then the cohort as one
annotated matrix, then per variable the effects with their bootstrap
intervals. Every signed effect carries a percentile bootstrap interval,
because on ten patients the interval usually spans zero and that is the
finding.

- [`stats_clinical_association()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_clinical_association.md)
  : Clinical variables against population frequencies and marker
  intensities
- [`fig_clinical_correlogram()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_correlogram.md)
  : Spearman correlation between the clinical variables themselves
- [`fig_clinical_landscape()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_landscape.md)
  : The cohort as one picture: populations, samples and every clinical
  variable
- [`fig_clinical_forest()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_forest.md)
  : Ordered effect sizes with bootstrap intervals for one clinical
  variable
- [`fig_clinical_trajectory()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_trajectory.md)
  : Per-patient trajectories across timepoints, split by outcome

## Batch structure

Measure the batch effect first, then decide. Correction is opt-in and is
refused when batch and study group overlap far enough that removing one
removes the other.

- [`lisi_score()`](https://bhagesh-h.github.io/cyRAVEN/reference/lisi_score.md)
  : Local Inverse Simpson's Index of a label over a neighbourhood graph
- [`batch_mixing_report()`](https://bhagesh-h.github.io/cyRAVEN/reference/batch_mixing_report.md)
  : Batch-mixing diagnostic with a permutation null and a confounding
  check
- [`marker_batch_drift()`](https://bhagesh-h.github.io/cyRAVEN/reference/marker_batch_drift.md)
  : Per-marker distributional drift between acquisition batches
- [`emd_1d()`](https://bhagesh-h.github.io/cyRAVEN/reference/emd_1d.md)
  : One-dimensional Earth Mover's distance between two samples
- [`correct_batch()`](https://bhagesh-h.github.io/cyRAVEN/reference/correct_batch.md)
  : Correct batch effects in the transformed marker matrix

## Diagnostics

Questions about the data and the gating that the analysis itself cannot
answer.

- [`stats_threshold_drift()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_threshold_drift.md)
  : Test each marker's per-sample threshold for a cohort difference
- [`write_run_manifest()`](https://bhagesh-h.github.io/cyRAVEN/reference/write_run_manifest.md)
  : Write the run manifest
- [`write_miflowcyt()`](https://bhagesh-h.github.io/cyRAVEN/reference/write_miflowcyt.md)
  : Write a MIFlowCyt-structured report of the run
- [`write_run_report()`](https://bhagesh-h.github.io/cyRAVEN/reference/write_run_report.md)
  : Write a single self-contained HTML report of everything a run
  produced

## Conformance across runs

The within-run check in thresholds_used.csv compares each sample against
its peers and catches one bad tube. It cannot catch a cohort that moved
as a whole, because the peer median moves with it. A baseline written
from an accepted run is the fixed reference that can.

- [`write_spec_baseline()`](https://bhagesh-h.github.io/cyRAVEN/reference/write_spec_baseline.md)
  : Write a conformance baseline from an accepted run
- [`read_spec_baseline()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_spec_baseline.md)
  : Read a conformance baseline
- [`specification_conformance()`](https://bhagesh-h.github.io/cyRAVEN/reference/specification_conformance.md)
  : Test this run against a baseline
- [`spec_fingerprint()`](https://bhagesh-h.github.io/cyRAVEN/reference/spec_fingerprint.md)
  : Canonical text form of a population specification

## Figures

- [`fig_absolute_counts_qc()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_absolute_counts_qc.md)
  : QC heatmap for –absolute-counts: one tile per sample x population

- [`fig_acquisition_qc()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_acquisition_qc.md)
  : Event rate across the acquisition, with the flagged intervals marked

- [`fig_batch_diagnostic()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_batch_diagnostic.md)
  : Batch diagnostic figure: the embedding coloured by batch, plus the
  LISI distribution against its permutation null

- [`fig_clinical_correlogram()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_correlogram.md)
  : Spearman correlation between the clinical variables themselves

- [`fig_clinical_detail()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_detail.md)
  : Per-population detail against one clinical variable

- [`fig_clinical_forest()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_forest.md)
  : Ordered effect sizes with bootstrap intervals for one clinical
  variable

- [`fig_clinical_heatmap()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_heatmap.md)
  : Heatmap of clinical association across populations or markers

- [`fig_clinical_landscape()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_landscape.md)
  : The cohort as one picture: populations, samples and every clinical
  variable

- [`fig_clinical_trajectory()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_trajectory.md)
  : Per-patient trajectories across timepoints, split by outcome

- [`fig_cohort_confusion()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_cohort_confusion.md)
  : Cohort composition heatmap ("confusion matrix")

- [`fig_density_by_sample()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_density_by_sample.md)
  : Per-sample (or per-group) density comparison over the shared
  embedding

- [`fig_detection_limits()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_detection_limits.md)
  : How many samples each population is actually measurable in

- [`fig_frequency_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_frequency_uncertainty.md)
  : Population frequencies with the uncertainty of the gate that
  produced them

- [`fig_functional_markers()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_functional_markers.md)
  : Functional-marker positivity figure

- [`fig_gate_strategy()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_gate_strategy.md)
  : Draw a proposed gating strategy, one panel per gate

- [`fig_gating_qc()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_gating_qc.md)
  : Per-gate QC figure: every applied gate shown per sample with its
  threshold

- [`fig_group_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_group_comparison.md)
  : Grouped abundance figure: one panel per population, bar + SD + every
  sample

- [`fig_group_volcano()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_group_volcano.md)
  : Between-group differences for every population in one figure

- [`fig_marker_grid()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_marker_grid.md)
  : Marker-expression grid over the shared embedding

- [`fig_marker_state()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_marker_state.md)
  : Differential-state figure, drawn by the baseline's own comparison
  figure

- [`fig_marker_umaps_by_group()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_marker_umaps_by_group.md)
  : One UMAP per marker, faceted by study group

- [`fig_multigraph_overlay()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_multigraph_overlay.md)
  : FlowJo-style Multigraph Overlay: every cluster x every marker, one
  peak per group

- [`fig_population_frequencies()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_population_frequencies.md)
  : Population frequency figure – one lettered panel per population

- [`fig_population_marker_heatmap()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_population_marker_heatmap.md)
  : Population x marker phenotype heatmap

- [`fig_population_ratios()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_population_ratios.md)
  :

  Population-ratio figure – same layout as fig_group_comparison(), one
  panel per ratio defined in the config's `ratios:` block

- [`fig_populations_by_batch()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_populations_by_batch.md)
  : Population abundance per acquisition batch

- [`fig_recon_diagnostics()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_recon_diagnostics.md)
  : Initial reconnaissance / QC diagnostic figure – run on EVERY batch

- [`fig_threshold_drift()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_threshold_drift.md)
  : Threshold-drift figure: per-sample cut for every marker, coloured by
  cohort

- [`fig_umap_overview()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_umap_overview.md)
  : Population + sample + covariate colouring panel set

- [`fig_umap_overview_by_group()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_umap_overview_by_group.md)
  : Combined + per-group UMAP small multiples

- [`fig_uncertainty_budget()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_uncertainty_budget.md)
  : Where each population's uncertainty comes from

- [`fig_unsupervised_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_unsupervised_clusters.md)
  : Unsupervised-clustering figure: the embedding by cluster, next to
  the same embedding by gate label

- [`theme_cyto()`](https://bhagesh-h.github.io/cyRAVEN/reference/theme_cyto.md)
  : Publication-oriented minimal theme WHY: consistent, legible defaults
  across every emitted figure; outward ticks and frameless legends read
  better in print than ggplot2 defaults.

## Inputs

The sample sheet is one row per FCS file, carrying its identity, its
subject’s attributes, its study variables and any externally measured
counts. It is read into the same structures the three separate loaders
below produce, by calling the same code, so the two routes cannot
diverge.

- [`read_samplesheet()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_samplesheet.md)
  : Read the unified sample sheet
- [`samplesheet_columns()`](https://bhagesh-h.github.io/cyRAVEN/reference/samplesheet_columns.md)
  : Reserved column names of the unified sample sheet
- [`load_patient_table()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_patient_table.md)
  : Read, filter, rename and translate the patient table
- [`load_sample_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_sample_map.md)
  : Read and validate the sample map
- [`load_absolute_counts()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_absolute_counts.md)
  : Load, normalize, parse and sample-match an –absolute-counts input

## Metadata and export

- [`export_flowjo_fcs()`](https://bhagesh-h.github.io/cyRAVEN/reference/export_flowjo_fcs.md)
  : Write UMAP-annotated FCS files for FlowJo

## Appearance and session state

- [`default_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_colors.md)
  : Every colour this script draws with, in one place
- [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md)
  [`set_fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md)
  : The colour palette every figure draws with
- [`reference_group()`](https://bhagesh-h.github.io/cyRAVEN/reference/reference_group.md)
  [`set_reference_group()`](https://bhagesh-h.github.io/cyRAVEN/reference/reference_group.md)
  : The cohort every comparison is drawn against
