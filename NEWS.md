# cyRAVEN 0.4.0

Seven additions, completing the feature backlog in `TODO.md`. Every file the
previous version wrote is still written with the same name, the same columns and
the same values, verified by running both versions on the same cohort and
comparing byte for byte: 34 of 35 outputs identical, none lost, and the one that
differs is `miflowcyt.md`, differing only in its own generation timestamp.

## Acquisition-time quality control

* Every number a run reports for a sample is derived from that sample's pooled
  events, which is correct only if the instrument was doing the same thing
  throughout the acquisition. The Time channel was read and then discarded, so
  nothing checked. A partial clog, a bubble or a drift in laser power makes a
  file two instruments over its run, and one threshold then suits neither half.
* New `acquisition_qc.csv`, `acquisition_qc_bins.csv` and `acquisition_qc.png`.
  The Time channel is binned into equal-width intervals and two quantities are
  tracked: the event rate, which a clog lowers and a bubble spikes, and each
  channel's median, which catches a shift that leaves the rate untouched. Both
  are judged by robust z against the file's own bins.
* Bins are equal in TIME rather than equal in count. A bin defined to hold a
  fixed number of events cannot show that the rate changed.
* Nothing is removed. `acquisition_qc_impact.csv` states how far each population
  would move if the flagged intervals were excluded, which is what turns a flag
  into a decision: a file with a visible clog that moves no population by more
  than its own gate uncertainty does not need re-acquiring.
  `--drop-unstable-events` performs the removal and is recorded in the manifest.
* The detector is written directly rather than wrapping PeacoQC, which would add
  a Bioconductor dependency for a binned robust location test. Emmaneel et al.
  2022, Cytometry A 101:325 is the reference for the approach.
* The figure is drawn at the end of the run rather than where it is computed.
  Drawing it early changed the rendering of nine existing figures; the data
  behind them was provably identical either way, but a new diagnostic that
  silently re-renders every published figure is not an addition.

## Fluorescence-minus-one controls

* Control handling was one unstained tube per panel, used as the reference for
  every marker. That is the right control for where autofluorescence ends and
  the wrong one for where a marker's background ends in a panel, because it
  cannot show spillover. An FMO is the same panel with one reagent left out, so
  its distribution in that channel is the negative population under the
  spreading the real samples experience.
* Two optional sample-map columns, `fmo_for` and `control_group`. The first
  names the markers a file controls for; the second confines a control to the
  batch it was acquired in, since a reagent lot changes between batches.
* New threshold source `fmo_q995`, taking the same place in the resolution chain
  as the unstained control and named apart from it because the two are different
  experiments supporting different claims.
* New `fmo_agreement.csv`, which is the reason to have the feature rather than
  just the control. It reports the distance between the derived cut and its
  FMO-anchored equivalent in units of that threshold's own uncertainty. Within
  about one they agree to the precision either can claim, and the derived cut is
  corroborated by an independent experiment. Beyond about three one of them is
  wrong: far above the FMO the cut is discarding signal, far below it is calling
  spillover positive.

## Per-sample threshold overrides

* `--config` could pin a threshold for the whole run, which corrects one tube by
  applying one number to every sample and so reintroduces the fixed-coordinate
  bias the package exists to remove. A `sample_overrides:` block now corrects one
  sample and one marker.
* New source `manual`, distinct from `config`: the first says a named person
  moved this one cut for a stated reason, the second says the assay declares this
  cut everywhere.
* `thresholds_used.csv` gains `override_reason` and `override_by`, and the run
  manifest lists every override. Both are absent on a run that declares none, so
  that run writes the table it always wrote.
* An override matching no sample or marker in the cohort is reported rather than
  ignored. Silence there means the analyst believes a cut was corrected, the run
  says nothing, and the uncorrected number is published.
* `specification_conformance.csv` reports a hand-set marker as `manually set`
  rather than as agreement. A marker made to match the baseline by hand did not
  conform.

## Spillover spreading

* New `spreading_pairs.csv` and `spreading_receivers.csv`. Spreading is the one
  optical effect that removes a density minimum without moving anything else, so
  it is the most common reason a marker that should resolve does not, and the run
  could previously only report `quantile_fallback` with no cause.
* For each ordered channel pair, the spread of the receiver's negative population
  is compared between cells negative and positive for the source. Restricting to
  the receiver's own negatives is what makes it spreading rather than biology:
  co-expression moves the positive cells, not the width of the negatives.
* The actionable join is to the fallback rate. A marker that both fails to
  resolve in most samples and receives substantial spreading is reported as a
  panel design problem, which no gating strategy fixes. On the demonstration
  cohort CD38 and CCR7 are both flagged that way.
* This is a ranking from the samples in hand, not a spillover spreading matrix:
  the published SSM is computed from single-stain controls, which this package is
  not given.

## Run report

* New `report.html`, written by default. A run writes several dozen files and the
  reading order is the product: each check can invalidate the ones after it, and
  a directory listing presents a failed staining QC and a headline p-value as
  equals. The report puts them in order, states what each section answers, and
  links every claim to the file it came from. A run that excluded samples says so
  above every result.
* Written as HTML directly rather than through rmarkdown, which would need
  pandoc: a report the documented execution path cannot produce would be worse
  than none. Figures are referenced rather than embedded, so the report is read
  in the results directory beside them.

## Opt-in changes to reported numbers

Both of these change existing numbers, so both default to the previous
behaviour and are recorded in the manifest.

* `--subsample rare` weights the embedding draw by inverse local density. A
  population at 0.3% contributes about 60 cells out of 20,000 under a uniform
  draw, which is enough to be a smudge and not enough to be a cluster, so the
  unsupervised cross-check cannot recover it. That makes the sampling the binding
  constraint on this package's own falsification claim for rare populations.
  `cells_umap.csv` gains `sampling_weight` so anything computed from the embedded
  cells can be weighted back to the true composition.
* `--calibration-beads` with `--calibration-values` converts channel units to
  MESF or ERF before any threshold is derived. Arbitrary units are comparable
  within an instrument and never between two, which is the ceiling on what
  `--baseline` can establish: it can say a cohort moved and not whether the assay
  or the cytometer did. A channel whose fit does not hold is left in instrument
  units and said to be, because a calibration nobody checked converts an honest
  arbitrary number into a dishonest absolute one.

# cyRAVEN 0.3.0

## Counting uncertainty and detection limits

* Every frequency now carries a second uncertainty beside the gate one: what it
  gets from the number of events behind it. The budget in 0.2.0 perturbed
  thresholds and nothing else, so a population of twelve events and one of twelve
  thousand, both sitting behind the same clean valley, were reported with the
  same uncertainty. Displacing a well-separated cut moves neither of them much.
  One of those numbers is worth acting on.
* New `counting_uncertainty()`. The interval is Wilson rather than the textbook
  `sqrt(p(1-p)/n)`, which goes to zero as p does and so claims perfect knowledge
  of a population nobody observed. The two agree to three figures once a
  population has a few hundred events, so nothing is given up in the common case
  to fix the boundary.
* `population_frequencies.csv` gains `n_parent_events`,
  `u_counting_pct_points`, `u_total_pct_points`, `pct_lo_total`,
  `pct_hi_total`, `lod_pct`, `loq_pct` and `detection`. `u_pct_points` still
  means gate placement alone and still holds the value it held before, and
  `pct_lo`/`pct_hi` are still built from it. Only the new columns are new.
* `group_comparison_stats.csv` gains `total_u_pct_points` and
  `difference_over_total_u`. `difference_over_gate_u` is untouched. The two
  separate for a rare population, where the cut can be well placed and the
  frequency still rest on too few events: the first ratio passes and the second
  does not.
* The limits follow the clinical event-count convention, `--lod-events` (20) and
  `--loq-events` (50), against the parent-gate events THIS run saw. That means
  `--max-events-per-file` raises every limit in proportion, which is the honest
  reading: a population under the limit of a subsample may be well resolved in
  the whole file, and the answer is to raise the cap rather than to believe the
  limit.
* New `detection_limits.png` counts, per population, how many samples clear each
  limit. A population mostly below the limit of quantification cannot be fixed by
  changing the gating, because the numerator is small for want of cells. One
  split across the limit is the dangerous case, since its group difference can be
  driven by which samples happened to clear it.
* Counting is deliberately absent from `uncertainty_budget.csv`, which answers
  which THRESHOLD a population's uncertainty comes from. Counting is not a
  threshold, and adding a term to that table would have changed a published
  figure to say something the file does not mean.
* The two components are not strictly independent, since displacing a cut also
  changes the count. Quadrature treats them as though they were, which is the
  same approximation the GUM makes for the marker terms; it is stated in the
  documentation rather than left to be discovered.

## MIFlowCyt report

* New `miflowcyt.md`, written by default beside the run manifest. MIFlowCyt is
  the ISAC reporting standard and is checked at submission by Cytometry A, Nature
  and PLOS. Its hardest section to write by hand is the one this package already
  knows in full -- the gating specification, the transform and its derived
  parameters, and how every threshold was obtained -- and its instrument section
  was sitting unused in memory, because `read_fcs_resolved()` keeps the whole
  keyword block and the pipeline read four keys out of it.
* Sections 3 and 4, instrument and data analysis, are completed from the keyword
  block and from the run: cytometer, serial, acquisition software, date, operator
  and event counts per file, then per panel a detector table carrying `$PnN`,
  `$PnS`, voltage, range and whether amplification was linear or logarithmic.
* Sections 1 and 2, experiment intent and specimen biology, are marked
  **TO BE COMPLETED** rather than omitted. An absent section reads as one that
  did not apply; a marked one reads as work outstanding, which is what it is. The
  same rule governs individual keywords: one the file does not carry is written
  as "not recorded in the FCS file", because its absence is a fact about the
  acquisition and is often the reason an analysis could not be done.
* Compensation is reported as applied, absent, or mixed. Mixed is a real and
  dangerous state -- some files carrying a `$SPILLOVER` matrix and others not --
  and both of the other two answers would be false for it.
* `--no-miflowcyt` skips it.

## Per-marker batch drift

* `--batch-column` now also writes `marker_batch_drift.csv` and
  `threshold_batch_drift.csv`. iLISI says whether the batches mix in the shared
  embedding, which is the right question and names nothing anyone can act on. A
  flagged marker names a reagent lot, a detector voltage, or a laser that was
  serviced between runs.
* `threshold_batch_drift.csv` cost nothing to add: `stats_threshold_drift()` was
  always generic in its grouping vector, and the pipeline had only ever handed it
  the study group. It is the same test, grouped by batch.
* `marker_batch_drift.csv` is the part that is new. A threshold is one number per
  sample, so the test above only sees drift that moves the CUT. A marker can
  change its spread, grow a tail, or lose the separation between its modes while
  the density minimum between them stays exactly where it was, and no threshold
  statistic can detect that. New `marker_batch_drift()` compares the
  distributions themselves.
* The statistic is the one-dimensional Earth Mover's distance, new `emd_1d()`:
  the L1 distance between two empirical quantile functions. No binning choice, no
  dependency, and it is in the units of the analysis scale, so it is reported as
  a distance rather than as a score. Since those units differ per marker it is
  also divided by the marker's own pooled MAD, which is what makes two markers
  comparable; `emd_over_mad` at or above 0.5 is flagged.
* Its subsampling borrows and returns the RNG stream, for the reason in 0.2.0.

## Fixes

* The RNG regression test was testing `withr`, not cyRAVEN. Its fixture used
  `withr::local_seed()`, and withr 3.0.2 restores the generator's position but
  not its state vector, so the fixture leaked a changed stream into a test whose
  entire purpose is to prove nothing leaks. The fixture now saves and restores by
  hand, the same way the package does. The package's own guard was correct
  throughout and is unchanged; verified directly by calling
  `threshold_uncertainty()` on a pre-computed vector and comparing
  `.Random.seed` and the next five draws either side.

# cyRAVEN 0.2.0

Three additions. All three are additive in the strict sense: every file the
previous version wrote is still written, with the same name, the same columns in
the same order, and the same values. Verified by running both versions on the
same cohort and comparing byte for byte, which is also how the RNG hazard below
was ruled out.

## Gate placement uncertainty

* Every reported frequency now carries a standard uncertainty saying how far it
  moves when the thresholds behind it move. Two components, combined in
  quadrature per the GUM convention: resampling the parent-gate events, and
  sweeping the settings `density_valley()` itself takes. The second is the
  honest form of the claim that automation removed analyst discretion. It did
  not, it fixed it, and this is how much was fixed.
* The parent gate is a term in every population's budget, because displacing the
  CD45 cut moves cells into and out of the denominator all of them are expressed
  against. On the development cohort it is the largest median contributor, which
  is what the operator-variability literature reports for manually gated data as
  well.
* `group_comparison_stats.csv` gains `difference_over_gate_u`. Below 1, the
  groups differ by less than the distance the cut itself moves. A p-value does
  not disclose that.
* New: `threshold_uncertainty.csv`, `uncertainty_budget.csv`,
  `frequency_uncertainty.png`, `uncertainty_budget.png`. On by default;
  `--no-uncertainty` restores the previous output exactly.
* The reported threshold is never touched. Perturbation runs on copies, so the
  numbers are identical with the analysis on or off.
* The RNG is borrowed and returned. `run_cyraven()` seeds once and STEP 6's UMAP
  cell selection draws from that one stream, so a new step consuming draws would
  have changed which cells are embedded and silently redrawn every UMAP in the
  run. Every entry point saves and restores `.Random.seed`; the embedding is
  bit-identical either way, and a test asserts it.
* `bootstrap_valley_rate` records how often a resample finds the cut at all. A
  threshold that comes and goes is not imprecise, it is unresolved, and it looks
  identical to a solid one in `thresholds_used.csv`.

## Conformance against a baseline

* New `--write-baseline` and `--baseline`. The existing peer check in
  `thresholds_used.csv` compares each sample against the others in its own run,
  which catches one bad tube. It cannot catch a cohort that moved as a whole,
  because the leave-one-out peer median moves with it. A baseline written from an
  accepted run is the fixed reference that can.
* `threshold_scale_qc.csv` keeps its within-run meaning. Repointing it at a
  baseline would have kept the filename and the column names while changing what
  the numbers mean, so the baseline comparison goes to
  `specification_conformance.csv` instead.
* A changed transform is reported once as "not comparable" rather than as every
  marker having drifted, and a redefined population is reported as a different
  measurement rather than a moving one.
* `--fail-on-drift` turns the verdict into an exit code for scheduled runs. It
  raises after the run is marked completed, because every output was still
  produced.

## Labels from another tool

* New `--external-labels`, which takes a CSV of cell labels this package did not
  produce -- a cyCONDOR clustering, for instance -- and learns a gating strategy
  for each one. `learn_convex_gate()` already accepted an arbitrary label; what
  was missing was a way to supply one, a way to find out whether the result
  transfers, and a file the instrument can read.
* The join is on sample and event index, never position. Tools subsample
  independently, so row *i* of one table is not row *i* of another and a
  positional join mislabels every cell without erroring.
* New `gate_transferability()` refits the whole strategy with one donor withheld
  and scores it on that donor, once per donor. The existing held-out metric
  reserves CELLS, which come from the same donors acquired in the same tubes, so
  it can look excellent on a gate that fails on the next patient. Report the
  minimum across donors, not the mean.
* New `--export-gates` writes ISAC Gating-ML 2.0 plus a plain table of vertices,
  both in the linear units the FCS file stores. A gate that stays in R cannot be
  sorted on.
* Polygon edges are subdivided before inversion. An edge that is straight on the
  analysis scale is a curve in linear units, so inverse-transforming only the
  corners would emit a different region from the one that was fitted and
  validated.
* No new dependencies. The XML is written directly rather than through a
  `GatingSet`, which keeps the export available in the container this package
  tests in.

# cyRAVEN 0.1.0

First packaged release. Previously a set of command-line R scripts; the analysis
is unchanged, the delivery is not.

## Standalone coverage

* New `--transform logicle` alongside the existing arcsinh path. Logicle is the
  flow-cytometry convention: linear near zero, so compensated negative values
  stay on scale rather than piling against an axis limit. Parameters are pooled
  per PANEL, not fitted per file, so every sample stays on one ruler and
  cross-sample medians remain comparable.
* The transform choice is not cosmetic. On a seven-colour T-cell panel the
  top-level populations agree to 0.1% between methods while the CCR7/CD45RA
  memory subsets move by tens of percent, because under arcsinh the CCR7
  threshold resolves to a density valley in some samples and a quantile fallback
  in others. Logicle finds a valley in every sample. `thresholds_used.csv`
  records which happened.
* New `--correct-batch`, which aligns each marker across batches by monotone
  quantile mapping. Being monotone means it cannot reorder cells within a batch,
  so it corrects location and spread without inventing structure.
* Batch correction REFUSES when Cramer's V between batch and group exceeds
  `--batch-max-cramers-v`. Where the two overlap that far, removing the batch
  and removing the finding are the same operation.
  `--force-batch-correction` overrides and is recorded in the run manifest.

## Learned gating strategies

* New: `--explain-clusters` derives a gating strategy for any cluster the
  population specification does not describe. `cluster_gate_agreement()` could
  already report that a cluster matched nothing; it could not say what the
  cluster *was*, which left the finding as a cluster number and a cell count --
  nothing anyone could act on at the instrument. `explain_cluster()` learns a
  short sequence of two-marker gates that selects those cells and reports how
  well each performs.
* Gates are convex polygons -- intersections of half-planes -- rather than the
  axis-aligned rectangles a conjunction of one-dimensional thresholds can
  express. Diagonal boundaries (CD4/CD8, CD14/CD16, FSC-A/FSC-H) cannot be
  captured by a rectangle without either admitting contaminants or rejecting
  real cells; `derive_singlet_band()` was already a hand-written special case of
  exactly this, and a polygon is its general form.
* No new dependencies. The objective's gradient is derived in closed form and
  the fit is 24 parameters, so `stats::optim`'s L-BFGS-B solves it directly with
  no automatic-differentiation package; `stats::prcomp` initialises and
  `grDevices::chull` tightens, both already imported.
* Every reported metric is computed on cells held out of the fit, with the
  in-sample value printed beside it. Eight free half-planes fitted to a few
  thousand cells can memorise them, so a resubstitution F1 is optimistic by
  construction and the size of that gap is information the reader needs.
* Strictly descriptive. A proposed gate says where cells sit in marker space,
  not that they are a real population; it carries no p-value, and the scored
  populations, frequencies and tested results are identical whether this runs or
  not.

## Packaging

* Reorganised into a standard R package: documented, namespaced, tested, and
  installable with `R CMD INSTALL`. The command-line entry point survives as
  `system.file("scripts", "cyraven.R", package = "cyRAVEN")` and calls the same
  `run_cyraven()` that is now exported.
* Mutable globals (`COLORS`, `REFERENCE_GROUP`) replaced by a private package
  environment with accessors: `fcs_colors()`/`set_fcs_colors()` and
  `reference_group()`/`set_reference_group()`. A namespace is locked after
  loading, so the previous pattern could not work in a package; assigning into
  the global environment instead is something a package must never do.
* All progress reporting moved from `cat()` to `message()`, so it can be
  suppressed, captured, and redirected. Verbosity is set once per session with
  `options(cyRAVEN.verbose = )` rather than threaded through every call.
* The figure raster ceiling is now `options(cyRAVEN.max_raster_px = )` instead of
  a hard-coded constant.

## Performance

All three changes below are exact: they compute the same numbers by a cheaper
route, and are covered by tests that assert bit-identity rather than closeness.

* **Cofactor derivation is roughly two orders of magnitude faster.** The
  bisection evaluated `asinh()` over the entire background vector on each of its
  40 iterations, around 264 million transcendental calls per sample at 300,000
  events. Because `asinh(./m)` is strictly increasing it maps order statistics to
  order statistics, and R's type-7 quantile interpolates between positions that
  depend only on `n` and `p`, never on `m`. The four order statistics the
  interquartile range needs are therefore selected once, and each iteration costs
  four scalar `asinh()` calls.
* **SOM node updates are O(n) per epoch instead of O(nodes x n).** The per-node
  `colSums(X[mapping == j, ])` loop is replaced by a single `rowsum()` pass,
  removing about 100 million redundant comparisons per `--cluster` run.
* **LISI neighbour selection uses a partial sort.** Only the identity of the k
  nearest neighbours affects the Simpson index, so the full `order()` of every
  distance row was unnecessary.

## Memory

* Each sample's raw expression matrix is released as soon as its transformed
  matrix exists, rather than being held for the whole run. Both were previously
  alive simultaneously during population scoring, which is exactly where an
  uncapped run was killed by the OOM reaper. Unstained controls are held back
  until every sample that could reference them has been scored, and
  `--keep-exprs` still retains everything.

## Fixes

* Long plot subtitles and captions are wrapped to the figure width. `ggplot2`
  does not wrap them, so a long one was drawn on a single line and silently
  clipped at the device edge, mid-sentence, with no warning.
* Heatmap tile labels choose black or white per tile by WCAG contrast ratio
  against the fill they sit on. A single fixed grey was unreadable on the dark
  end of a sequential palette, and which tiles were affected depended on the
  data.
* `--cofactor-from-first-sample` no longer collides with `--cofactor` through R's
  partial matching of `$` on lists. `optparse` omits `NULL`-default options from
  the list it returns, so `opt$cofactor` matched the longer name and returned
  `FALSE`, which was then used as a cofactor of zero: every marker became
  `asinh(x/0) = Inf` and every sample failed staining QC with a message that
  pointed nowhere near the cause.
