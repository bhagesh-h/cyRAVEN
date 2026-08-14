# cyRAVEN 1.1.0

Parametric tests beside the rank tests, a feasibility check before either, the
channel control that stops a cohort splitting into panels, a full-size UMAP per
marker, parallel reading, and one fix for a failure that was silent and severe.

Every option, output file, column name and value version 1.0.0 produced is
produced identically, with the single deliberate exception in "A sample that
fails staining QC is no longer assumed to be a control" below. That one changes
numbers, because the numbers it changes were wrong.

## A sample that fails staining QC is no longer assumed to be a control

* `staining_verdict()` reads an **absent** `is_control` column as "might be a
  control", because only `is_control = FALSE` positively asserts that a file is
  a biological sample. On a sheet without that column, any sample whose CD45
  gate found no density minimum was promoted to the unstained reference for its
  whole panel.
* Every threshold in that panel then came from the 99.5th percentile of that one
  sample. When it was really a stained sample, every cut landed at the top of a
  real distribution and every population collapsed. Measured on a real cohort:
  T cells at 0.034% of CD45+ against 4.4% once corrected, with no error raised
  and no warning printed. The frequency table looked like a finding.
* A control is now believed only when the sheet declares one. Absence of
  evidence is not evidence of an unstained tube, and a cohort with no controls
  is the ordinary case.
* When samples would have been promoted, the run now names them, says they are
  being treated as failed samples instead, and points at both remedies.
* `--discover-controls` restores the previous behaviour.
* Where a control reference IS used, the log now states that every threshold in
  the panel is a 99.5th percentile of it, and what a wrong reference looks like
  in `thresholds_used.csv`.
* If you have results from an earlier version on a cohort with no declared
  control, check `thresholds_used.csv`. A `source` of `control_q995` on most
  markers means those frequencies came from this failure and should be re-run.

## Knowing what to write in the config

* `--list-channels` prints one row per acquisition parameter: index, `$PnN`,
  `$PnS`, and the symbol the run resolves it to, which is the name a population
  specification has to use. On a spectral panel those differ, because `$PnS`
  reads `CD45 : SparkUV-387 - Area` and the run uses `CD45`. Reads headers only.
* It reads every file rather than the first, and reports whether they all carry
  the same panel, naming the markers that differ.

## A channel that splits the cohort

* `--ignore-channels` drops named channels before anything else, including
  before the panel fingerprint.
* This matters more than it sounds. The fingerprint is the set of resolved
  marker names, so a channel present in some files and absent in others makes
  them separate panels, each with its own cofactor, embedding and thresholds.
  Autofluorescence channels from spectral unmixing do this routinely: they are
  not stains, and how many are written depends on the sample. Twelve comparable
  files can become seven panels of one or two files each.
* Names match whole and case-insensitively; a name containing `*` is a glob, so
  `[AF color*` catches all of them and `CD16` does not catch `CD161`.

## Testing what is testable

* `design_feasibility.csv` is written whenever a group column resolves. One row
  per group with `n_samples`, `n_donors`, `will_be_tested`, `valid_unpaired` and
  a reason.
* It catches two failures that are otherwise invisible. A group below the
  minimum is skipped silently, and a skipped comparison is indistinguishable
  from one that ran and found nothing. Worse, a donor contributing to more than
  one group makes an unpaired test treat repeated measures as independent, and
  that test runs and produces a well-formed p-value answering a question nobody
  asked.
* `--no-group-tests` keeps the grouping and skips the between-group tests.
  Figures stay split and coloured by group and every per-sample quantity is
  still reported. For a cohort worth describing but too small, or too
  confounded, to test. Diagnostics are unaffected.
* `--min-group-n` exposes the minimum group size, which was fixed at 3.

## The t-test and ANOVA, with their assumptions attached

* `parametric_tests.csv` holds Welch's t-test and Student's t-test for two
  groups, Welch's ANOVA and the classical one-way ANOVA for more, with Cohen's
  *d* and eta squared.
* `posthoc_tests.csv` holds Games-Howell, Tukey HSD and Dunn for every pair.
* All of it is computed on arcsine square root transformed percentages, the
  standard variance-stabilising transform for proportion data, and both scales
  are reported because a difference in arcsine units is uninterpretable.
* Every row carries `shapiro_p`, `brown_forsythe_p` and a `recommended` column.
  Normality is tested on within-group residuals rather than pooled values, which
  are bimodal whenever the groups genuinely differ.
* These sit beside the rank tests and never replace them. No existing table
  changes. `--no-parametric` turns them off.

## Figures are written on white, not on transparency

* Every figure now carries an opaque white background. `ggsave()` takes its
  background from the theme, and the theme left `plot.background` unset, so PNGs
  were written RGBA with a transparent ground. **Measured on one run: 11 of 20
  figures.**
* That is fine inside a white report and unreadable anywhere else. Expand a
  figure in a dark viewer, or drop it on a dark slide, and black axis text sits
  on whatever is behind it.
* Set in three places: `safe_ggsave()`, which every one of the 43 figure writes
  passes through; `theme_cyto()`, so the plot object itself carries a white
  ground however it is rendered; and the raw `png()` device the multigraph
  overlay draws on, which does not go through `ggsave()` at all. A caller
  passing its own `bg` still wins.

## Legends sit top-left on the large multi-panel figures

* `functional_markers.png`, `population_frequencies.png`,
  `population_ratios.png`, `group_comparison.png`, `marker_state.png` and
  `umap_multigraph_overlay.png` put their key at the top left instead of the
  right.
* patchwork collects the guide once for a whole composition and centres it
  vertically on the right. On a tall grid that is thousands of pixels from the
  top: `functional_markers.png` runs to 26 panels and six thousand pixels, so
  finding out what a colour means meant scrolling away from the panels and back.
* **The standalone UMAP figures are unchanged.** They are single panels where
  the key already sits beside what it describes, and `threshold_drift.png` keeps
  its right-hand key for the same reason.

## Group colours are separated hues, and no two groups share one

* Study, cohort and timepoint groups take green, red, purple, amber, blue, pink
  in that order, with green anchoring the reference. Each is roughly a sixth of
  the colour wheel from the last.
* **A fourth group used to be drawn in the reference's green.** The palette held
  three colours and was recycled with `rep_len()`, so with four or more groups
  two of them came out identical, with nothing said about it. Past the palette
  the run now generates additional hues instead of repeating. The same bug and
  the same fix applied to the multigraph overlay, where overlaid curves sharing
  a colour are unreadable rather than merely ambiguous.
* Small non-study category sets no longer get red, orange and yellow. Those are
  the first three entries of the population palette, which is ordered by hue, so
  taking the first *n* handed two or three categories three neighbouring
  colours: the hardest triple to separate at bar-chart size and nearly identical
  under red-green colour vision deficiency. Up to six categories the picks are
  now spread across the palette instead.
* `fig_group_comparison()` passes the reference it resolved rather than relying
  on the global default, so a figure's colours no longer depend on what an
  earlier call happened to leave set.

## Long panel keys no longer overrun the panel beside them

* A panel's key is its y-axis title, and a y-axis title is rotated, so the space
  it has is the panel's height rather than its width. Population names fit. The
  keys `functional_markers.png` uses are `block: population, marker`, reaching
  60 characters, and they ran out of the panel onto the neighbour's axis numbers
  and tag letter.
* Long titles now wrap onto two or three lines, with the left margin growing to
  match, so the tag letters stay legible. Short titles are untouched, so
  `population_frequencies.png` and `population_ratios.png` are unchanged.
* Lines are balanced rather than wrapped at a fixed width: wrapping a
  47-character title at 30 leaves one full line and a stub, which reads badly
  rotated.

## A full-size UMAP for every marker

* `marker_umaps_by_group/` holds full-size UMAPs per marker, written by every
  run. `umap_<marker>.png` pools every sample, and
  `umap_<marker>_by_<category>.png` is written beside it for **every category
  present**, not only the one named by `--group-column`.
* Naming a group column decides what is tested. It says nothing about what is
  worth seeing, and a cohort usually carries several categories: timepoint,
  infection focus, phenotype, sex. Before this, a run grouped by timepoint drew
  a single merged panel for infection focus even though the column was in the
  sheet and already panelled on `umap_overview.png`, and an ungrouped run drew
  no split at all.
* Pooled and split are both kept, because they answer different questions.
  Splitting asks whether a marker sits differently between categories; pooling
  asks where the marker is at all, which no split figure can show, because each
  facet holds only a subset of the cells.
* What is never faceted: numeric columns, which would give one panel per
  distinct value; identifiers such as `sample_id` and `patient_id`, which are
  not categories and would multiply the folder by the donor count; and
  single-level columns, since one facet is not a comparison. Past four
  categories the extras are named in the log rather than dropped silently, and
  `--group-column` moves one to the front.
* Neither existing figure covers this. `umap_markers.png` colours by intensity
  but pools the groups and shrinks each marker into a grid cell.
  `umap_overview_by_group.png` splits by group but colours by population or
  covariate, never by intensity. "Is this marker brighter, or somewhere else,
  in one group" needs both at once.
* Within one file the colour scale is shared across the facets and clipped to
  the 1st-99th percentile of the cells drawn, so a colour difference between two
  panels is a real intensity difference. It is not shared between files, for the
  reason `fig_marker_grid()` gives.
* `--no-marker-group-umaps` turns it off. Two PNGs per marker on a grouped run is
  a folder you may not want on a panel of forty.

## Batch correction per cell type, not per file

* `--batch-method cluster` fits one alignment map per marker **per cell type**
  instead of one per marker over the whole file. `cytonorm` is accepted as a
  synonym, because that is the name the method is known by. `quantile`, the
  whole-file map, stays the default and is unchanged.
* Why it matters: a detector shift moves a bright population and a dim one by
  different amounts, so a single map fitted to the pooled distribution
  over-corrects one and under-corrects the other. A per-file map also moves with
  the biology, so it can remove the difference it was meant to preserve.
* It is fitted in two stages, and the second stage is not optional. Clustering
  the raw matrix lets a large batch shift become the dominant variance, so the
  clusters turn out to be the batches, every cluster holds one batch and the
  correction silently does almost nothing. Measured on a synthetic three-batch
  shift: per-cluster alignment fitted that way left a mean between-batch gap of
  1.296 against whole-file alignment's 0.003. Fitting the clustering on a
  whole-file-aligned copy and the per-cluster maps on the original values brings
  it to 0.014, while leaving the true between-cell-type separation at 3.996
  against a true 3.994, where whole-file alignment inflates it to 4.078.
* `--batch-cluster-k` sets the number of cell types, default 10.
* No new dependency. The clustering reuses `run_unsupervised_clusters()`, which
  is already seeded, already stream-safe, and already falls back to a built-in
  SOM when FlowSOM is absent.
* **The Cramer's V refusal is evaluated before the method is consulted.** Both
  methods are refused on identical evidence, and `--force-batch-correction`
  remains the only way past. A better alignment algorithm does not make a
  confounded design correctable.
* `batch_correction.csv` is new. The decision was previously held in memory and
  never written, so no output said whether a correction had run, which method
  fitted it, or why it was refused. Written whenever `--correct-batch` is set,
  including on a refusal.

## `--list-channels` now names the fix, not just the symptom

* When files carry different panels, it reports which channels vary and, where
  the cause is mechanical, prints the exact `--ignore-channels` pattern that
  merges the cohort.
* Two causes are detected: autofluorescence and unmixing artefacts, which are
  not stains and whose count varies per acquisition, and a real marker present
  in a minority of files. Anything else is reported as a genuine panel
  difference and left alone, because that is a decision about the experiment.
* This was previously something you had to know to look for. A fragmented cohort
  raises no error: every file loads and every table is written, while each panel
  quietly gets its own cofactor, thresholds and embedding.

## Speed

* `--read-threads N` reads N files at once. Reading is the slowest stage on a large
  cohort and files are independent, so this scales close to linearly. Peak
  memory scales with it too, since each worker holds one expression matrix.
* Everything after reading stays sequential, because cofactors, thresholds and
  the embedding are derived across samples from a seeded generator that has to
  be consumed in one order for a run to reproduce. Results do not depend on the
  thread count.

## Fixes

* A UMAP marker panel whose marker is constant over the cells drawn no longer
  fails the run. The 1st and 99th percentile came back equal, the colour scale
  was zero width, ggplot collapsed its two identical breaks into one while both
  labels remained, and the run died with "`breaks` and `labels` have different
  lengths" after most of its output had been written.
* `--paired-column` without `--condition-column` was documented as sufficient in
  places. It is not, and never was: the pipeline logs `paired test skipped` and
  writes no paired output.

# cyRAVEN 1.0.0

Explore mode, the statistics catalogue, and the input format. Every option,
output file, column name and value an earlier run produced is produced
identically; nothing is renamed and nothing is removed.

## Explore mode

* `--explore` runs unsupervised discovery over **every eligible channel** --
  scatter and viability included. ignoring the population specification and the
  parent gate. It is what the declared path structurally cannot do: find a
  population nobody declared, and look outside the parent gate.
* Everything it writes lands in `<outdir>/explore/`, including its own
  self-contained `explore_report.html`. No existing output changes. That is
  verified rather than asserted: the demonstration cohort is run with and
  without the flag and every top-level deliverable compared byte-for-byte.
* `--explore-only` needs no specification at all. `--dir` is the whole input.
  For a panel you have not written a config for yet.
* A cluster-level quality gate, the one the unsupervised family uses: cluster
  coarsely, then judge whole clusters by their marker profile, so no per-event
  cutoff is invented. Calls debris, dead and saturated, and records the basis
  for each in `explore_qc_clusters.csv`.
* Cells are equalised per sample after the gate, so a faceted panel's density
  is comparable between groups rather than reflecting group size.
* `explore_suggested_spec.yaml` turns the clusters into a draft specification to
  curate and run back through the supervised path.

## `--maybe-learn`

* The only thing that lets the two analyses see each other, and it is off by
  default. An explore run that used the declared thresholds is no longer a blind
  unsupervised run, and a check on a specification has to be independent of the
  specification to be worth anything.
* With it, the declared run lends explore its per-sample thresholds, so a
  cluster gets a phenotype measured against each sample's own cut --
  `CD19+ HLA-DR+ CD3-`, with a percentage behind every call. instead of one
  read off a colour scale. Its staining QC verdicts and its batch/group
  confounding verdict travel too, the latter into the cluster statistics, so a
  q-value cannot be read without it.
* In the other direction, explore writes `spec_gaps.csv`: declared populations
  spanning several clusters, and clusters no population covers. The first is the
  one worth the flag. a population whose total is flat while a subset inside it
  moves is invisible to any amount of testing on the declared table.

## The statistics catalogue

* `statistical_methods.csv`, on every run, names every commonly reported method
  in the immunophenotyping and cytometry literature. Student's and Welch's *t*,
  one-way and two-way ANOVA with Tukey, repeated measures, Mann-Whitney,
  Kruskal-Wallis with Dunn, chi-squared, Bonferroni, and diffcyt's edgeR,
  limma-voom and GLMM. and says whether this run computed it and why.
* `normality_tests.csv` is the evidence for the choice: Shapiro-Wilk per
  population per group, Brown-Forsythe for equal variances across groups. Its
  `interpretation` column states the thing a p-value alone hides. that
  Shapiro-Wilk on 4 to 10 donors has almost no power, so a non-significant
  result is not evidence of normality, which is itself the argument for the
  rank tests.
* What deliberately does not happen: running every test and reporting all the
  p-values. At these group sizes that is p-hacking with extra steps.

## Input format and the run report

Two changes to how a run is specified and how its result is read, and a version
number that says the interface is now settled. Every option, output file, column
name and value that 0.4.0 produced is produced identically; nothing is renamed
and nothing is removed.

## One CSV and one YAML

* A run needed up to three tables: a sample map keyed by filename, a patient
  table keyed by patient identifier, and an absolute-count export keyed by
  whichever of the two the counting instrument happened to write. Each is a
  different shape with a different key, and the joins between them are where
  cohorts get mislabelled: a patient present in one and absent from another
  produces no error, only an empty covariate panel or a silently dropped count.
* `--samples` takes one sheet with a row per acquired file, carrying every fact
  about that file: what it is, whose it is, which group and batch it belongs to,
  and any externally measured counts in `count.<Population>` columns. One key,
  one file, no joins. `--config` continues to carry the analysis. Between them
  the two files specify a run completely, and the config's new `samples:` section
  can name which column is the group and which the batch so that no further flag
  is needed.
* The split is deliberate. Anything varying per sample belongs in the CSV;
  anything that is one decision for the whole study belongs in the YAML. Putting
  an analysis choice in the CSV would repeat it on every row and invite the rows
  to disagree about it.
* The sheet is a different way to supply the same facts, not a different
  analysis. The reader splits it into the same three structures the pipeline
  always consumed and applies the same coercions by calling the same code:
  `normalise_patient_columns()` was factored out of `load_patient_table()` so the
  two routes cannot diverge in how they parse a date or translate a value. The
  claim is measured rather than asserted: one study run both ways produced 31 of
  31 tables and 22 of 22 figures byte-identical.
* `--sample-map`, `--patient-table` and `--absolute-counts` still work exactly as
  before and are not deprecated. They cannot be combined with `--samples`,
  because two sources of truth for one fact is what the format removes.
* One hazard the three-file shape did not have. A subject attribute is a property
  of a patient, but the sheet has a row per FILE, so a patient with several
  acquisitions repeats it. Two rows can therefore disagree about that patient's
  sex. A reader taking the first value would silently pick one, so
  `read_samplesheet()` reports every disagreement as
  `subject / column: value vs value` and stops. A blank is not a disagreement and
  is filled from the rows that carry a value.
* `--write-samples` writes a sheet covering every file in the input directory,
  with the filename-derived identifiers filled in. Identifiers are never inferred
  from plate order at run time: a wrong guess mislabels a patient and nothing
  downstream reveals it, so a file the sheet does not cover is a fatal error.
* Annotated templates ship as `inst/examples/samples_template.csv` and
  `inst/examples/analysis_template.yaml`.

## Validation before the run rather than during it

* `--check` reads FCS headers and the two input files, reports what the run would
  do, and exits without analysing anything. It names every marker resolved from
  the files, every specification entry matching none of them, whether the sheet
  covers every file, the levels of the group column and their sizes, the number
  of batches, and the study variables available.
* Everything it catches was already knowable before the first event was read, and
  was already being reported: partway through a run that costs many minutes on a
  real cohort. A marker name not matching `$PnS` exactly is the most common cause
  of an empty frequency table, and it was surfacing after the reading stage
  rather than before it.
* It reads keyword blocks rather than event matrices, so it costs seconds
  regardless of cohort size.

## The report is one self-contained file

* `report.html` embedded nothing. It referenced the figures beside it, so moving
  it produced a page of broken images, and a result that cannot survive being
  moved is not a record.
* Every figure is now embedded at full resolution and every table in full. The
  file references nothing, loads no font or script from a network, and works from
  a `file://` path. It can be attached to an email, put in a supplement or
  archived on its own.
* Base64 encoding is written out in `R/base64.R` rather than taken from a
  package: the CRAN options would each add a dependency to an import list that is
  deliberately short and a container image that is pinned, for fifteen lines with
  no edge cases beyond padding. It is vectorised over the raw vector, because a
  per-byte loop over a 6 MB PNG is minutes in R.
* Sections are collapsible, with a sidebar indexing every section, figure and
  table. Figures share one display box whatever their native aspect ratio, zoom
  on click with keyboard control, and download at full resolution from the same
  bytes that are already embedded, so nothing is stored twice. Tables are
  searchable, sortable by column, shown 10, 50, 100 or all rows, and export to
  CSV exactly as filtered and sorted.
* Tables are carried as JSON and rendered in the browser rather than written as
  markup. Over pre-rendered rows the export would re-parse the DOM and the search
  would hide rather than filter; carrying the data means all three operations
  read one array, so what is exported is what was filtered.
* Section headings state what the section reports rather than posing a question,
  and each carries a description of what its figures and tables show and how to
  read them.
* No output can be silently omitted. The named sections put files in reading
  order; anything a section does not name is collected under "Further outputs"
  rather than dropped. A table above 8 MB is named with its row count instead of
  being embedded, which in practice means only the per-cell exports, whose row
  count is the number of cells; the limit is
  `options(cyRAVEN.report_table_max_mb = )`.
* The file is large, tens of megabytes on a full run, and its size is logged.
  That is the cost of self-containment and is stated rather than discovered.

## A failed run explains itself

* A run that failed left the outputs written before it stopped, a manifest marked
  `failed`, and an R error on stderr that whoever opens the results directory
  afterwards may never have seen. Reconstructing what happened meant reading a
  stack trace out of a container log.
* A failed run now writes `report.html` too, with the diagnosis first: what
  stopped it, which stage it reached, the log leading up to that point, what the
  message means in the vocabulary of the analysis rather than of R, and the
  specific next action. Everything produced before the failure is embedded below
  it, because the partial output is usually where the evidence is: a
  `gating_qc.png` written before the failure often shows the cause directly.
* Twelve failure modes carry an interpretation, including a sheet that does not
  cover every file, subject rows that disagree, a path that does not exist inside
  the container, memory exhaustion, an absent optional package, and a marker name
  that resolves to nothing. "Subscript out of bounds" is not a useful thing to
  tell a cytometrist; "the specification names a marker this panel does not
  contain" is the same fact in terms they can act on.
* The report is written by the entry point rather than by an exit handler,
  because an exit handler cannot see the condition that ended the call and
  "something failed" is not a diagnosis. The original error is re-raised
  immediately afterwards, so the exit status and the message a caller sees are
  unchanged, and a failure inside the reporting cannot replace the real error
  with one about reporting.
* The run log is now kept in a bounded in-memory buffer as well as being written
  to stderr, which is what lets the report show what the run was doing when it
  stopped.

## Documentation

* New Inputs article covering both routes, every reserved column, the count
  columns, the config sections, the templates, and the errors the format can
  raise with what each means.
* The Get started article is organised by use case, each with the full command
  sequence: a first run on a new cohort, the full analysis, a cohort that must
  fit in memory, comparing against an accepted baseline, and driving it from R.
* The diagnostics article opens with `--check`, since every other check in it
  applies to a run that has already happened.
* README, the options reference, the output article and the Claude skill all
  document the new flags, the report's behaviour and the failure report.
# cyRAVEN 0.4.0

Seven additions, completing the feature backlog in `dev/TODO.md`. Every file the
previous version wrote is still written with the same name, the same columns and
the same values, verified by running both versions on the same cohort and
comparing byte for byte: 34 of 35 outputs identical, none lost, and the one that
differs is `miflowcyt.md`, differing only in its own generation timestamp.

## Height-only acquisitions

* The reading stage kept area channels only, because a marker recorded as both
  area and height would otherwise be counted twice in the embedding. Applied to a
  file that records height and no area, which older instruments and several
  clinical archives do, that rule resolved zero markers. The run then failed
  several stages later with every population reported UNAVAILABLE, which points
  at the specification rather than at the panel.
* `read_fcs_resolved()` now falls back to pulse height when a file contains no
  area channel, and says so. The rule exists to avoid double-weighting a marker
  that appears twice; where nothing appears twice there is nothing to avoid. The
  substitution is announced rather than silent because height understates a
  bright wide event, so thresholds derived from such a file are not
  interchangeable with those from an area acquisition of the same panel.
* The scatter gate is defined on area and stops without it, so a height channel
  now stands in there too. That alone would have made `derive_singlet_band()`
  divide a channel by itself: a ratio identically one, a MAD of zero, a band of
  zero width, and every event discarded. The singlet gate now requires FSC-H and
  FSC-A to be distinct columns and skips with a log line when they are not.

## Demonstration cohort

* The worked example now uses the graft-versus-host disease dataset distributed
  with flowCore (Brinkman et al. 2007, Biol Blood Marrow Transplant 13:691),
  replacing two CytoTrol control acquisitions partitioned into pseudo-samples
  with randomised group labels.
* The previous cohort could demonstrate calibration and nothing else: its groups
  were randomised, so the correct answer was always that nothing differs, and its
  panel carried neither CD45 nor a viability dye. The replacement has five
  patients, seven visits as genuine acquisition batches, a real clinical contrast
  between GvHD grades, CD45 in the panel, a Time channel, and 2,205 to 66,105
  events per file.
* It requires no download. The data ship inside a package cyRAVEN already depends
  on, so the example is reproducible offline and cannot break when a repository
  moves or a certificate expires.
* Both groups are transplant recipients, so the contrast is disease severity
  rather than disease against health, and the grades are unbalanced across
  patients. `inst/scripts/demo_data.R` states both, and the confounding
  diagnostic is the check to read before interpreting the contrast.

## Documentation

* The demonstration configuration now declares two `functional_blocks:` and one
  entry under `ratios:` in addition to its five populations, so a default run
  writes all 22 figures the package can produce from a cohort with no volumetric
  counting, rather than the 20 that the `populations:` block alone reaches.
  `absolute_counts.png` and `absolute_counts_qc.png` are the two it cannot, and
  the worked example says why and documents the measurement they need instead of
  filling the gap with invented counts.
* `functional_blocks:` and `ratios:` had no syntax reference. Both are now
  documented in the Gating article and in the skill, including the requirement
  that decides whether a block is meaningful: a marker must not be read inside a
  gate its own threshold helped draw, because a percent positive pinned at 100 in
  every sample has zero variance and measures the definition rather than the
  biology.
* The `--absolute-counts` input format is documented in the Output article: the
  sheet layout, how the unit label is read, and the two-key sample matching with
  its ambiguity rule.
* The worked example carries every figure the run produces, each with what it
  measures and what this cohort shows, and each with alt text.
* New `options.Rmd`, a reference for all 83 command-line options with their
  defaults and consequences. An audit found 26 of them documented nowhere in the
  README, the articles or the skill. It also states the convention that decides
  which options are on by default: analyses that only add output are, and the
  five that change numbers a previous run reported are not.
* Maintainer working files moved to `dev/`, which is excluded from the build, the
  container context and version control. `dev/README.md` records which documents
  are authoritative where they disagree.

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
