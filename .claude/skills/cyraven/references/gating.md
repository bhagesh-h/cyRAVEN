# Gating

Three stages: exclude non-cellular and non-viable events, place per-sample marker
thresholds, evaluate a declarative specification against those thresholds.

## 1. The hierarchy

Four gates precede any marker evaluation. Each is derived from the data, not
transferred between samples.

**Scatter.** The lower forward-scatter bound sits at the deepest minimum of the
log₁₀ FSC-A density. Sub-cellular debris forms a distinct low-scatter mode; the
minimum between it and intact cells is the boundary.

**Singlets.** Coincident events cross the interrogation point over a longer
interval than single cells, inflating pulse area relative to height, so their
FSC-H:FSC-A ratio is low. Events are kept above median − *k*·MAD of that ratio,
computed inside the scatter gate; one-sided, since a high ratio is no sign of a
doublet. `--singlet-mad-k`, default 3.

**Viability.** Amine-reactive dyes enter only cells that have lost membrane
integrity, so dead cells are the bright tail and are excluded at the dye's
density minimum. Panels without a viability marker skip this gate and the
omission is logged. `--viability-marker` overrides auto-detection. The live share
of single cells is in `staining_qc.csv` (`pct_live_of_singlets`), and
`--min-live-pct` fails samples below the assay's minimum; there is no default
because no published one exists. A mostly dead sample reads as mostly live,
since the dye's tallest mode is then the dead cells; `pct_live_if_dominant_dead`
gives the live share for that case where the dye has a minimum below its
dominant mode (a second live mode near zero, common on whole blood, produces one
too: check the dye's plot).

**CD45.** Leukocytes are selected on CD45, expressed across haematopoietic
lineages and absent from erythrocytes and platelets. Where CD45 is not in the
panel, all viable events become the parent and a warning is raised, because
percent-of-leukocytes and percent-of-viable are not interchangeable denominators.

Every population frequency is a percentage of the CD45⁺ parent.

## 2. Per-sample thresholds

### Why per sample

Staining index varies between samples through reagent lot, fluorochrome
degradation, time between staining and acquisition, and autofluorescence. A
threshold transferred between samples is mis-specified in both directions:
conservative in dim samples, permissive in bright ones. That bias is systematic
and correlated with acquisition order, which is the configuration most likely to
be confounded with study group.

### How

Within each sample, for each marker, the density of parent-gate cells is
evaluated and the threshold placed at the minimum separating negative from
positive. The resolution order is fixed:

1. **`config`**: an explicit value in the YAML always wins.
2. **`valley`**: the density minimum. It counts only where it is at least 30%
   as deep as the lower of its two peaks and the density falls significantly
   into it and rises significantly out of it, each slope judged against
   counting noise (SiZer-style, Bonferroni-corrected 5%), so the sparse tail of
   a single mode or a flat shoulder does not produce one. A valley falling *below*
   the control-derived value is rejected, because it sits inside the background
   distribution and would call noise positive.
3. **`control_q995`**: the 99.5th percentile of an unstained control tube, where
   one is declared through `is_control` in the sample map. This is the path for
   unimodal markers.
4. **`tailgate`**: no minimum and no control. The cut sits 3.5 robust standard
   deviations above the dominant mode, which it takes to be the negative one.
   Flagged `needs_review`. A marker most of the parent carries, with no
   minimum, is cut above its positive mode and reads near 0%: declare its cut
   under `thresholds:`, give the population a `parent:` in which it separates,
   or declare the side (`thresholds: CD127: {tail: lower}`, as openCyto's
   `side` argument does). The side is not guessed from the distance to zero, because
   autofluorescence and spillover spread lift negative modes there too. The
   spread is the smaller of the mirrored lower half's MAD and the mode's own
   half-maximum width on the cut side, so a second negative hump (monocytes
   and lymphocytes with different background, as on PBMC) does not widen it.
5. **`quantile_fallback`**: a fixed quantile, only where no tail gate can be
   placed, flagged `needs_review`. It calls the same share positive in every
   sample, which is why it is the last resort.

Two further sources exist. **`manual`** is a per-sample override declared in the
config's `sample_overrides:` block, carrying `reason` and `set_by`; it outranks
everything, and is distinct from `config` because `config` applies one number to
every sample. **`fmo_q995`** takes the place of `control_q995` when a
fluorescence-minus-one control is declared for that marker through the sample
map's `fmo_for` column, which is the better reference because an FMO shows the
negative population under the spillover the real samples experience.
**`fmm_q995`** marks a tube that leaves out several reagents at once: the other
missing fluorochromes are not spreading into the channel either, so its negative
is narrower than a true FMO's and the cut it anchors is permissive. A
**`control_q995`** cut from an unstained tube is permissive for the same reason,
on every channel that receives spreading. Under `--fmo-sets-cut` a true FMO also
lowers a valley that sits above it (`fmo_q995_valley_lowered`), so dim positives
between the FMO's edge and the valley count as positive.

A control-anchored cut keeps the bare source name where the parent indicated no
cut at all, and gains `_valley_rejected` where it did and the control put the
boundary higher. Either way `derived_threshold` and `derived_source` record what
would have been used without the control.

### When a control cannot be one

Before a declared control anchors anything, `control_cut()` checks it, because
two failure modes both raise cuts and delete populations.

A far quantile finds the negative's upper edge only when the control holds
nothing else. A tube with about one percent of bright events puts the 99.5th
percentile in that tail: on real minus-multiple data one channel read 1.0 at the
95th percentile and 5.6 at the 99.5th, a cut above the full stain's own 99th
percentile, so nothing would have been positive. The cut is also placed from the
negative mode and the spread of the half below it, and the lower of the two is
taken. A clean control still gives the quantile.

A minus control must be dimmer in the channel it leaves out. One that is
brighter is not that experiment, whatever the sample sheet says. Comparing the
two negative modes in units of the sample's own spread catches it; the channel
is listed in `controls_refused.csv` and the threshold comes from the sample's
own data instead.

**Check `controls_refused.csv` before reading a run that declares controls.** A
minus-multiple tube commonly leaves several of its declared channels no dimmer
than the full stain, because those markers were barely expressed to begin with.
Those channels gained nothing from the control.

### When a marker never resolves

`spreading_receivers.csv` says why. Spreading widens a channel's negative
population without moving it, which fills in the valley a threshold would sit in.
A marker that both falls back in most samples and receives substantial spreading
is reported as a panel design problem, and no gating strategy recovers it.

`fmo_agreement.csv` checks the cuts that did resolve: the distance between the
derived cut and its FMO equivalent, in units of that cut's own uncertainty.
Within about one they agree; beyond about three, a cut above the FMO is
discarding signal and one below is calling spillover positive.

### Reading the result

`thresholds_used.csv`, one row per sample per marker.

`source` summarises *how* the cut was obtained. A `tailgate` or
`quantile_fallback` means the marker did not separate positive from negative in
that sample. Frequencies from a fallback are not invalid, but they carry no
evidence of separation. **A marker that falls back across most samples is not
resolving in that panel, and no gating strategy will recover it**. Change the
panel or the transform.

`source` does not say how well determined the value is. Two cuts both recorded as
`valley` can differ by an order of magnitude in that respect: one in a wide empty
gap, the other on a shallow dip a slightly different histogram would have missed.
`threshold_uncertainty.csv` supplies that quantity. See
`interpretation.md` §3.

`threshold_scale_qc.csv` reports, per panel and marker, the median threshold
across samples and the robust *z* of each deviation. Flagged rows are the first
candidates for review. This check is **within run**: it finds one deviant tube
among sound ones and is blind to a cohort that moved together, because the
leave-one-out peer median moves with it. For that case use `--write-baseline` and
`--baseline`.

## 3. Writing the specification

### Syntax

The `populations:` block of the `--config` YAML. Each population is a map of
marker to direction.

```yaml
populations:
  CD4 T cells:
    CD3: above
    CD4: above
    CD8: below

  Central memory CD4:
    CD3: above
    CD4: above
    CD8: below
    CCR7: above
    CD45RA: below

  Classical monocytes:
    CD3: below
    CD14: above
    CD16: below
```

An event joins a population when it satisfies every declared direction against
*that sample's* thresholds. Every requirement is ANDed with the CD45⁺ parent.

Populations are labelled most-specific-first, so a cell satisfying both a broad
and a deeper definition receives the deeper label.

### Directions

`above`: value above threshold.
`below`: value below threshold.
`intermediate`: between the threshold and an upper bound derived per sample from
a second density minimum inside the positive fraction. Where no second minimum
exists the population is reported UNAVAILABLE rather than merged into the bright
fraction. The canonical case is monocyte subsetting: classical CD14⁺⁺,
intermediate CD14⁺CD16⁺, non-classical CD14^dim.

`any_of`: a reserved key holding a nested map satisfied when **any** member is.

```yaml
  NK cells:
    CD3: below
    any_of:
      CD56: above
      CD16: above
```

`parent:`: a reserved key naming another population. This population's own
markers are cut on the parent's events rather than across all CD45⁺ cells
(openCyto's hierarchy), the parent's conditions are inherited, and the cuts go
to `thresholds_within_parents.csv`. A parent under 200 events in a sample leaves
its children UNAVAILABLE there.

```yaml
  Tregs:
    parent: CD4 T cells
    CD25: above
    CD127: below
```

It is the right tool where the subset is a minority of the parent (Tregs in CD4
T cells). Where the parent is mostly positive for the marker, check
`thresholds_within_parents.csv`: monocyte subsets cut within monocytes scored
worse against manual gating on one whole-blood cohort than the CD45⁺ cuts.

Scatter channels work like any other marker. `SSC-A: above` is how granulocytes
are declared, since dense cytoplasmic granules scatter light strongly and no
lineage antibody is needed. SSC-A is thresholded on a log₁₀ scale by the same
mode-finding machinery.

### How to build one

1. **Enumerate the panel.** `thresholds_used.csv` from any run lists every
   resolved marker. The `$PnS` keyword of any file in the batch carries the same
   information. Names must match `$PnS` exactly. This is the single most common
   cause of an empty frequency table.
2. **Declare lineage-level populations first** and confirm they score plausibly
   before adding subsets. An error in the CD3 gate propagates to every T cell
   subset beneath it.
3. **Pass it with `--config`.** The same file optionally carries threshold
   overrides, colour assignments and metadata translations.
4. **Read the run log.** A population whose markers are not all present is
   reported UNAVAILABLE and scored as absent.

`system.file("config", "config_cohorts.yaml", package = "cyRAVEN")` shows the
file structure. Replace its populations with those of the panel in use; it is a
worked example for one myeloid-lymphoid panel, not a default.

### Markers read inside a population, and ratios between populations

Two further blocks live in the same `--config` file. Both are opt-in.

`functional_blocks:` reports a marker's intensity and percent positive *within*
an already-defined population, rather than using it to define one, and
`functional_markers_stats.csv` tests those values between groups.

```yaml
functional_blocks:
  CD33 on gated myeloid subsets:
    markers: [CD33]
    populations: [Granulocytes, Monocytes]
  activation:
    markers: [CD69, HLA-DR, CD25]
    require: CD3          # populations whose definition has CD3: above
  exhaustion:
    markers: [LAG-3, TIM-3, PD-1]
    exclude: [Classical monocytes]   # everything else
```

Scope resolves in order: explicit `populations:`, then `require:`, then all
populations minus `exclude:`. When authoring one, check that no marker is scoped
to a population its own threshold helped define. CD15 read inside a
CD15-positive gate is 100 percent in every sample: zero variance, undefined
p-value, and a result that measures the definition rather than the biology.
Declaring this block replaces the built-in blocks rather than adding to them.

`ratios:` declares a ratio of two populations, written to
`population_ratios.csv` and tested like any abundance.

```yaml
ratios:
  gran_lymph:
    label: "Granulocyte:lymphocyte ratio"
    numerator: Granulocytes
    denominator: Lymphocytes
```

There is no default, because a ratio naming populations a panel does not contain
would be silently meaningless.

### Do not fit the specification to the data

A specification fitted to the samples it is then tested against cannot be
falsified by them, and falsification is what threshold drift, phenotype
concordance and gate-cluster concordance are for.

Data-driven proposal exists and is kept separate. `--explain-clusters` derives
two-marker gate geometry for any cluster the specification does not describe,
scored on held-out cells, and writes `cluster_gate_proposals.csv`. Promoting a
proposal to a named population is a manual decision, and the next run reads an
edited specification.

## 4. Transform

Fluorescence spans four to five decades. Untransformed, the negative population
compresses against the axis and mode separation cannot be resolved. Every
threshold is placed on a transformed scale, and the choice decides where it
lands.

### arcsinh: the default

The cofactor sets the width of the quasi-linear region near zero and is estimated
per panel by bisection until the background interquartile range reaches a target,
about 1.6 standard deviations of a negative population (a heuristic, not a
published criterion). The conventional value of 5 comes from mass cytometry and
over-expands the background band on spectrally unmixed fluorescence data, so it
is estimated rather than assumed. `derive_cofactor_pooled()` estimates across
samples so one weakly stained file does not set the scale for the batch, and on
scatter-gated single cells so debris does not set it
(`--cofactor-from-all-events` restores the older estimate).

`--per-marker-cofactor` gives each marker its own cofactor from its negative
population (values below zero, mirrored), the per-channel scale the literature
recommends. Opt-in: it changes every cut. Against manual gating on one
whole-blood cohort it ranked samples no better (mean ρ 0.64 for both), with
larger worst-sample errors and one more staining-QC failure, and on a PBMC
cohort it emptied CD4 T cells in 3 of 30 samples and regulatory T cells and NKT
cells in 6 each. `--transform logicle`, per channel too, scored the same as the
shared cofactor on the whole-blood cohort and emptied regulatory T cells in 4
PBMC samples against 1.

`--cofactor` fixes it. `--cofactor-from-first-sample` restores the older
single-file derivation.

### logicle

`--transform logicle`. The automatic logicle rule with linearisation width
*w* = (*m* - log₁₀(*t*/|*r*|))/2, where *r* is the fifth percentile of the
negative population. `--logicle-m` sets the decades, default 4.5.

Both transforms are linear near zero and both represent compensated negative
values; arcsinh is defined for every real number. What logicle adds is a linear
region whose width is fitted per channel to that channel's negatives. Its top of
scale *t* is the channel's declared range, not the brightest event sampled.

**Parameters are pooled across the panel, not fitted per file.** Per-file fitting,
which some tools do, gives each sample an independent scale and invalidates
cross-sample comparison of medians.

### The choice is not cosmetic

On a seven-colour T cell panel, lineage-level populations agreed within 0.1%
between transforms while CCR7 and CD45RA memory subsets differed by tens of
percent.

The mechanism is visible in `thresholds_used.csv`. Under arcsinh the CCR7
threshold resolved to a density minimum in some samples and a quantile fallback
in others, about 2.5 units apart on that scale. Under logicle a minimum was
resolved in every sample.

**For dim markers without clear bimodality, the transform decides whether a
stable threshold exists at all.** Check the `source` column before quoting any
memory-subset frequency.

`--transform none` bypasses transformation for pre-transformed input.

## 5. Compensation

Emission spectra overlap, so signal from one dye registers in detectors assigned
to others. Uncorrected, this produces apparent positivity in channels the cell
does not express.

The spillover matrix is determined at acquisition from single-stain controls and
stored in the FCS keyword block. cyRAVEN applies it where present and reports its
absence. Spectral instruments write already-unmixed data with no matrix, so
`maybe_compensate()` detects rather than assumes, because applying compensation twice is
as damaging as omitting it.

Matrix construction is out of scope. It belongs in acquisition software, where
the single-stain controls can be inspected.

## 6. Gates out

`--export-gates` writes learned gate geometry as ISAC Gating-ML 2.0, on the
analysis scale through its declared `fasinh` or `logicle` transformation, each
dimension named by the files' `$PnN` (a `--spillover` matrix is written into the
file and its compensated dimensions named by their markers), so it can be
executed at the instrument or read by Cytobank, FlowRepository or FlowJo. A
table of vertices in the linear units the FCS file stores sits beside it.

In that table polygon edges are subdivided before inversion. An edge that is
straight on the analysis scale is a curve in linear units, so inverse-transforming
only the corners would emit a different region from the one that was fitted and
validated.

The reverse direction is absent: cyRAVEN does not read gates from a FlowJo `.wsp`
workspace. Writing Gating-ML directly needs no dependency; parsing a proprietary
workspace needs `CytoML` and `flowWorkspace`.
