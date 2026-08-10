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
interval than single cells, inflating pulse area relative to height. Events are
kept within median ± *k*·MAD of the FSC-H:FSC-A ratio computed inside the scatter
gate. `--singlet-mad-k`, default 3.

**Viability.** Amine-reactive dyes enter only cells that have lost membrane
integrity, so dead cells are the bright tail and are excluded at the dye's
density minimum. Panels without a viability marker skip this gate and the
omission is logged. `--viability-marker` overrides auto-detection.

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
2. **`valley`**: the density minimum. A valley falling *below* the
   control-derived value is rejected, because it sits inside the background
   distribution and would call noise positive.
3. **`control_q995`**: the 99.5th percentile of an unstained control tube, where
   one is declared through `is_control` in the sample map. This is the path for
   unimodal markers.
4. **`quantile_fallback`**: a fixed quantile, flagged `needs_review`.

### Reading the result

`thresholds_used.csv`, one row per sample per marker.

`source` is a three-valued summary of *how* the cut was obtained. A
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
per panel by bisection until the background interquartile range reaches a target.
The conventional value of 5 comes from mass cytometry and over-expands the
background band on spectrally unmixed fluorescence data, so it is estimated
rather than assumed. `derive_cofactor_pooled()` estimates across samples so one
weakly stained file does not set the scale for the batch.

`--cofactor` fixes it. `--cofactor-from-first-sample` restores the older
single-file derivation.

### logicle

`--transform logicle`. The automatic logicle rule with linearisation width
*w* = (*m* − log₁₀(*t*/|*r*|))/2, where *r* is the fifth percentile of the
negative population. `--logicle-m` sets the decades, default 4.5.

The quasi-linear region near zero accommodates compensated negative values, which
arcsinh cannot represent.

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

`--export-gates` writes learned gate geometry as ISAC Gating-ML 2.0 plus a table
of vertices, both in the linear units the FCS file stores, so they can be
executed at the instrument or read by Cytobank, FlowRepository or FlowJo.

Polygon edges are subdivided before inversion. An edge that is straight on the
analysis scale is a curve in linear units, so inverse-transforming only the
corners would emit a different region from the one that was fitted and validated.

The reverse direction is absent: cyRAVEN does not read gates from a FlowJo `.wsp`
workspace. Writing Gating-ML directly needs no dependency; parsing a proprietary
workspace needs `CytoML` and `flowWorkspace`.
