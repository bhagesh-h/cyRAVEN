# cyRAVEN: Gated Multi-Sample Flow Cytometry Embedding and Differential Analysis

An end-to-end workflow for batches of flow cytometry standard (FCS)
files. Gating thresholds are derived from each sample's own marker
density rather than transferred from fixed coordinates, populations are
scored from a declarative specification held in YAML, and one shared
'UMAP' embedding is built per marker panel. Differential abundance and
differential state are tested on per-sample summaries so that donors,
not cells, are the replicates, following the aggregation strategy of
Weber and others (2019)
[doi:10.1038/s42003-019-0415-5](https://doi.org/10.1038/s42003-019-0415-5)
. Compositional constraints on population percentages are handled with
the centred log-ratio, batch structure is quantified against a
permutation null rather than silently corrected, and an unsupervised
self-organising map clustering provides a cross-check that can
contradict the gate specification.

## Overview

`cyRAVEN` takes a directory of FCS files and produces a gated, embedded,
statistically tested result set. The pipeline runs in nine stages, each
of which is also callable on its own:

1.  **Read** –
    [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)
    resolves `$PnN`/`$PnS` to marker symbols and applies the spillover
    matrix when the file carries one.

2.  **Transform** –
    [`derive_cofactor_pooled()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_cofactor_pooled.md)
    derives the arcsinh cofactor from the data rather than assuming one.

3.  **Gate** –
    [`apply_gate_hierarchy()`](https://bhagesh-h.github.io/cyRAVEN/reference/apply_gate_hierarchy.md)
    derives scatter, singlet, viability and CD45 gates per sample;
    [`density_valley()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_valley.md)
    places each marker threshold at that sample's own bimodal split.

4.  **Score** –
    [`score_populations()`](https://bhagesh-h.github.io/cyRAVEN/reference/score_populations.md)
    evaluates a declarative population specification (see the `config`
    vignette).

5.  **Embed** –
    [`run_umap()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_umap.md)
    builds one shared embedding per marker panel from a size-balanced
    subsample.

6.  **Test** –
    [`stats_group_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_group_comparison.md)
    for abundance,
    [`stats_marker_state()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_marker_state.md)
    for marker state, both on per-sample values.

7.  **Diagnose** –
    [`batch_mixing_report()`](https://bhagesh-h.github.io/cyRAVEN/reference/batch_mixing_report.md),
    [`stats_threshold_drift()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_threshold_drift.md),
    [`stats_confounding()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_confounding.md)
    and
    [`cluster_gate_agreement()`](https://bhagesh-h.github.io/cyRAVEN/reference/cluster_gate_agreement.md).

8.  **Explain** –
    [`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)
    runs the gating stage backwards: given cells the specification does
    not describe, it learns a two-marker gating strategy that selects
    them. Descriptive only; it proposes gates and never alters a scored
    population.

9.  **Report** – figures and CSVs, plus
    [`write_run_manifest()`](https://bhagesh-h.github.io/cyRAVEN/reference/write_run_manifest.md).

[`run_cyraven()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_cyraven.md)
runs all of it.
`system.file("scripts", "cyraven.R", package = "cyRAVEN")` is a
command-line front end to the same function.

## Statistical stance

Replicates are samples, never cells. Anything computed over pooled cells
is labelled descriptive and carries no p-value, because the number of
cells is a property of acquisition rather than of the design. See the
[`vignette("statistics", package = "cyRAVEN")`](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.md).

## Options

- `cyRAVEN.verbose`:

  one of `"none"`, `"inform"` (default) or `"debug"`; controls progress
  messaging. All messaging goes through
  [`message()`](https://rdrr.io/r/base/message.html), so
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) also
  works.

- `cyRAVEN.max_raster_px`:

  hard ceiling on any figure's pixel dimension (default 30000), below
  the graphics device limit.

## See also

Useful links:

- <https://bhagesh-h.github.io/cyRAVEN/>

- <https://github.com/bhagesh-h/cyRAVEN>

- Report bugs at <https://github.com/bhagesh-h/cyRAVEN/issues>

## Author

**Maintainer**: Bhagesh Hunakunti <bhunakun@uni-bonn.de>
([ORCID](https://orcid.org/0000-0002-5957-8005)) \[copyright holder\]

Authors:

- Bhagesh Hunakunti <bhunakun@uni-bonn.de>
  ([ORCID](https://orcid.org/0000-0002-5957-8005)) \[copyright holder\]
