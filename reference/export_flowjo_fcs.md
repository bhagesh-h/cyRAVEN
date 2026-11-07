# Write UMAP-annotated FCS files for FlowJo

Write UMAP-annotated FCS files for FlowJo

## Usage

``` r
export_flowjo_fcs(
  cells,
  outdir,
  concat = TRUE,
  groups = TRUE,
  group_col = "cohort",
  log = message
)
```

## Arguments

- cells:

  data.frame/data.table of per-cell rows, or a path to cells_umap.csv.
  Requires sample_id, umap_1, umap_2, population_label.

- outdir:

  directory to create and write into.

- concat:

  also write the concatenated \_ALL_SAMPLES.fcs.

- groups:

  also write one concatenated *GROUP*.fcs per study group.

- group_col:

  column holding the study group. Defaults to `cohort`, the name the
  pipeline gives the grouping variable in cells_umap.csv; absent that
  column the group files are skipped with a note rather than failing,
  since an ungrouped run is legitimate.

- log:

  message sink; the main pipeline passes its own log_msg so the export's
  progress lands in the same run log as everything else.

## Value

character vector of written FCS paths, invisibly.
