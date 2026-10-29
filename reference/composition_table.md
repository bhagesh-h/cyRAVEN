# Share of the embedded cells per label, with its spread across donors

Each embedded cell carries ONE label, its most specific population (or
its cluster), so the shares sum to 100. That is why the numbers differ
from population_frequencies.csv, where the gates overlap – "CD3 pos"
there includes the CD4 and CD8 T cells, and here it holds only the CD3+
cells no deeper population claimed. Under `--subsample rare` the cells
carry a `sampling_weight`, and the shares are weighted back to the true
composition.

## Usage

``` r
composition_table(
  cells,
  label_col = "population_label",
  phenotype = NULL,
  group_col = NULL,
  colours = NULL,
  name_of = NULL,
  levels = NULL,
  marker_pct = NULL
)
```

## Arguments

- cells:

  Embedded cell table with `sample_id` and the label column.

- label_col:

  Label column.

- phenotype:

  Optional named character, label -\> defining markers.

- group_col:

  Optional study-group column in `cells`.

- colours:

  Optional named character, label -\> colour, carried into the table so
  a row can be matched to the figure.

- name_of:

  Optional named character, label -\> a second name (explore mode's
  matched subset).

- levels:

  Optional label order. Every label in it gets a row, at 0% where
  absent, and rows keep this order – so the per-timepoint tables and
  figures list the cell types in the same order as the pooled one.

- marker_pct:

  Optional numeric matrix, labels x markers, of the percent of each
  label positive for each marker; appended as `pct_pos_<marker>`.

## Value

data.frame, largest share first and the catch-all last.
