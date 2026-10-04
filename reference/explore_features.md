# Eligible explore features for a panel

Every fluorescence marker the panel resolves, plus scatter, minus
anything excluded. Unlike
[`select_umap_features()`](https://bhagesh-h.github.io/cyRAVEN/reference/select_umap_features.md)
there is no lineage preference: the point of explore mode is to assume
nothing about which channels matter.

## Usage

``` r
explore_features(
  panel_markers,
  scatter_names = character(0),
  prefer = NULL,
  exclude = NULL
)
```

## Arguments

- panel_markers:

  Character vector of marker names in the panel.

- scatter_names:

  Character vector of scatter channel names.

- prefer:

  Optional explicit list; when given it replaces the default entirely
  and only its intersection with what is available is used.

- exclude:

  Channels to drop.

## Value

Character vector of feature names.
