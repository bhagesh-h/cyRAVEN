# One large UMAP coloured by population, with the names on the map

WHAT: a single panel, every population in its own colour, its name
printed at its median position, and a key beside the panel giving each
population's share of the embedded cells.

## Usage

``` r
fig_umap_populations(
  cells,
  outfile,
  panel_label = "",
  width = 11,
  height = 8.5,
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  Embedded cell table with `umap_1`, `umap_2`, `population_label`.

- outfile:

  Path to write the figure to.

- panel_label:

  Panel name for the title. Default `""`.

- width, height:

  Size in inches. Defaults `11` and `8.5`.

- dpi:

  Resolution. Default `300`.

- colors:

  Named list of colours. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.

## Details

WHY: umap_overview.png answers several questions at once, sample mixing,
sex, age, cohort, and to fit six panels in one image each is small, the
points are specks and nothing is labelled. The first question anyone
asks of a UMAP is which island is which population, and that needs one
panel large enough to read, with the answer written on it rather than
matched by eye against a key.
