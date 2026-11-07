# The marker scaling the embedding and the clustering use

WHAT: `--umap-scale` when it names a method. Under `auto`, the default:
no scaling where the transform already puts every marker on one scale
(logicle, or arcsinh with per-marker cofactors), and the percentile
range otherwise.

## Usage

``` r
umap_scale_method(opt)
```

## Arguments

- opt:

  The run's option list.

## Value

"none", "range", "robust" or "zscore".

## Details

WHY NO SCALING ON A SHARED SCALE. Logicle puts every channel on the same
4.5-decade display scale by construction, and per-marker cofactors give
every marker's negative population the same width, so the markers are
already comparable. Dividing each by its spread (MAD) then only
reweights them by how narrow their negative population is: a mostly
negative marker such as Vd1 has a small MAD and gets its noise
multiplied. Measured on 40,000 cells of one cohort, embedding unscaled
logicle values kept 97.0% of each cell's 15 nearest neighbours in its
own cluster; median/MAD scaling kept 94.8%, and merged islands that
unscaled data keeps apart. FlowSOM itself defaults to no scaling.

WHY THE PERCENTILE RANGE, NOT MEDIAN/MAD, UNDER ONE PANEL COFACTOR. One
cofactor does not put markers on a common scale, so some scaling is
needed; but the MAD of a mostly negative marker is the width of its
negative population, so median/MAD scaling inflates exactly the
background and spreading noise that distort embeddings (Mazza et al.
2018), and does so more the rarer the positive cells are. Dividing by
the 0.5th-99.5th percentile range instead gives every marker the same
span from its negatives to its brightest cells, whatever share is
positive. `robust` remains available by name.
