# Save a trained UMAP model plus the scaling needed to reproduce its input

The scaling parameters travel WITH the model. Re-deriving median/MAD
from a new batch would scale it against itself, so the projection would
land in a space subtly different from the one the model was trained in –
a silent error that looks like a batch effect.

## Usage

``` r
save_umap_model(model, path, scale_params, features, meta = list())
```

## Arguments

- model:

  A trained uwot model.

- path:

  File path.

- scale_params:

  Scaling constants that must accompany the model for a projection to
  land in the same space.

- features:

  Character vector of feature names the model was trained on.

- meta:

  Named list of metadata stored alongside the model. Default
  [`list()`](https://rdrr.io/r/base/list.html).
