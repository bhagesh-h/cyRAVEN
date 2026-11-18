# Resolve `parent:` declarations in a population specification

WHY. A cut on the whole CD45+ parent sits wherever the lineages that
dominate that mixture put a valley, which is not where a subset is
defined. Regulatory T cells are CD25+ CD127lo WITHIN CD4+ T cells
(Santegoets et al. 2015), and the CD25 valley across all leukocytes is
not the one that separates them from activated conventional CD4 T cells
(Seddiki et al. 2006); monocyte subsets are split by CD14 and CD16
within monocytes (Abeles et al. 2012). Gating templates express this as
a hierarchy, each gate derived on its parent's events (openCyto; Finak
et al. 2014).

## Usage

``` r
resolve_population_parents(spec)
```

## Arguments

- spec:

  population specification as read from YAML.

## Value

the specification with every `parent:` resolved, carrying attributes
`cut_scope` and `parents`; unchanged when no population names a parent.

## Details

WHAT. A population may name another in its `parent:` key:

    Tregs:
      parent: CD4 T cells
      CD25: above
      CD127: below

Its own markers' cuts are then derived on the parent's events (see
[`derive_scoped_cuts()`](https://bhagesh-h.github.io/cyRAVEN/reference/derive_scoped_cuts.md)),
and it is scored within the parent. The parent's conditions are merged
into the child here, so the definition every other part of the package
reads is the complete one; attribute `cut_scope` records, per population
and marker, which parent's events that marker's cut comes from (NA for
the CD45+ cut), and `parents` the declared parent of each population.

A child may not re-declare a marker its parent requires in a different
direction, and parent and child cannot both carry an `any_of` group:
both would leave the child's relation to its parent undefined.
