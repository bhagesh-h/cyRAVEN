# Immune subset definitions resolvable from a given marker panel

Proposes the subsets the supplied markers can actually define, in the
same format as
[`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md):
a named list of populations, each a list of marker -\> direction
requirements.

## Usage

``` r
immune_subset_spec(markers, parent = "CD45")
```

## Arguments

- markers:

  Character vector of marker names present in the panel.

- parent:

  Marker gating the parent population, added to every definition so
  subsets sit inside the same parent as the lineages. `NULL` to omit.

## Value

Named list of population definitions; empty when the panel resolves
none.

## Examples

``` r
immune_subset_spec(c("CD3", "CD4", "CD8", "CD38", "HLA-DR"))
#> $`CD4 T cells`
#> $`CD4 T cells`$CD3
#> [1] "above"
#> 
#> $`CD4 T cells`$CD4
#> [1] "above"
#> 
#> $`CD4 T cells`$CD8
#> [1] "below"
#> 
#> 
#> $`CD8 T cells`
#> $`CD8 T cells`$CD3
#> [1] "above"
#> 
#> $`CD8 T cells`$CD8
#> [1] "above"
#> 
#> $`CD8 T cells`$CD4
#> [1] "below"
#> 
#> 
#> $`T cells`
#> $`T cells`$CD3
#> [1] "above"
#> 
#> 
#> $`CD4 T cells CD38+ HLA-DR+ (activated)`
#> $`CD4 T cells CD38+ HLA-DR+ (activated)`$CD3
#> [1] "above"
#> 
#> $`CD4 T cells CD38+ HLA-DR+ (activated)`$CD4
#> [1] "above"
#> 
#> $`CD4 T cells CD38+ HLA-DR+ (activated)`$CD8
#> [1] "below"
#> 
#> $`CD4 T cells CD38+ HLA-DR+ (activated)`$CD38
#> [1] "above"
#> 
#> $`CD4 T cells CD38+ HLA-DR+ (activated)`$`HLA-DR`
#> [1] "above"
#> 
#> 
#> $`CD8 T cells CD38+ HLA-DR+ (activated)`
#> $`CD8 T cells CD38+ HLA-DR+ (activated)`$CD3
#> [1] "above"
#> 
#> $`CD8 T cells CD38+ HLA-DR+ (activated)`$CD8
#> [1] "above"
#> 
#> $`CD8 T cells CD38+ HLA-DR+ (activated)`$CD4
#> [1] "below"
#> 
#> $`CD8 T cells CD38+ HLA-DR+ (activated)`$CD38
#> [1] "above"
#> 
#> $`CD8 T cells CD38+ HLA-DR+ (activated)`$`HLA-DR`
#> [1] "above"
#> 
#> 
```
