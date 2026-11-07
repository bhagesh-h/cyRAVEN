# The defining markers of each declared population, as a short string

"CD3+ CD4+ CD127- CD25+" from `list(CD3 = "above", CD4 = "above", ...)`.
Intermediate and bright requirements read `CD14int` and `CD56hi`, and an
`any_of` disjunction is written in brackets: `(CD16+ or CD56+)`.

## Usage

``` r
population_phenotype(spec)
```

## Arguments

- spec:

  Population specification, as
  [`default_population_spec()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_population_spec.md).

## Value

Named character vector, one string per population.
