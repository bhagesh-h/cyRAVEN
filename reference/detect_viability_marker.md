# Locate the viability dye among the available markers WHY: the dye differs per panel (Zombie-NIR, Zombie-Violet, L/D, 7-AAD, DAPI, PI, Live/Dead...). Detect by pattern so no panel needs code changes, but let the user name it explicitly when the pattern is ambiguous.

Locate the viability dye among the available markers WHY: the dye
differs per panel (Zombie-NIR, Zombie-Violet, L/D, 7-AAD, DAPI, PI,
Live/Dead...). Detect by pattern so no panel needs code changes, but let
the user name it explicitly when the pattern is ambiguous.

## Usage

``` r
detect_viability_marker(markers, explicit = NULL)
```

## Arguments

- markers:

  Character vector of marker names to use.

- explicit:

  The explicit.
