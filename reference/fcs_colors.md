# The colour palette every figure draws with

`fcs_colors()` returns the active palette; `set_fcs_colors()` replaces
it, merging over the built-in defaults so a partial specification does
not blank out the keys it omits.

## Usage

``` r
fcs_colors()

set_fcs_colors(colors)
```

## Arguments

- colors:

  Named list of colours; defaults to the package palette. See
  `fcs_colors()`.

## Value

`fcs_colors()` returns the active named list of colours.
`set_fcs_colors()` returns the previous value invisibly, so it can be
restored with [`on.exit()`](https://rdrr.io/r/base/on.exit.html).

## Details

Held in a private package environment rather than a global variable, so
that re-theming a run cannot leak into another session and the state can
be restored (`set_fcs_colors(NULL)`).

`NULL` to restore the defaults.

## See also

[`default_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/default_colors.md)
for the full set of keys and what each controls.

## Examples

``` r
old <- set_fcs_colors(list(gate_highlight = "#0072F0"))
fcs_colors()$gate_highlight
#> [1] "#0072F0"
set_fcs_colors(old)
```
