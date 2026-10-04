# Composite panel key for population x marker tables

WHY: fig_group_comparison() panels on `population`. A marker only means
something inside the population it was measured in – CD16 in monocytes
and CD16 in NK cells are different biology – so the panel unit is the
PAIR. Mirrors fx_panel_key() in the baseline, which solved the same
problem for functional blocks.

## Usage

``` r
mfi_panel_key(mfi)
```

## Arguments

- mfi:

  Per-sample x population x marker summary table.
