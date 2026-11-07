# Merge a –config `colors:` block over the defaults WHY modifyList() and not replacing wholesale: a config that overrides only `gate_highlight` should not silently blank out every other colour to NULL.

Merge a –config `colors:` block over the defaults WHY modifyList() and
not replacing wholesale: a config that overrides only `gate_highlight`
should not silently blank out every other colour to NULL.

## Usage

``` r
apply_color_config(cfg_colors)
```

## Arguments

- cfg_colors:

  The cfg colors.
