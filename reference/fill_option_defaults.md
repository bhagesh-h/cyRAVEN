# Fill an option list from the command-line defaults

WHY THIS EXISTS.
[`run_cyraven()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_cyraven.md)
is documented – in the README, in the vignette, in its own examples – as
taking a PARTIAL list: name the two or three options that matter and let
the rest default. That only worked when the list arrived from
`optparse`, which supplies every default itself. Called directly with a
short list, every unnamed option was NULL, and the failures that
produced were silent and remote from their cause: `singlet_mad_k`
missing made the singlet band `median +/- NULL * mad`, which selects no
events, so every marker threshold resolved over an empty vector to NA,
and the run reported that every population's markers were "not in panel"
– naming the panel for a fault in the caller's option list.

## Usage

``` r
fill_option_defaults(opt)
```

## Arguments

- opt:

  Named list, possibly partial.

## Value

`opt` with every unset option filled from its command-line default.

## Details

Taking the defaults from `build_option_list()` rather than restating
them means there is one place a default is written down, and the
programmatic and command-line paths cannot drift apart.
