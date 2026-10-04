# One line saying what a table contains

WHY IT IS A LOOKUP AND NOT DERIVED. A column list is not a description:
a reader who meets `compositional_concordance.csv` needs to be told that
it classifies each result against both parameterisations, and no amount
of inspecting its headers says that. The lines below are the same
one-liners the `outputs` vignette carries, so the report and the
documentation cannot drift into saying different things about the same
file.

## Usage

``` r
report_table_note(file)
```

## Arguments

- file:

  table file name, with or without a directory.

## Value

A single string, or "" when the table is not known.

## Details

A run with several marker panels writes
`group_comparison_stats_<panel>.csv`, so a name that is not found is
tried once more with its last underscore segment removed. An unknown
table gets no line rather than a vague one – a description that could be
about any table is worse than none, and its absence is a visible prompt
to add it here.
