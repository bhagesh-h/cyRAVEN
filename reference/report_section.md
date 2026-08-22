# Build one report section, embedding its figures and tables

Build one report section, embedding its figures and tables

## Usage

``` r
report_section(
  outdir,
  id,
  title,
  description,
  figures = character(0),
  tables = character(0),
  body = NULL,
  open = FALSE
)
```

## Arguments

- outdir:

  the results directory.

- id:

  anchor id, unique per section.

- title:

  section heading, stated as what the section reports.

- description:

  what the figures and tables below show and how to read them.

- figures:

  figure filenames, embedded when present.

- tables:

  table filenames, embedded when present.

- body:

  optional raw HTML inserted before the figures.

- open:

  whether the section starts expanded.

## Value

list(html, nav, n_fig, n_tab, bytes, used)
