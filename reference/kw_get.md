# Pull one keyword out of an FCS keyword list, trying several spellings

Pull one keyword out of an FCS keyword list, trying several spellings

## Usage

``` r
kw_get(kw, ...)
```

## Arguments

- kw:

  keyword list from
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)

- ...:

  candidate keyword names, in preference order

## Value

the first non-empty value as a length-1 character, or NA
