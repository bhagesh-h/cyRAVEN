# Do these marker symbols match any ignore pattern

Each pattern is matched as a whole name, case-insensitively. A pattern
containing `*` is treated as a glob, so `[AF color*` catches all three
autofluorescence channels while `CD16` matches CD16 and not CD161.
Whole-name matching is the point: substring matching would silently drop
markers nobody named.

## Usage

``` r
ignore_channel_match(sym, patterns)
```

## Arguments

- sym:

  Character vector of resolved marker symbols.

- patterns:

  Character vector of patterns, or NULL.

## Value

Logical vector the length of `sym`.
