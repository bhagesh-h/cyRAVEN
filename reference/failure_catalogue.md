# Known failure modes, and what each one means for the input

Each entry is a regular expression matched against the error message,
with the meaning and the action to take. Order matters: the first match
wins, so specific patterns precede general ones.

## Usage

``` r
failure_catalogue()
```

## Value

list of list(pattern, meaning, action)
