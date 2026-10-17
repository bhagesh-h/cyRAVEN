# Read the multiplier out of a count column's own header text

WHY THIS IS NEITHER OPTIONAL NOR GUESSED SILENTLY. "PBMC count x 10^6"
and "PBMC count" differ by a factor of a million. Getting it wrong
produces a table that is internally consistent, passes every downstream
check, and is wrong by six orders of magnitude. So the multiplier is
read from the header wherever one is stated, and where none is stated
the values are taken at face value with a NOTE naming the header it
looked at, which the user can confirm against a file they can open.

## Usage

``` r
total_counts_scale(txt)
```

## Arguments

- txt:

  Header text of the count column.

## Value

list(scale, label)
