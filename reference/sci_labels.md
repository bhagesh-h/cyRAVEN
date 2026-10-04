# Axis labels in the 1x10^7 form these figures conventionally use

WHY NOT scales::label_scientific(): it renders "1e+07", which is code
notation. WHY NOT plain numerals: absolute cell concentrations span 10^3
to 10^7 across populations, so full numerals make the tick column wider
than the panel. Returns an expression vector so the exponent renders as
a true superscript; a plain "10^7" string would print the caret
literally.

## Usage

``` r
sci_labels(x)
```

## Arguments

- x:

  A vector of values.
