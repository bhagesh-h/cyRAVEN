# Functional marker blocks, transcribed from the document's four scoping rules

Scope is expressed as a RULE, not a hardcoded population list, so that
adding a population to the config automatically scores it under the
right blocks: exclude = populations to omit (everything else is
included) require = only populations whose definition requires this
marker "above"

## Usage

``` r
default_functional_blocks()
```

## Details

Document rules, verbatim in intent: exhaustion (LAG-3/TIM-3/PD-1/CTLA-4)
– "in each cell subset except for monocytes CD14, CD16 and CD14 int CD16
int" homing (CCR5/CCR9/CXCR3) – "in CD3 pos subsets" activation
(CD69/HLA-DR/CD25/CD57/CD38) – "in CD3 pos subsets" monocyte
(HLA-DR/BTN2A2/BTN3A1/2/3) – "in CD14 pos, CD16 pos, and CD14 int CD16
int subsets"
