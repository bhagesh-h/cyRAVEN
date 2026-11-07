# Count the events below a channel's lower limit, by flowAI's rule

flowAI removes "outliers in the lower limit" of the dynamic range
(Monaco et al. 2016): among a channel's negative values, those at or
below their median minus 3.5 \*
[`mad()`](https://rdrr.io/r/stats/mad.html) / 0.6745, which with R's
[`mad()`](https://rdrr.io/r/stats/mad.html), already scaled by 1.4826,
is 5.2 robust SDs. PeacoQC's `RemoveMargins()` removes events at the
channel's minimum instead (Emmaneel et al. 2022).

## Usage

``` r
lower_limit_outliers(ex, cols)
```

## Arguments

- ex:

  Expression matrix, after compensation.

- cols:

  Column indices judged: the marker channels. Scatter is not, as flowAI
  excludes it by default.

## Value

named integer, events at or below the limit per channel with any, with
attribute `events`, the number of distinct events (one event can be
below the limit in several channels, as flowAI's total_UNIQUE counts).

## Details

COUNTED, NOT REMOVED. In compensated data a channel's most negative
values are mostly spreading: cells bright in a channel that spills into
this one, correctly centred and widely spread. Simulated with 10% of
cells bright in a source that widens the receiver five-fold, the rule
removed 2% of events, every one of them bright in the source, and took
that population from 10.0% to 8.1%. Removing them would bias the
frequencies this package reports; the count shows a lower tail the
analyst should look at, an electronic baseline fault for instance,
without deciding for them.
