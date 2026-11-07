# Run an expression with the global RNG stream left exactly as it was found

WHY THIS EXISTS.
[`run_cyraven()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_cyraven.md)
seeds once and every later draw comes from that one stream: the
embedding's cell selection, the clustering, the bootstraps. A step that
spends draws therefore does not merely get its own random numbers, it
shifts every draw after it, and the visible symptom is that adding an
unrelated flag silently moves the UMAP. `geom_jitter()` draws, so a
figure is enough to do it – on a two-panel explore run, drawing panel
1's figure would re-roll panel 2's embedding.

## Usage

``` r
without_spending_draws(expr)
```

## Arguments

- expr:

  Expression to evaluate.

## Details

The package's contract is that a run without `--total-counts` produces
byte-identical output to one with it, apart from the files the flag
adds. Restoring the stream is what makes that true. Same idiom as
[`clin_boot_ci()`](https://bhagesh-h.github.io/cyRAVEN/reference/clin_boot_ci.md)
and `run_flowsom()`.
