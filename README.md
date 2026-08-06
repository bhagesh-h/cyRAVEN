# cyRAVEN <img src="man/figures/logo.png" alt="cyRAVEN logo" width="200" align="right"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/pkgdown.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Documentation: <https://bhagesh-h.github.io/cyRAVEN/>**

An R package for counting immune cell types in blood samples and comparing those
counts between groups of people.

### The problem it solves

A flow cytometer measures several proteins on each of a few million individual
cells and writes the result to a file. To turn that into "this patient has 18%
CD4 T cells", somebody has to decide, for every protein, which cells count as
positive and which do not. That decision is called gating, and it is normally
made by a person dragging shapes on a screen.

That step is where most of the error comes from. Two experts gating the same
files can differ by a third in the population size they report, because staining
brightness varies between samples and everyone draws the line somewhere slightly
different.

cyRAVEN does the gating from rules you write down in advance, places each cut
using the sample's own data, and then spends most of its effort trying to show
you where those rules failed.

### What you give it

- A folder of `.fcs` files, one per sample, straight off the cytometer.
- A spreadsheet saying which file is which sample and which study group it
  belongs to.
- A short text file naming the cell types you expect and the proteins that define
  them, for example a CD4 T cell is `CD3: above`, `CD4: above`, `CD8: below`.
- Optionally a table of patient age, sex and other clinical details.

### What it does with it

1. Reads each file and corrects for the optical bleed between colour channels.
2. Removes debris, clumped cells and dead cells using the standard scatter,
   singlet and viability steps.
3. For every protein in every sample, finds the dip between the dim and bright
   populations and puts the cut-off there. This is the step normally done by
   hand.
4. Labels each cell by checking it against your list of cell types.
5. Draws a UMAP, a two-dimensional map where cells with similar protein levels
   sit near each other, using the same number of cells from every sample so no
   one sample dominates the picture.
6. Compares your groups using one number per person rather than one per cell,
   with Wilcoxon and Kruskal-Wallis tests and correction for multiple testing.
7. Runs a series of checks designed to catch its own mistakes.

### What you get back

About 20 tables and 20 figures. The ones most people use:

| File | Contents |
|---|---|
| `population_frequencies.csv` | The percentage each cell type makes up, per sample |
| `population_marker_mfi.csv` | How brightly each protein is expressed in each cell type, per sample |
| `group_comparison_stats.csv` | One test per cell type, with effect size and a corrected p-value |
| `thresholds_used.csv` | Every cut-off it chose, and whether it found a clean dip or had to guess |
| `umap_overview.png` | The map, coloured by cell type and by sample |
| `run_manifest.txt` | Software versions and settings, so the run can be repeated |

Optional settings add a second, independent clustering that can disagree with
your gates, suggested gates for cell populations you did not list, batch
correction, and files you can open in FlowJo.

## What is different about it

Most cytometry software groups the cells first and works out afterwards what each
group was. cyRAVEN goes the other way: you say what you expect to find, and the
software tries to prove you wrong. Four consequences follow.

### Each sample gets its own cut-off

Staining brightness varies between samples. Reagents age, cells sit for different
lengths of time, and one patient's sample simply glows less than another's. A
cut-off copied from one sample to the next will therefore misclassify cells in
both directions.

cyRAVEN looks at each protein in each sample separately, finds the dip between
the dim and bright groups of cells, and cuts there. A sample that stained weakly
gets a lower cut-off and keeps its positive cells.

When no clean dip exists, it says so rather than pretending. The `source` column
of `thresholds_used.csv` reads `valley` when a real dip was found and
`quantile_fallback` when there was none, which tells you which numbers to trust.

### People are counted, not cells

A cytometer keeps running until someone stops it, so one file might hold 200,000
cells and the next 4 million. If you do statistics on individual cells, your
sample size is really a measure of how long the machine ran, and a single
patient whose sample was recorded for longer can produce an apparently
overwhelming result on their own.

cyRAVEN reduces every sample to one number per cell type before any test runs, so
the sample size is the number of people in the study. Where a number genuinely
cannot be reduced that way, it is reported with an effect size and no p-value.

### A second opinion that can disagree

With one extra setting, cyRAVEN also groups the cells by similarity without
looking at your definitions at all, then compares the two answers.

This catches a specific and common failure. Suppose your CD4 T cell gate reports
0.3%, which looks like these patients have almost no CD4 T cells. If the
independent grouping finds a large cluster of CD4-bright cells sitting under the
label "other", the cells are clearly there and your cut-off is in the wrong
place. The frequency table alone could never tell you that, and software that
only groups cells has no prior claim to contradict.

### It refuses to correct a batch effect that is really your result

Samples processed on different days can differ for purely technical reasons, and
there are standard methods to remove that. The catch is that patients are often
run on different days from controls, so "the batch" and "the difference you are
studying" become the same thing. Correcting then quietly deletes the finding and
leaves a clean-looking picture.

cyRAVEN measures how far the two overlap before touching anything. Past a
threshold it stops and tells you why, rather than returning a result that looks
better than the data supports.

## Installing

```r
# flowCore comes from Bioconductor, so add that repository first
install.packages("BiocManager")
BiocManager::install("flowCore")

remotes::install_github("bhagesh-h/cyRAVEN")
```

`FlowSOM` is optional. Without it the clustering step falls back to the
implementation built into the package and nothing else changes.

## A first run

```r
library(cyRAVEN)

run_cyraven(list(
  dir             = "data/",
  outdir          = "results/",
  recursive       = TRUE,
  sample_map      = "data/sample_map.csv",
  config          = "data/my_panel.yaml",
  group_column    = "cohort",
  reference_group = "Healthy controls"
))
```

The same thing from a shell:

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts","cyraven.R",package="cyRAVEN"))')" \
  --dir data/ --outdir results/ --recursive \
  --sample-map data/sample_map.csv \
  --config data/my_panel.yaml \
  --group-column cohort --reference-group "Healthy controls"
```

Anything that only adds a file is on by default. Anything that changes a number
you already had is off by default:

```bash
  --transform logicle       # flow-standard transform instead of arcsinh
  --cluster                 # unsupervised clustering and the gate cross-check
  --explain-clusters        # learn a gating strategy for undescribed clusters
  --batch-column run_date   # measure the batch effect
  --correct-batch           # and correct it, if the measurement allows
  --auto-subcluster-k       # pick the subcluster count from the data
  --save-umap-model         # so the next batch lands on the same axes
```

Read `recon_diagnostics.png` and `gating_qc.png` before you read any number the
run produced. They show you the gates it drew.

## Where to read more

The detail moved out of this file and into the articles, so that this page stays
short enough to finish.

| Article | Contents |
|---|---|
| [Getting started](https://bhagesh-h.github.io/cyRAVEN/articles/cyRAVEN.html) | The nine stages one by one: read and compensate, derive the cofactor, gate, score populations, embed, test, diagnose, explain clusters, write the manifest. Each with the function that does it and runnable examples of `derive_cofactor()`, `density_valley()` and `stats_group_comparison()`. |
| [Gating and populations](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html) | The four-step gate hierarchy and what happens when CD45 or a viability dye is missing. Why thresholds are per sample and how to read the `source` column. Writing a `populations:` block, with the YAML, the meaning of `above`, `below` and `intermediate`, and a four-step recipe. arcsinh against logicle, and the CCR7 example where the choice moves memory subsets by tens of percent. |
| [Checking the result](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.html) | Nine checks in reading order: the two gate figures, staining QC, the phenotype heatmap, threshold drift by group, the four gate against cluster patterns and what each means, learned gates for undescribed clusters, confounding, iLISI batch mixing with the Cramér's V refusal rule, and the run manifest. |
| [What lands in results/](https://bhagesh-h.github.io/cyRAVEN/articles/outputs.html) | Every file the pipeline can write, grouped as read-first, populations and intensity, embedding, statistics, diagnostics, clustering and learned gates. Says what each column means, which flag produces it, and why `count` is an event count and not a cell number. |
| [Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.html) | Why n is donors and not cells, with a worked example. Rank tests instead of limma and when to switch at about 15 samples per group. The compositional problem and what CLR does and does not fix. Confounders diagnosed rather than adjusted. Benjamini-Hochberg within families, and why differences not ratios on arcsinh. |
| [Using cyRAVEN with cyCONDOR](https://bhagesh-h.github.io/cyRAVEN/articles/with-cycondor.html) | A table of which tool answers which question. Three worked questions: is my frequency drop real, what is that unexpected cluster, how do I hold up across sites. Code to drive both from one sample map, the four reasons their UMAPs differ, and why FlowSOM returns 8 clusters where Phenograph returns 24 on the same cells. |
| [What it does not do](https://bhagesh-h.github.io/cyRAVEN/articles/scope.html) | Nine deliberate omissions with the reason and the trigger for reconsidering each: Harmony and CytoNorm, Phenograph, diffcyt, per-marker cofactors, pseudotime, label-transfer classifiers, reading `.wsp` workspaces, cell-level differential expression, and building a compensation matrix. |

## Citing

```r
citation("cyRAVEN")
```

If you use the unsupervised clustering, cite FlowSOM as well. If you use the
per-sample aggregation for differential state, cite diffcyt. Both are named in
the citation output.

## Licence

GPL-3. See [LICENSE](LICENSE).
