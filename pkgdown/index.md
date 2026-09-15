<div class="cy-hero">
<div class="cy-hero-text">
<h1>cyRAVEN</h1>
<p class="cy-hero-lede">
A directory of FCS files becomes population frequencies, marker expression, a
shared UMAP and between-group statistics. Every threshold is placed inside the
sample it applies to, and every frequency carries the uncertainty of the cut and
of the event count behind it.
</p>
<p class="cy-hero-actions">
<a class="cy-btn cy-btn-primary" href="articles/cyRAVEN.html">Get started</a>
<a class="cy-btn" href="articles/gallery.html">See every figure</a>
<a class="cy-btn" href="articles/about.html">Why it works this way</a>
</p>
<p class="cy-hero-meta">
R package &middot; GPL-3 &middot; runs through Docker on a pinned R 4.4.3, a dated
CRAN snapshot and Bioconductor 3.20
</p>
</div>
</div>

<div class="cy-cards">
<div class="cy-card">
<h3>Gates derived, not transferred</h3>
<p>Each threshold is placed at a density minimum inside its own sample. Copying
one set of coordinates between samples does not remove operator variance, it
converts it into a bias that tracks staining intensity.</p>
</div>
<div class="cy-card">
<h3>Uncertainty that propagates</h3>
<p>Every cut is resampled from the events it came from and re-derived over the
settings that placed it. A frequency separated by a clean gap and one sitting on
a shoulder are not reported to the same precision.</p>
</div>
<div class="cy-card">
<h3>A specification that can be refuted</h3>
<p>Explore mode clusters every eligible channel without reference to what you
declared, which is the only way to find a population nobody named. Six
diagnostics exist to contradict the specification.</p>
</div>
</div>

<div class="cy-quote">
The unit of replication is the donor, not the event.
</div>

## Quick setup

Docker is the supported execution path, for a numerical reason rather than
convenience: `uwot` is a stochastic embedder whose output varies with its own
version and with the BLAS beneath it, and every gate is placed at a kernel
density minimum, so a different density implementation moves thresholds and
therefore the frequencies that would be published. The image pins R 4.4.3, a
dated CRAN snapshot and Bioconductor 3.20.

### 1. Get the image

Either pull the published one or build it from source. They are the same image;
pick by whether you intend to change the code.

**Pull.** A download rather than a compile, and the image every published run
here was produced with.

```bash
docker pull bhagesh/cyraven:1.0.0
```

Then give it the short local name every command below uses:

```bash
docker tag bhagesh/cyraven:1.0.0 cyraven:1.0.0
```

That second command is only for readability. Use `bhagesh/cyraven:1.0.0`
directly instead if you prefer.

**Build.** For a modified source tree, an air-gapped host, or if you would
rather not run a third-party image. The build context is the repository root,
the directory holding `DESCRIPTION`.

```bash
git clone https://github.com/bhagesh-h/cyRAVEN.git
```

```bash
cd cyRAVEN
```

```bash
docker build -f inst/scripts/Dockerfile -t cyraven:1.0.0 .
```

The first build compiles the dependency stack and takes 15 to 25 minutes. It
finishes by running `--help` and printing every package version, so a broken
image fails at build time rather than during an analysis.

Both routes pin R 4.4.3, the same dated CRAN snapshot and Bioconductor 3.20, so
they are numerically interchangeable. Confirm what you have:

```bash
docker run --rm --entrypoint Rscript cyraven:1.0.0 -e 'packageVersion("cyRAVEN")'
```

```
[1] '1.0.0'
```

### 2. Run the demonstration cohort

Nothing is downloaded: the data ship inside a package cyRAVEN already depends on,
so the example is reproducible offline and cannot break when a repository moves
or a certificate expires.

Four commands, one block each. Copy and run them one at a time.

**1. Make the folders.**

```bash
mkdir -p demo results
```

**2. Write the cohort, its sample sheet and its config.**

```bash
docker run --rm -v "$PWD/demo:/demo" \
 --entrypoint Rscript cyraven:1.0.0 \
  /opt/cyraven/src/inst/scripts/demo_data.R /demo
```

Ends with:

```
wrote samples.csv, panel.yaml and sample_map.csv to /demo
cohort: GvHD grade 1 n=7; GvHD grade 3 n=28
```

**3. Validate the inputs in seconds, without analysing anything.**

```bash
docker run --rm -v "$PWD/demo:/data:ro" -v "$PWD/results:/results" \
  cyraven:1.0.0 --dir /data/fcs \
 --samples /data/samples.csv --config /data/panel.yaml \
 --group-column cohort --reference-group "GvHD grade 1" \
 --batch-column visit --outdir /results --check
```

It should report 35 files, four resolved markers, that every marker the
specification names is present, that the sheet covers every file, and the group
sizes. Nothing is written and nothing is analysed. Fix anything it reports before
running step 4, which costs minutes rather than seconds.

**4. Run it.**

```bash
docker run --rm -v "$PWD/demo:/data:ro" -v "$PWD/results:/results" \
  cyraven:1.0.0 --dir /data/fcs \
 --samples /data/samples.csv --config /data/panel.yaml \
 --group-column cohort --reference-group "GvHD grade 1" \
 --batch-column visit --cluster --outdir /results
```

Every path inside a flag is a path inside the container, and `--outdir` must fall
within a mounted volume or the output is discarded when the container exits. On
Windows PowerShell substitute `${PWD}` for `$PWD`; on Git Bash prefix with
`MSYS_NO_PATHCONV=1`.

This writes 24 figures and 35 tables from 35 samples across five allogeneic
transplant recipients (Brinkman et al., 2007, *Biol Blood Marrow Transplant*
13:691; Artistic-2.0).

### 3. Read the result

Open `results/report.html`. It carries every figure and table the run produced,
embedded in one self-contained file that references nothing beside it, in the
order the outputs have to be read. Figures zoom and download at full resolution;
each table is its own toggle, carrying a line saying what it holds, and is
searchable, sortable, paged and exportable to CSV.

Inspect `recon_diagnostics.png` and `gating_qc.png` before any quantity derived
from them: a threshold placed on a distribution shoulder rather than a density
minimum is visible there and in no downstream table.

<img src="reference/figures/demo_gating_qc_example.png" alt="Two samples from the demonstration cohort, five markers each, with the derived threshold drawn on the distribution it came from" width="100%"/>

Two samples from that run, five markers each. The full figure carries one row per
sample; this is the top of it.

Each panel is one marker in one sample, and the dashed line is the cut that was
applied to it. Read it column by column:

- **`valley` in black** is a threshold placed at a real density minimum. `CD45` in
  the top row sits just right of the peak, separating the leukocytes from
  everything below them.
- **`quantile_fallback (REVIEW)` in red** means no minimum was found, so the cut
  came from a percentile instead. The second row has three of them. Those
  populations are still reported, with the fallback recorded in
  `thresholds_used.csv`, and they are the ones to check before quoting a number.

Comparing the two rows is the point. The `CD14` cut is in a different place in
each, because these are different samples with different staining, and the
correct cut genuinely differs. A fixed coordinate copied between them would be
right for at most one.

### 4. Your own data

Two files besides the FCS directory: one CSV with a row per file, and one YAML
declaring what to score.

**A sheet with a row for every file, which you then fill in.**

```bash
docker run --rm -v "$PWD/data:/data" cyraven:1.0.0 \
 --dir /data/fcs --recursive --write-samples /data/samples.csv
```

**The annotated config template.**

```bash
docker run --rm --entrypoint sh cyraven:1.0.0 -c \
  'cat /usr/local/lib/R/site-library/cyRAVEN/examples/analysis_template.yaml' \
  > data/analysis.yaml
```

Edit both, run `--check` until it reports no problems, then run. The complete
sequence for every task, in Docker and in R, with the situation each option is
for, is in
[Running cyRAVEN](https://bhagesh-h.github.io/cyRAVEN/articles/cyRAVEN.html#commands-and-every-option).

### A repeated-measures study

If the sample sheet carries a `timepoint` column, `--split-by-timepoint` writes
the whole figure set once per visit into `by_timepoint/<level>/`, beside the
pooled figures:

```bash
docker run --rm -v "$PWD/data:/data" cyraven:1.0.0 \
 --dir /data/fcs --recursive --samples /data/samples.csv \
 --config /data/analysis.yaml --split-by-timepoint --outdir /data/results
```

The embedding, the gating thresholds and the statistics are **not** recomputed
per visit -- only the rows drawn are restricted. A position on one visit's UMAP
is therefore the same position on another's, and a threshold is the same cut.
Running the pipeline separately per visit would give three embeddings computed
from different cells, which cannot be compared with each other at all.

Paired views that need every visit at once -- per-patient trajectories, subset
balance, marker intensity by visit -- are written to the run directory as
`timepoint_*.png` whether or not the flag is set, because splitting them by
visit is what destroys them.

## Documentation

Every page below is also available as a single PDF, one chapter per page, with a
table of contents and continuous page numbering:
**[cyRAVEN-manual.pdf](https://bhagesh-h.github.io/cyRAVEN/cyRAVEN-manual.pdf)**.
It is rebuilt with the site, so it never lags behind these pages. Every column of
every reference table is sized to the width its own content needs, so the wide
ones read upright rather than being squeezed until the columns touch, and each
table and figure stays on the same page as the heading above it.

New to flow cytometry?
[Flow cytometry for dummies](https://bhagesh-h.github.io/cyRAVEN/articles/flow-cytometry.html)
explains what an event is, what a cut is, why `SSC-A` sits in a population
definition next to two antibodies, and what the `-A` on a channel name means.
Every other page assumes it.

Before you trust a result,
[Known limitations](https://bhagesh-h.github.io/cyRAVEN/articles/limitations.html#known-limitations)
collects every caveat in one place: what the 1.0.0 control fix means for results
produced by an earlier version and how to tell whether yours were affected,
where the automated gates have and have not been validated, and which comparisons
this package will not make. It also records what was deliberately left out and
why, so the same arguments do not get relitigated.

The documentation is six chapters, read in order. Each one assumes the chapter
before it.

| Chapter | What it covers |
|---|---|
| [About cyRAVEN](https://bhagesh-h.github.io/cyRAVEN/articles/about.html) | What the pipeline does and why it is built this way: the variance manual gating introduces, the ten stages a run moves through, and why unsupervised discovery is kept at arm's length from the specification it checks |
| [About flow cytometry](https://bhagesh-h.github.io/cyRAVEN/articles/flow-cytometry.html) | What an event is, what a cut is, why `SSC-A` sits in a population definition beside two antibodies, and what the `-A` on a channel name means. Nothing in it is specific to cyRAVEN |
| [Get started](https://bhagesh-h.github.io/cyRAVEN/articles/cyRAVEN.html) | Install it, run the demonstration cohort, then build your own two input files. Carries the sample sheet and config reference, the gating specification, and every one of the 116 options |
| [Advanced](https://bhagesh-h.github.io/cyRAVEN/articles/advanced.html) | Beyond a first run: unsupervised discovery, the diagnostics in the order they have to be read, every output file, the statistics and their assumptions, and using cyRAVEN with cyCONDOR or from Claude Code |
| [Gallery](https://bhagesh-h.github.io/cyRAVEN/articles/gallery.html) | Every figure a run produces, as unmodified output, with what each one measures and what that run shows |
| [Limitations](https://bhagesh-h.github.io/cyRAVEN/articles/limitations.html) | Every caveat in one place, and nine excluded methods with the reasoning and the condition under which each becomes appropriate |
| [Function reference](https://bhagesh-h.github.io/cyRAVEN/reference/index.html) | Every exported function, grouped by stage |

**These pages are the website, not the install.** The package sources carry all
six chapters, but they are not built into the installed copy, so
`browseVignettes("cyRAVEN")` returns nothing inside the image and on a plain
install. Building them would add several megabytes of rendered HTML and figures
to a package whose footprint is already deliberate, and the site is always
current with `main`. Function help is installed and complete, `?run_cyraven`
and 192 other topics work offline. For the prose, use the links above or build
the vignettes yourself from a source checkout with
`devtools::build_vignettes()`.

One consequence worth knowing if you check the package yourself: because
`inst/doc` is absent, `R CMD check` reports two warnings,
`Files in the 'vignettes' directory but no files in 'inst/doc'` and
`Directory 'inst/doc' does not exist`. Both are consequences of not building
vignettes rather than defects, and every other check, including all 248 tests,
passes.

## Installation without Docker

Reproducibility of numerical output is not guaranteed outside the pinned
environment, for the reasons given under [Quick setup](#quick-setup).

```r
install.packages("BiocManager")
BiocManager::install("flowCore")
remotes::install_github("bhagesh-h/cyRAVEN")
```

```r
library(cyRAVEN)

run_cyraven(list(
  dir             = "demo/fcs",
  outdir          = "results",
  samples         = "demo/samples.csv",
  config          = "demo/panel.yaml",
  group_column    = "cohort",
  reference_group = "GvHD grade 1"
))
```

`FlowSOM` is an optional dependency. When absent, self-organising map clustering
falls back to the implementation within the package and all other behaviour is
unchanged.

## Citation

```r
citation("cyRAVEN")
```

Cite FlowSOM (Van Gassen et al., 2015, *Cytometry A* 87:636) when using the
unsupervised clustering, and diffcyt (Weber et al., 2019, *Commun Biol* 2:183)
for the aggregation strategy underlying the differential state tests.

## A note on how parts of this were written

Commit messages and parts of the documentation, this README, the vignettes and
the roxygen comments, were drafted with Anthropic's Claude, mostly Claude Opus 5
through [Claude Code](https://claude.com/claude-code), with earlier sections
written using Opus 4.1 and Sonnet 4.5.

The analysis code and its results are the author's, and every number quoted in
the documentation comes from a run rather than from a draft. Prose written this
way can still describe the code incorrectly. If something here does not match
what the package does, please
[open an issue](https://github.com/bhagesh-h/cyRAVEN/issues), a documentation
error is a bug and is worth reporting as one.

## Licence

GPL-3. The full text is in
[LICENSE](https://github.com/bhagesh-h/cyRAVEN/blob/main/LICENSE), and on the
documentation site under [License](https://bhagesh-h.github.io/cyRAVEN/LICENSE-text.html).
