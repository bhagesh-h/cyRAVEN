# cyRAVEN <img src="man/figures/logo.png" alt="cyRAVEN logo" width="200" align="right"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/bhagesh-h/cyRAVEN/actions/workflows/pkgdown.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Documentation: <https://bhagesh-h.github.io/cyRAVEN/>**

An R package for supervised immunophenotyping of multi-sample flow cytometry
data, with per-sample gate derivation and donor-level differential abundance and
differential state testing.

A directory of FCS files becomes population frequencies, marker expression, a
shared UMAP and between-group statistics. Every threshold is placed within the
sample it applies to, every frequency carries the uncertainty of the cut and of
the event count behind it, and six diagnostics test whether the gating
specification held.

## Why

Manual gating is the dominant source of technical variance in multi-sample
immunophenotyping: operators gating identical files report population sizes
differing by approximately 32%, and analyst subjectivity accounts for up to 78%
of technical variability once more than one person is involved (Cadwell et al.,
2021, *PDA J Pharm Sci Technol* 75:33). Fixed gate coordinates transferred
between samples do not remove this variance; they convert it into a systematic
bias that tracks staining intensity.

cyRAVEN removes the analyst from threshold placement and quantifies the residual
uncertainty rather than concealing it. Each cut is resampled from the events it
was derived from and re-derived over the settings that placed it, and the
resulting spread is propagated to every population that reads it, so a frequency
separated by a clean gap and one sitting on a shoulder are not reported to the
same apparent precision.

One constraint governs the design:

> The unit of replication is the donor, not the event.

Event counts are set by acquisition duration rather than by study design. Every
test aggregates to one value per sample before estimation; quantities that cannot
be aggregated are reported with effect sizes and no p-value.

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

**Pull** — a download rather than a compile, and the image every published run
here was produced with.

```bash
docker pull bhagesh/cyraven:1.0.0
docker tag bhagesh/cyraven:1.0.0 cyraven:1.0.0
```

The second line is only so that every command below can say `cyraven:1.0.0`. Use
`bhagesh/cyraven:1.0.0` directly instead if you prefer.

**Build** — for a modified source tree, an air-gapped host, or if you would
rather not run a third-party image. The build context is the repository root,
the directory holding `DESCRIPTION`.

```bash
git clone https://github.com/bhagesh-h/cyRAVEN.git
cd cyRAVEN
docker build -f inst/scripts/Dockerfile -t cyraven:1.0.0 .
```

The first build compiles the dependency stack and takes 15 to 25 minutes. It
finishes by running `--help` and printing every package version, so a broken
image fails at build time rather than during an analysis.

Both routes pin R 4.4.3, the same dated CRAN snapshot and Bioconductor 3.20, so
they are numerically interchangeable. Confirm what you have with
`docker run --rm --entrypoint Rscript cyraven:1.0.0 -e 'packageVersion("cyRAVEN")'`.

### 2. Run the demonstration cohort

Nothing is downloaded: the data ship inside a package cyRAVEN already depends on,
so the example is reproducible offline and cannot break when a repository moves
or a certificate expires.

```bash
mkdir -p demo results

# write the cohort, its sample sheet and its config
docker run --rm -v "$PWD/demo:/demo" \
  --entrypoint Rscript cyraven:1.0.0 \
  /opt/cyraven/src/inst/scripts/demo_data.R /demo

# validate the inputs in seconds, without analysing anything
docker run --rm -v "$PWD/demo:/data:ro" -v "$PWD/results:/results" \
  cyraven:1.0.0 --dir /data/fcs \
  --samples /data/samples.csv --config /data/panel.yaml \
  --group-column cohort --reference-group "GvHD grade 1" \
  --batch-column visit --outdir /results --check

# run
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

This writes 22 figures and 31 tables from 35 samples across five allogeneic
transplant recipients (Brinkman et al., 2007, *Biol Blood Marrow Transplant*
13:691; Artistic-2.0).

### 3. Read the result

Open `results/report.html`. It carries every figure and table the run produced,
embedded in one self-contained file that references nothing beside it, in the
order the outputs have to be read. Figures zoom and download at full resolution;
tables are searchable, sortable, paged and exportable to CSV.

Inspect `recon_diagnostics.png` and `gating_qc.png` before any quantity derived
from them: a threshold placed on a distribution shoulder rather than a density
minimum is visible there and in no downstream table.

<img src="man/figures/demo_gating_qc.png" alt="Per-sample marker densities with the derived threshold marked on each" width="100%"/>

### 4. Your own data

Two files besides the FCS directory: one CSV with a row per file, and one YAML
declaring what to score.

```bash
# a sheet with a row for every file, which you then fill in
docker run --rm -v "$PWD/data:/data" cyraven:1.0.0 \
  --dir /data/fcs --recursive --write-samples /data/samples.csv

# the annotated config template
docker run --rm --entrypoint sh cyraven:1.0.0 -c \
  'cat /usr/local/lib/R/site-library/cyRAVEN/examples/analysis_template.yaml' \
  > data/analysis.yaml
```

Edit both, run `--check` until it reports no problems, then run. The complete
sequence for every task, in Docker and in R, with the situation each option is
for, is in
[Running cyRAVEN](https://bhagesh-h.github.io/cyRAVEN/articles/usage.html).

## Documentation

**Code and setup**

| Page | Content |
|---|---|
| [Running cyRAVEN](https://bhagesh-h.github.io/cyRAVEN/articles/usage.html) | Every command in Docker and R, which variation to use in which situation, and the complete reference for all 83 options |
| [Inputs](https://bhagesh-h.github.io/cyRAVEN/articles/inputs.html) | The sample sheet and the config: every column, the templates, and the errors the format can raise |
| [Gating specification](https://bhagesh-h.github.io/cyRAVEN/articles/gating.html) | Declaring populations, functional blocks and ratios; per-sample thresholding and the `source` column; arcsinh against logicle |
| [Claude skill](https://bhagesh-h.github.io/cyRAVEN/articles/claude-skill.html) | Installing and using the bundled Claude Code skill, which executes through Docker by default |
| [Interoperability](https://bhagesh-h.github.io/cyRAVEN/articles/with-cycondor.html) | Method selection by question type, and handing a cyCONDOR clustering to cyRAVEN to obtain an executable gate |
| [Function reference](https://bhagesh-h.github.io/cyRAVEN/reference/index.html) | Every exported function, grouped by stage |

**Science and output**

| Page | Content |
|---|---|
| [How it works](https://bhagesh-h.github.io/cyRAVEN/articles/pipeline.html) | The ten pipeline stages with the function implementing each, and executable examples of cofactor estimation, density minimum detection and group comparison |
| [Diagnostics](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.html) | The checks in reading order, from validating inputs before the run through gate inspection, staining QC, threshold drift, gate uncertainty, batch structure and conformance |
| [Output files](https://bhagesh-h.github.io/cyRAVEN/articles/outputs.html) | Every file the pipeline writes, its columns, the flag producing it, and why event counts are not cell counts |
| [Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.html) | Sample-level aggregation; rank tests against moderated *t*; the compositional constraint; covariate diagnosis against adjustment; multiplicity; differences expressed in units of uncertainty |
| [Worked example](https://bhagesh-h.github.io/cyRAVEN/articles/figures.html) | All 22 figures from a run on public data, each with what it measures and what that run shows |
| [Scope](https://bhagesh-h.github.io/cyRAVEN/articles/scope.html) | Nine excluded methods with the reasoning and the condition under which each becomes appropriate |

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

## Licence

GPL-3. See [LICENSE](LICENSE).
