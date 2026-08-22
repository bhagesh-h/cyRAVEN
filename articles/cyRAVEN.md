# Get started

## What cyRAVEN does

A directory of FCS files becomes population frequencies, marker
expression, a shared UMAP and between-group statistics. Populations are
declared in advance; every marker threshold is derived inside the sample
it applies to; and the declaration is then tested against unsupervised
structure recovered from the same cells.

There are two ways to use it, and they answer different questions.

|  | **Supervised**, the default | **Explore mode**, `--explore` |
|----|----|----|
| You supply | A specification of the populations to score | Nothing beyond the FCS files |
| It returns | Named populations, with uncertainty and group statistics | Numbered clusters over every eligible channel |
| It can find | Only what you declared | A population nobody declared |
| Runs on | The parent gate | Every event in the file |

Explore mode writes everything into `<outdir>/explore/` and changes no
existing output, so adding `--explore` to a run is always safe.
`--explore-only` skips the declared analysis entirely and needs no
specification, which is the mode for a panel you have not written one
for yet. The two stay isolated unless you pass `--maybe-learn`; see
[Explore mode](https://bhagesh-h.github.io/cyRAVEN/articles/explore.md)
for what that exchanges, and [the design
record](https://bhagesh-h.github.io/cyRAVEN/articles/design-explore.md)
for why it is off by default.

The design addresses one failure mode. Manual gating is the dominant
source of technical variance in multi-sample immunophenotyping:
operators gating identical files report population sizes differing by
about 32%, and analyst subjectivity accounts for up to 78% of technical
variability once more than one person is involved (Cadwell et al., 2021,
*PDA J Pharm Sci Technol* 75:33). Transferring fixed gate coordinates
between samples does not remove that variance; it converts it into a
bias that tracks staining intensity.

cyRAVEN removes the analyst from threshold placement and quantifies what
remains. Each cut is resampled from the events it was derived from and
re-derived across the settings that placed it, and the resulting spread
is propagated to every population that reads it, so a frequency
separated by a clean gap and one sitting on a shoulder are not reported
to the same apparent precision.

One constraint governs everything below.

> The unit of replication is the donor, not the event.

Event counts are set by acquisition duration rather than by study
design. Every test aggregates to one value per sample before estimation;
quantities that cannot be aggregated are reported with effect sizes and
no p-value.

## Setup

Docker is the supported path, for a numerical reason rather than
convenience: `uwot` is a stochastic embedder whose output varies with
its own version and with the BLAS beneath it, and every gate is placed
at a kernel density minimum, so a different density implementation moves
thresholds and therefore the frequencies that would be published. The
image pins R 4.4.3, a dated CRAN snapshot and Bioconductor 3.20.

Two ways to get it. They produce the same image; pick by whether you
intend to change the code.

**Pull.** A download rather than a compile, and the image the worked
example and every published run here was produced with.

``` bash
docker pull bhagesh/cyraven:1.0.0
```

``` bash
docker tag bhagesh/cyraven:1.0.0 cyraven:1.0.0
```

The `docker tag` line only exists so every command on this site can say
`cyraven:1.0.0`. Substitute `bhagesh/cyraven:1.0.0` throughout instead
if you prefer.

**Build.** For a modified source tree, an air-gapped host, or if you
would rather not run a third-party image.

``` bash
git clone https://github.com/bhagesh-h/cyRAVEN.git
```

``` bash
cd cyRAVEN
```

``` bash
docker build -f inst/scripts/Dockerfile -t cyraven:1.0.0 .
```

The build context is the repository root, the directory holding
`DESCRIPTION`. The first build compiles the dependency stack and takes
15 to 25 minutes. It ends by running `--help` and printing every package
version, so a broken image fails at build time rather than during an
analysis.

Both pin R 4.4.3, the same dated CRAN snapshot and Bioconductor 3.20, so
they are numerically interchangeable. Check what you have with:

``` bash
docker run --rm --entrypoint Rscript cyraven:1.0.0 \
  -e 'cat(as.character(packageVersion("cyRAVEN")), R.version.string)'
```

Prints:

    1.0.0 R version 4.4.3 (2025-02-28)

Local installation, for development or where Docker is unavailable, is
in the
[README](https://github.com/bhagesh-h/cyRAVEN#installation-without-docker).

## A complete run in five minutes

No data of your own required, and nothing downloaded: the demonstration
cohort ships inside a package cyRAVEN already depends on.

Four commands, one block each. Copy and run them one at a time.

**1. Make the folders.**

``` bash
mkdir -p demo results
```

**2. Write the cohort, its sample sheet and its config.**

``` bash
docker run --rm -v "$PWD/demo:/demo" \
  --entrypoint Rscript cyraven:1.0.0 \
  /opt/cyraven/src/inst/scripts/demo_data.R /demo
```

**3. Validate the inputs, in seconds, without analysing anything.**

``` bash
docker run --rm -v "$PWD/demo:/data:ro" -v "$PWD/results:/results" \
  cyraven:1.0.0 --dir /data/fcs \
  --samples /data/samples.csv --config /data/panel.yaml \
  --group-column cohort --reference-group "GvHD grade 1" \
  --batch-column visit --outdir /results --check
```

**4. Run it, then open `results/report.html`.**

``` bash
docker run --rm -v "$PWD/demo:/data:ro" -v "$PWD/results:/results" \
  cyraven:1.0.0 --dir /data/fcs \
  --samples /data/samples.csv --config /data/panel.yaml \
  --group-column cohort --reference-group "GvHD grade 1" \
  --batch-column visit --cluster --outdir /results
```

On Windows PowerShell substitute `${PWD}` for `$PWD`; on Git Bash prefix
with `MSYS_NO_PATHCONV=1`.

This writes 24 figures and 35 tables from 35 samples across five
allogeneic transplant recipients. `report.html` carries all of them
embedded in one self-contained file, in the order they have to be read.

## Your own data

Two files besides the FCS directory: one CSV with a row per file, and
one YAML declaring what to score.

A sheet with a row for every file, which you then fill in:

``` bash
docker run --rm -v "$PWD/data:/data" cyraven:1.0.0 \
  --dir /data/fcs --recursive --write-samples /data/samples.csv
```

The config template:

``` bash
docker run --rm --entrypoint sh cyraven:1.0.0 -c \
  'cat /usr/local/lib/R/site-library/cyRAVEN/examples/analysis_template.yaml' \
  > data/analysis.yaml
```

Then `--check` until it reports no problems, then run. The full
sequence, and every option with the situation it is for, is in [Commands
and every
option](https://bhagesh-h.github.io/cyRAVEN/articles/usage.md).

## Where to go next

New to flow cytometry? Start with [Flow cytometry for
dummies](https://bhagesh-h.github.io/cyRAVEN/articles/flow-cytometry-for-dummies.md),
which explains what an event is, what a cut is, why `SSC-A` appears in a
population definition next to two antibodies, and what the `-A` on a
channel name means. Every other page assumes it.

**Code and setup**

- [Running
  cyRAVEN](https://bhagesh-h.github.io/cyRAVEN/articles/usage.md): every
  command in Docker and R, which option to use when, and the complete
  option reference.
- [Inputs](https://bhagesh-h.github.io/cyRAVEN/articles/inputs.md): the
  sample sheet, the config, and the templates.
- [Gating
  specification](https://bhagesh-h.github.io/cyRAVEN/articles/gating.md):
  declaring populations, functional blocks and ratios.

**Science and output**

- [How it
  works](https://bhagesh-h.github.io/cyRAVEN/articles/pipeline.md): the
  ten stages and the function implementing each.
- [Diagnostics](https://bhagesh-h.github.io/cyRAVEN/articles/diagnostics.md):
  the checks that can invalidate a result, in the order they have to be
  read.
- [Explore
  mode](https://bhagesh-h.github.io/cyRAVEN/articles/explore.md):
  unsupervised discovery, run beside the declared analysis or on its
  own, and what `--maybe-learn` lets the two tell each other.
- [Why explore mode is built that
  way](https://bhagesh-h.github.io/cyRAVEN/articles/design-explore.md):
  the design record, and where cyRAVEN sits against FlowSOM, cyCONDOR
  and diffcyt.
- [Output
  files](https://bhagesh-h.github.io/cyRAVEN/articles/outputs.md): every
  file a run writes, with its columns.
- [Statistics](https://bhagesh-h.github.io/cyRAVEN/articles/statistics.md):
  test selection and its assumptions.
- [Worked
  example](https://bhagesh-h.github.io/cyRAVEN/articles/figures.md):
  every figure from the demonstration cohort, with its interpretation.
- [Scope](https://bhagesh-h.github.io/cyRAVEN/articles/scope.md): what
  cyRAVEN does not do, and what to use instead.
- [Known
  limitations](https://bhagesh-h.github.io/cyRAVEN/articles/known-limitations.md):
  every caveat collected in one place, how to tell whether a result from
  an earlier version was affected by the control fix, and what was
  deliberately left out with the reason.
