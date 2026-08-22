<!--
Source of truth for the Docker Hub description at
https://hub.docker.com/r/bhagesh/cyraven

Two fields on Docker Hub:
  Short description  the one line below the repository name, 100 characters
  Full description   everything under the "Overview" heading, Markdown

To update: edit this file, then either paste the two sections into the web UI, or
PATCH them (needs jq and a Docker Hub personal access token in $DH_TOKEN):

  jq -n --arg d "$SHORT" --arg f "$FULL" '{description:$d, full_description:$f}' \
    | curl -sS -X PATCH -H "Authorization: Bearer $DH_TOKEN" \
        -H "Content-Type: application/json" --data-binary @- \
        https://hub.docker.com/v2/repositories/bhagesh/cyraven/

Keep it in step with README.md. The short description has a 100-character limit.

Every code block holds ONE command. A reader copies a block whole, and a block
holding four commands runs all four back to back with no chance to read what the
first one said.
-->

## Short description

Supervised flow-cytometry immunophenotyping: per-sample gates, uncertainty, donor-level statistics.

## Full description

# cyRAVEN

Supervised immunophenotyping for multi-sample flow cytometry. A directory of FCS
files becomes population frequencies, marker expression, a shared UMAP and
between-group statistics.

Every threshold is placed inside the sample it applies to, at that sample's own
density minimum, rather than copied from a template. Every frequency carries the
uncertainty of the cut and of the event count behind it. Six diagnostics test
whether the gating specification held.

- **Source:** https://github.com/bhagesh-h/cyRAVEN
- **Documentation:** https://bhagesh-h.github.io/cyRAVEN/
- **Manual as one PDF:** https://bhagesh-h.github.io/cyRAVEN/cyRAVEN-manual.pdf
- **Licence:** GPL-3

## Tags

| Tag | Contents |
|---|---|
| `1.0.0` | cyRAVEN 1.0.0. Pin this for anything you intend to publish |
| `latest` | The most recent release. Currently identical to `1.0.0` |

`linux/amd64`, about 675 MB compressed.

### If you pulled this image before, re-pull it

The `1.0.0` tag has been republished. An earlier image under the same tag read an
**absent** `is_control` column in the sample sheet as "this file might be an
unstained control", because only `is_control = FALSE` positively asserts that a
file is a biological sample. On a sheet without that column, any sample whose
CD45 gate found no density minimum was promoted to the unstained reference for
its whole panel, and every threshold there became the 99.5th percentile of that
one sample. Where it was really a stained sample, every population collapsed,
with no error raised and no warning printed.

**The version number cannot tell you which image you have, so check the output.**
Open `thresholds_used.csv` from the run. If the `source` column reads
`control_q995` on most markers, those thresholds came from a reference; find out
which sample supplied it before reading a single frequency. A cohort that
declared no control should show only `valley` and `quantile_fallback`.

The current image believes a control only when the sample sheet declares one, and
names any sample it declines to promote. `--discover-controls` restores the old
behaviour for anyone who needs it.

## Quick start

Five commands, one block each. Copy and run them one at a time.

**1. Get the image.**

```bash
docker pull bhagesh/cyraven:1.0.0
```

**2. See what the flags are.**

```bash
docker run --rm bhagesh/cyraven:1.0.0 --help
```

**3. Make the folders.**

```bash
mkdir -p demo results
```

**4. Write the demonstration cohort, its sample sheet and its config.**

```bash
docker run --rm -v "$PWD/demo:/demo" --entrypoint Rscript bhagesh/cyraven:1.0.0 /opt/cyraven/src/inst/scripts/demo_data.R /demo
```

Ends with:

```
wrote samples.csv, panel.yaml and sample_map.csv
cohort: GvHD grade 1 n=7; GvHD grade 3 n=28
```

**5. Validate the inputs from the FCS headers alone.** Seconds, and analyses
nothing. It reports the markers resolved, any specification entry matching none
of them, whether the sheet covers every file, and the group sizes. Fix anything
it names before the next step, which costs minutes rather than seconds.

```bash
docker run --rm -v "$PWD/demo:/data:ro" -v "$PWD/results:/results" bhagesh/cyraven:1.0.0 --dir /data/fcs --samples /data/samples.csv --config /data/panel.yaml --group-column cohort --reference-group "GvHD grade 1" --batch-column visit --outdir /results --check
```

**6. Run it.**

```bash
docker run --rm -v "$PWD/demo:/data:ro" -v "$PWD/results:/results" bhagesh/cyraven:1.0.0 --dir /data/fcs --samples /data/samples.csv --config /data/panel.yaml --group-column cohort --reference-group "GvHD grade 1" --batch-column visit --cluster --outdir /results
```

That writes **24 figures and 35 tables** at the top level, full-size UMAPs per
marker in `marker_umaps_by_group/` (each marker pooled, and split by every
category the sample sheet carries), and a self-contained `results/report.html`:
one file carrying every figure and table, linking to nothing and loading nothing
from a network. Each table in the report is its own toggle, carrying a line
saying what it holds. Nothing is downloaded at any point; the demonstration
cohort ships inside a package cyRAVEN already depends on.

## Your own data

Two files beside the FCS directory: one CSV with a row per file, one YAML
declaring what to score.

**A sheet listing every file, which you then fill in.**

```bash
docker run --rm -v "$PWD/data:/data" bhagesh/cyraven:1.0.0 --dir /data/fcs --recursive --write-samples /data/samples.csv
```

**The annotated config template.**

```bash
docker run --rm --entrypoint sh bhagesh/cyraven:1.0.0 -c 'cat /usr/local/lib/R/site-library/cyRAVEN/examples/analysis_template.yaml' > data/analysis.yaml
```

Edit both, run with `--check` until it reports no problems, then run.

**Run `--list-channels` before you write the config.** It prints what each
channel resolves to, which is the name a population definition has to use, and on
a spectral panel that differs from what the file stores: `$PnS` reads
`CD45 : SparkUV-387 - Area` while the run uses `CD45`. Writing the raw string
produces a table of zeros.

```bash
docker run --rm -v "$PWD/data:/data:ro" bhagesh/cyraven:1.0.0 --dir /data/fcs --recursive --list-channels
```

It reads every file rather than the first, so it also catches a cohort splitting
into separate panels, which happens silently, because every file still loads and
every table is still written. Where the files disagree it names the fix rather
than the symptom, separating autofluorescence and unmixing artefacts from a real
marker stained in only some files, and printing the exact `--ignore-channels`
pattern that merges the cohort.

Every path inside a flag is a path **inside the container**, and `--outdir` must
fall within a mounted volume or the output is discarded when the container exits.
On Windows PowerShell use `${PWD}` for `$PWD`; on Git Bash prefix the command
with `MSYS_NO_PATHCONV=1`.

Full command reference, with the situation each option is for:
https://bhagesh-h.github.io/cyRAVEN/articles/usage.html

## Clinical variables

A severity score, a laboratory value or an outcome flag is neither the study
group nor a confounder. A confounder is screened to decide whether a group
difference can be believed; a clinical variable is the question. Name any sheet
column and it is tested against every population and every marker.

```bash
docker run --rm -v "$PWD/data:/data:ro" -v "$PWD/results:/results" bhagesh/cyraven:1.0.0 --dir /data/fcs --samples /data/samples.csv --config /data/analysis.yaml --clinical-columns sofa,survival_28d,crp_mg_l --outdir /results
```

The test follows the column's type: Spearman for numeric, Wilcoxon with Cliff's
delta for two levels, Kruskal-Wallis with epsilon-squared for more. Each variable
is tested on the unit it is a property of, so a score recorded at every draw is
tested per sample while a property of the person is collapsed to one value per
patient first. Counting three samples from one patient as three observations
inflates n and narrows the interval without one extra patient being recruited.

Every signed effect carries a percentile bootstrap interval. On ten patients most
of those intervals span zero, and that is the result rather than a shortcoming of
the figure.

## Why an image rather than an R install

Numerical, not convenience. `uwot` is a stochastic embedder whose output varies
with its own version and with the BLAS beneath it, and every gate is placed at a
kernel density minimum, so a different density implementation moves thresholds
and therefore the frequencies that would be published.

This image pins **R 4.4.3**, a dated CRAN snapshot and **Bioconductor 3.20**.
Every run writes a manifest recording package versions, the git commit, the full
invocation and every input file with its size and modification time.

cyRAVEN also installs as a plain R package, and that is documented, but
reproducibility of numerical output is only guaranteed inside this image.

Check what you have:

```bash
docker run --rm --entrypoint Rscript bhagesh/cyraven:1.0.0 -e 'cat(as.character(packageVersion("cyRAVEN")), R.version.string)'
```

```
1.0.0 R version 4.4.3 (2025-02-28)
```

## Building it yourself

Pull and build give the same image. Build when you have changed the source, are
on an air-gapped host, or would rather not run a third-party image.

```bash
git clone https://github.com/bhagesh-h/cyRAVEN.git
```

```bash
cd cyRAVEN
```

```bash
docker build -f inst/scripts/Dockerfile -t cyraven:1.0.0 .
```

The build context is the repository root, the directory holding `DESCRIPTION`.
15 to 25 minutes on a first build; it ends by running `--help` and printing every
package version, so a broken image fails at build time rather than during an
analysis.

## Notes

- **Run one container at a time.** The stages are already multi-threaded, so
  nothing is gained by running two, and a run killed by the out-of-memory killer
  writes no report at all.
- A run that **fails** still writes `report.html`, with the diagnosis first: what
  stopped it, which stage it reached, the log leading to that point, what the
  message means for the data, and every output produced before the failure.
- The whole documentation set is also one PDF, one chapter per page, every
  reference table sized to fit an upright page and every figure on the same page
  as its own heading:
  https://bhagesh-h.github.io/cyRAVEN/cyRAVEN-manual.pdf
