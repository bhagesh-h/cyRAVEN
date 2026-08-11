# Running cyRAVEN in Docker

The image is the reference environment. Everything below is the detail behind the
quickstart in SKILL.md.

## Why the container is the recommended path

Two dependencies decide numbers rather than convenience.

`uwot` is a stochastic embedder. The same `--seed` on a different `uwot` build
gives a visually similar but numerically different embedding.

Every gate comes from a kernel density estimate, so a different density
implementation moves a threshold, and a moved threshold changes the population
frequencies you would publish.

`rocker/r-ver:4.4.3` pins R and sets a dated CRAN snapshot through Posit Package
Manager, so `install.packages()` resolves to versions current at that R release
rather than today's. Bioconductor is pinned separately at 3.20 because `flowCore`
follows its own release cycle. `r-base:latest` would silently re-roll all of it
on a rebuild months later.

Versions the image reproduces:

```
R 4.4.3 · flowCore 2.18.0 · uwot 0.2.4 · data.table 1.18.4
ggplot2 4.0.3 · patchwork 1.3.2 · viridisLite 0.4.2 · scales 1.4.0
yaml 2.3.12 · optparse 1.8.2 · hexbin 1.28.5
```

## Build

From the package root, the directory containing `DESCRIPTION`, not
`inst/scripts/`:

```bash
docker build -f inst/scripts/Dockerfile -t cyraven:1.0.0 .
```

The build copies `DESCRIPTION` and `install_deps.R` before anything else, so the
expensive dependency layer caches independently of the analysis code. Editing
`R/` does not rebuild dependencies.

`install_deps.R` reads Imports and Suggests out of `DESCRIPTION` rather than
repeating the list, so the image cannot install a different set of packages from
the one the package declares.

The final build step runs `--help`, which exercises the whole option parser and
every dependency the CLI touches, then prints each package version. A broken
image fails here.

## Run

```bash
docker run --rm \
  -v "$PWD:/data:ro" \
  -v "$PWD/results:/results" \
  cyraven:1.0.0 \
  --dir /data --recursive \
  --config /data/panel.yaml \
  --sample-map /data/sample_map.csv \
  --patient-table /data/patient_table.csv \
  --group-column cohort --reference-group "Healthy controls" \
  --outdir /results
```

Points that catch people out:

- Every path in a flag is a **container** path. `--config /data/panel.yaml` is
  the mount, not the host directory.
- `--outdir /results` must be inside a mounted volume or the output disappears
  with the container. `/results` is declared as a `VOLUME` in the image, so
  omitting the mount silently writes into an anonymous volume.
- The data mount is read-only (`:ro`). Keep it that way; the pipeline never needs
  to write beside the FCS files.
- `--rm` removes the container on exit. Drop it only if you want to inspect the
  filesystem afterwards.

### Windows

PowerShell expands `${PWD}`, not `$PWD`:

```powershell
docker run --rm `
  -v "${PWD}:/data:ro" `
  -v "${PWD}/results:/results" `
  cyraven:1.0.0 `
  --dir /data --recursive --outdir /results
```

Git Bash on Windows rewrites paths beginning with `/`, which turns `/data` into
`C:/Program Files/Git/data`. Either use PowerShell, or prefix with a second
slash: `-v "$PWD://data:ro"`, or set `MSYS_NO_PATHCONV=1`.

Docker Desktop must have the drive shared under Settings, Resources, File
sharing, or the mount silently produces an empty directory.

## Resource tuning

`OMP_NUM_THREADS` defaults to 4 in the image. `uwot` and `data.table` are
OpenMP-threaded and would otherwise take every core on the host, which on a
shared machine is antisocial and on a laptop causes thermal throttling that makes
the run slower.

```bash
-e OMP_NUM_THREADS=8            # container-wide
--threads 8                     # the pipeline's own flag, UMAP only
```

Memory is the usual constraint on large cohorts. Each sample's raw expression
matrix is released as soon as its transformed matrix exists, but the embedding
still holds up to `--max-cells` rows. If the container is killed without an R
error, it was the OOM reaper:

```bash
--memory=16g                    # docker flag, makes the limit explicit
--max-events-per-file 300000    # pipeline flag, bounds peak memory
--max-cells 100000              # pipeline flag, bounds the embedding
```

`--max-events-per-file` samples evenly through the acquisition rather than taking
the first N, so it does not bias toward the start of the tube. It does lower
every parent-gate event count, which raises the detection limit of every rare
population proportionally.

## Iterating without rebuilding

```bash
docker run --rm \
  -v "$PWD:/src:ro" \
  -v "$PWD/data:/data:ro" \
  -v "$PWD/results:/results" \
  -e CYRAVEN_SOURCE=/src \
  cyraven:1.0.0 --dir /data --outdir /results
```

`entrypoint.sh` sees `CYRAVEN_SOURCE` and loads that tree with
`pkgload::load_all()` through `R_PROFILE_USER`, which runs before the script. The
mounted tree gets the same namespace semantics as an installed package (internal
functions resolve, `system.file()` finds `inst/`) without writing to the image's
library, which is read-only in most run configurations.

A rebuild is still required once the mounted source needs a package the image
does not have. R says so plainly at startup: "there is no package called ...".

## Docker-specific failures

**`Rscript: not a valid option`, or the container ignores your command.**
The image has an `ENTRYPOINT`, so `docker run cyraven:1.0.0 Rscript foo.R` passes
`Rscript foo.R` as *arguments to the entrypoint*. Override it:

```bash
docker run --rm --entrypoint Rscript cyraven:1.0.0 -e 'sessionInfo()'
docker run --rm -it --entrypoint bash cyraven:1.0.0
```

**`R CMD INSTALL -l /tmp/rlib: cannot cd to directory`.**
The library path must exist first. `mkdir -p /tmp/rlib` before installing.

**The run seems to hang, or is far slower than expected.**
Check for an orphaned container competing for CPU. Killing the shell that started
a `docker run` does not stop the container:

```bash
docker ps
docker kill <container-id>
```

**Output directory is empty after a successful run.**
`--outdir` pointed somewhere outside a mount. Compare the `--outdir` value
against the right-hand side of every `-v`.

**Permission denied writing results on Linux.**
The container writes as root by default, or as a UID with no rights to the
mounted directory. Either `--user "$(id -u):$(id -g)"` and pre-create the results
directory, or `chown` it afterwards.

## Checking what an image actually contains

```bash
docker run --rm --entrypoint Rscript cyraven:1.0.0 \
  -e 'cat(as.character(packageVersion("cyRAVEN")), "\n")'

docker run --rm cyraven:1.0.0 --help
```

The manifest written by every run (`run_manifest.txt`) records the R version, the
platform, the version of every package loaded at run time, and the git commit
where the code was a checkout. That file, not the image tag, is what ties a
results folder to the code that produced it.
