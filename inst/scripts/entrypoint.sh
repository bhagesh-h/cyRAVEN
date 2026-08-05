#!/bin/sh
# Runs the cyRAVEN command-line front end -- against the copy of the package
# installed into the image at build time by default, or against a source tree
# mounted at $CYRAVEN_SOURCE instead. That indirection is what lets the image and
# the package be updated independently: mount your checkout read-only (the same
# device already used for --config/--sample-map/--patient-table, applied to the
# code as well) and every edit takes effect on the next `docker run`, with no
# rebuild in between.
#
#   # installed copy
#   docker run --rm -v "$PWD:/data:ro" -v "$PWD/results:/results" \
#     cyraven:0.1.0 --dir /data ...
#
#   # mounted source tree
#   docker run --rm -v "$PWD:/src:ro" -v "$PWD/results:/results" \
#     -e CYRAVEN_SOURCE=/src \
#     cyraven:0.1.0 --dir /data ...
#
# A rebuild is still required when the mounted source starts needing a package
# the image doesn't have -- no amount of mounting can install a package that
# isn't there. R says so plainly at startup ("there is no package called ...").
#
# WHY BOTH BRANCHES END IN `Rscript <file> "$@"` AND NOT `Rscript -e '...' "$@"`:
# two reasons, both learned the hard way.
#   1. R inserts its own `--args` separator ahead of the arguments that follow
#      `-e`, so an explicit one arrives at the front end as a literal option and
#      optparse rejects it ("long flag \"args\" is invalid").
#   2. Without a real script file there is no `--file=` in commandArgs(), so
#      optparse cannot name the program and its --help output reads `%prog`.
# Handing Rscript an actual file avoids both. The source branch gets its package
# loaded through R_PROFILE_USER, which runs before the script does, so the
# front end is invoked identically either way.
set -e

if [ -n "${CYRAVEN_SOURCE}" ]; then
  # pkgload::load_all() gives the mounted tree the same namespace semantics as an
  # installed package -- internal functions resolve, system.file() finds inst/ --
  # without writing to the image's library, which is read-only in most run
  # configurations.
  profile="$(mktemp)"
  cat > "$profile" <<'PROFILE'
local({
  src <- Sys.getenv("CYRAVEN_SOURCE")
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("CYRAVEN_SOURCE needs the 'pkgload' package in the image", call. = FALSE)
  pkgload::load_all(src, quiet = TRUE, helpers = FALSE)
})
PROFILE
  R_PROFILE_USER="$profile"
  export R_PROFILE_USER
  exec Rscript "${CYRAVEN_SOURCE}/inst/scripts/cyraven.R" "$@"
fi

CYRAVEN_CLI="$(Rscript -e 'cat(system.file("scripts", "cyraven.R", package = "cyRAVEN"))')"
if [ -z "${CYRAVEN_CLI}" ] || [ ! -f "${CYRAVEN_CLI}" ]; then
  echo "cyRAVEN is not installed in this image, and CYRAVEN_SOURCE is unset." >&2
  exit 1
fi
exec Rscript "${CYRAVEN_CLI}" "$@"
