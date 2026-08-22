# The PDF manual, as CI runs it.
#
# WHY THIS IS A SCRIPT AND NOT A `run:` BLOCK. Three builds produced no PDF and
# the only evidence available was a step duration, because a failing step's log
# needs a token to read and continue-on-error reports it as success either way.
# Everything here prints what it found, and the workflow tees the whole thing to
# docs/manual-build.log, which deploys with the site. That turns an opaque CI
# failure into a URL.
#
# It also keeps the quoting sane: an R heredoc inside YAML inside a shell was
# where two of the earlier attempts went wrong.

say <- function(...) cat(..., "\n", sep = "")

say("R:        ", R.version.string)
say("wd:       ", getwd())
say("pandoc:   ", Sys.which("pandoc"))

# CRAN rather than the snapshot the Bioconductor step pinned. The snapshot
# carries tinytex 0.56, which asks for an asset upstream has renamed.
if (!requireNamespace("tinytex", quietly = TRUE) ||
    utils::packageVersion("tinytex") < "0.60")
  utils::install.packages("tinytex", repos = "https://cloud.r-project.org")
say("tinytex:  ", as.character(utils::packageVersion("tinytex")))

say("is_tinytex before: ", tinytex::is_tinytex())
if (!tinytex::is_tinytex()) {
  # Pinned, not "daily". The daily asset was renamed and the download 404s.
  r <- tryCatch(tinytex::install_tinytex(version = "2026.08"),
                error = function(e) e)
  if (inherits(r, "error")) say("INSTALL ERROR: ", conditionMessage(r))
}
say("is_tinytex after:  ", tinytex::is_tinytex())
say("tlmgr:    ", Sys.which("tlmgr"))
say("pdflatex: ", Sys.which("pdflatex"))

# latexmk installs a missing style file on demand, so this is a pre-fetch and one
# unresolvable name must not stop the build before it has tried.
r <- tryCatch(tinytex::tlmgr_install(
  c("pdflscape", "titlesec", "fancyhdr", "fvextra", "xurl", "ragged2e",
    "booktabs", "framed", "environ", "trimspaces", "koma-script",
    # float pins each figure where it was written instead of letting it drift
    # away from the heading above it; see preamble.tex and page-fit.lua.
    "float")),
  error = function(e) e)
if (inherits(r, "error")) say("TLMGR ERROR: ", conditionMessage(r))

say("=== building ===")
r <- tryCatch(source("pkgdown/manual/build-manual.R"), error = function(e) e)
if (inherits(r, "error")) {
  say("BUILD ERROR: ", conditionMessage(r))
} else {
  f <- "docs/cyRAVEN-manual.pdf"
  say("result: ", if (file.exists(f))
        paste0(round(file.size(f) / 1024^2, 1), " MB") else "NO FILE")
}
