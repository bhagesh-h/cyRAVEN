# =============================================================================
# THE WHOLE DOCUMENTATION SITE AS ONE PDF
# =============================================================================
#
# WHY THIS EXISTS. The site is fifteen articles that are read in an order, and a
# reader who wants the whole thing -- to read offline, to hand to a collaborator,
# to attach to a submission -- has no way to get it. Fifteen browser prints
# produce fifteen files with fifteen sets of page numbers and no cross-references.
#
# WHAT IT PRODUCES. `cyRAVEN-manual.pdf` in docs/, one chapter per site page, in
# the reading order the navbar uses rather than alphabetically. The site links to
# it, so it is downloadable from the same place it is generated.
#
# HOW IT ASSEMBLES. Every vignette is knitted to markdown FIRST, in the package's
# own environment, so the chunks that evaluate (statistics, pipeline) run exactly
# as they do for the site. The markdown files are then concatenated with a
# preamble and rendered once by pandoc. Rendering each vignette separately and
# merging PDFs would lose the continuous page numbering and the single table of
# contents, which are the two reasons to want one file.
#
# WHY THE HEADINGS ARE SHIFTED. Each vignette starts its own headings at `#`.
# Concatenated unchanged that yields fifteen level-one headings and no chapter
# structure, so every `#` in a vignette body is pushed down one level and the
# chapter title is inserted above it. The article's own title becomes the
# chapter.

CHAPTERS <- c(
  # Reading order, matching the navbar. Not alphabetical and not the order
  # list.files() returns: the site tells a reader to read diagnostics before
  # results, and a manual that reversed them would contradict it.
  "cyRAVEN"                    = "Get started",
  "flow-cytometry-for-dummies" = "Flow cytometry for dummies",
  "usage"                      = "Commands and every option",
  "inputs"                     = "Inputs, the sample sheet and config",
  "gating"                     = "Gating specification",
  "claude-skill"               = "Driving it from Claude Code",
  "with-cycondor"              = "Using it with cyCONDOR",
  "pipeline"                   = "How it works, the ten stages",
  "diagnostics"                = "Diagnostics, in reading order",
  "explore"                    = "Explore mode, unsupervised discovery",
  "design-explore"             = "Why explore mode is built that way",
  "outputs"                    = "Every output file",
  "statistics"                 = "Statistics",
  "figures"                    = "Worked example, every figure",
  "scope"                      = "Scope and boundaries",
  "known-limitations"          = "Known limitations")

pkg_root <- normalizePath(".", winslash = "/")
vig_dir  <- file.path(pkg_root, "vignettes")
out_dir  <- file.path(pkg_root, "docs")
work     <- file.path(tempdir(), "manual")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(cyRAVEN))

# Figures land beside the markdown, and pandoc resolves relative paths from the
# document's directory, so knitting happens with the working directory set to
# vignettes/ and the images are copied next to the assembled file.
strip_yaml <- function(lines) {
  if (length(lines) && grepl("^---\\s*$", lines[1])) {
    close_at <- which(grepl("^---\\s*$", lines))[2]
    if (!is.na(close_at)) lines <- lines[-seq_len(close_at)]
  }
  lines
}

# Push every heading down one level so the chapter title is the only `#`.
# Fenced code blocks are skipped: a shell comment `# 1. write the cohort` is not
# a heading, and demoting it would corrupt the command.
demote <- function(lines) {
  in_fence <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("^\\s*```", lines[i])) in_fence <- !in_fence
    if (!in_fence && grepl("^#{1,5} ", lines[i])) lines[i] <- paste0("#", lines[i])
  }
  lines
}

# CROSS-REFERENCES BECOME INTERNAL LINKS. The articles link to each other by
# file name -- "see the [Diagnostics article](diagnostics.html)" -- which is
# correct on the site and dead in a PDF, where diagnostics.html does not exist.
# There are around sixty of them. Each is rewritten to the anchor pandoc will
# generate for that chapter's heading, so the link jumps to the chapter instead
# of failing.
#
# The slug rule is pandoc's own: lowercase, drop everything that is not
# alphanumeric, space or hyphen, then spaces to hyphens. Implemented here rather
# than guessed at, because a wrong slug is a link that silently goes nowhere.
slug <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9 -]", "", x)
  gsub(" +", "-", trimws(x))
}
ANCHOR <- stats::setNames(paste0("#", vapply(CHAPTERS, slug, character(1))),
                          paste0(names(CHAPTERS), ".html"))

relink <- function(lines) {
  for (html in names(ANCHOR)) {
    # Any fragment on the end is dropped: it names a section anchor from the
    # other page, which is not the anchor that section has here.
    pat <- sprintf("\\]\\(%s(#[^)]*)?\\)", gsub("[.]", "[.]", html))
    lines <- gsub(pat, sprintf("](%s)", ANCHOR[[html]]), lines)
  }
  lines
}

# RAW <img> BECOMES A MARKDOWN IMAGE. knitr emits an HTML <img> tag rather than
# markdown whenever a chunk sets out.width, which the worked-example article does
# on every one of its chunks. pandoc passes raw HTML through untouched and then
# DROPS it when the target is LaTeX, so that chapter came out with all 22 of its
# figures missing and no warning anywhere -- the text still read as though the
# figures were there.
#
# The alt attribute becomes the caption. These alt texts are full descriptive
# sentences written for screen readers, which is exactly what a figure caption
# should say, and without one a PDF figure is an unlabelled picture.
unhtml_img <- function(lines) {
  has <- grepl("<img ", lines, fixed = TRUE)
  if (!any(has)) return(lines)
  lines[has] <- vapply(lines[has], function(ln) {
    src <- sub('.*<img[^>]*src="([^"]*)".*', "\\1", ln)
    alt <- if (grepl('alt="', ln)) sub('.*alt="([^"]*)".*', "\\1", ln) else ""
    if (identical(src, ln)) return(ln)   # not an img after all, leave alone
    sprintf("![%s](%s)", alt, src)
  }, character(1), USE.NAMES = FALSE)
  lines
}

body <- character(0)
for (nm in names(CHAPTERS)) {
  rmd <- file.path(vig_dir, paste0(nm, ".Rmd"))
  if (!file.exists(rmd)) {
    message("  SKIP ", nm, ": no such vignette")
    next
  }
  md <- file.path(work, paste0(nm, ".md"))
  old <- setwd(vig_dir); on.exit(setwd(old), add = TRUE)
  knitr::knit(rmd, output = md, quiet = TRUE)
  setwd(old)
  txt <- relink(unhtml_img(demote(strip_yaml(readLines(md, warn = FALSE)))))
  # \clearpage rather than a page break inside the chapter: a chapter that starts
  # halfway down a page reads as a section of the one before it.
  body <- c(body, "", "\\clearpage", "",
            paste("#", CHAPTERS[[nm]]), "", txt, "")
  message("  chapter: ", CHAPTERS[[nm]], " (", length(txt), " lines)")
}

# Figures produced by the evaluated chunks, plus the static ones the articles
# include, both resolved relative to the assembled document.
for (d in c("figure", "figures")) {
  src <- file.path(vig_dir, d)
  if (dir.exists(src)) file.copy(src, work, recursive = TRUE)
}
if (dir.exists(file.path(work, "figure")))
  file.copy(file.path(work, "figure"), work, recursive = TRUE)

version <- as.character(read.dcf(file.path(pkg_root, "DESCRIPTION"))[1, "Version"])
header <- c(
  "---",
  'title: "cyRAVEN"',
  'subtitle: "Supervised immunophenotyping of multi-sample flow cytometry data"',
  sprintf('author: "Version %s"', version),
  sprintf('date: "%s"', format(Sys.Date(), "%d %B %Y")),
  "output:",
  "  pdf_document:",
  "    toc: true",
  "    toc_depth: 2",
  "    number_sections: true",
  "    latex_engine: xelatex",
  "    includes:",
  "      in_header: preamble.tex",
  "documentclass: report",
  "geometry: margin=2.5cm",
  "colorlinks: true",
  "linkcolor: RavenOrange",
  "urlcolor: RavenOrange",
  "toccolor: black",
  "---",
  "")

writeLines(c(header, body), file.path(work, "manual.md"))
file.copy(file.path(pkg_root, "pkgdown", "manual", "preamble.tex"), work,
          overwrite = TRUE)
file.copy(file.path(pkg_root, "pkgdown", "manual", "landscape-tables.lua"), work,
          overwrite = TRUE)
logo <- file.path(pkg_root, "man", "figures", "logo.png")
if (file.exists(logo)) file.copy(logo, work, overwrite = TRUE)

old <- setwd(work); on.exit(setwd(old), add = TRUE)
rmarkdown::render(
  "manual.md",
  output_file = "cyRAVEN-manual.pdf",
  # The Lua filter is what puts wide tables on their own landscape page; see
  # landscape-tables.lua for why it is a filter rather than a LaTeX package
  # option.
  output_options = list(pandoc_args = c("--lua-filter=landscape-tables.lua")),
  quiet = TRUE)
setwd(old)

file.copy(file.path(work, "cyRAVEN-manual.pdf"),
          file.path(out_dir, "cyRAVEN-manual.pdf"), overwrite = TRUE)
size_mb <- round(file.size(file.path(out_dir, "cyRAVEN-manual.pdf")) / 1024^2, 1)
message("wrote docs/cyRAVEN-manual.pdf (", size_mb, " MB, ",
        length(CHAPTERS), " chapters)")
