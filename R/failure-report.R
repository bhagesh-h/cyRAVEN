# SECTION 13b -- FAILURE DIAGNOSIS
# =============================================================================
#
# WHY THIS FILE EXISTS. A run that fails used to leave a results directory
# holding whatever it had written before it stopped, a manifest marked "failed",
# and an R error on stderr that the person opening the directory afterwards may
# never have seen. Reconstructing what happened then means reading a stack trace
# out of a container log, which is the point at which most users conclude the
# tool is broken rather than the input.
#
# So a failed run writes the same report a successful one writes, with a
# diagnosis at the top: what stopped it, where, the log leading up to it, what
# the error means in the vocabulary of the analysis rather than of R, and the
# specific next action. Everything the run did manage to produce is embedded
# below it, because the partial output is usually where the evidence is: a
# gating_qc.png written before the failure often shows the cause directly.
#
# WHY THE INTERPRETATIONS ARE A TABLE RATHER THAN GENERATED. Each entry pairs a
# message pattern with what that message means for the DATA. "subscript out of
# bounds" is not a useful thing to tell a cytometrist; "the specification names a
# marker this panel does not contain" is the same fact stated in terms they can
# act on. That mapping is knowledge about this pipeline and has to be written
# down; nothing can infer it from the condition object.

#' Known failure modes, and what each one means for the input
#'
#' Each entry is a regular expression matched against the error message, with
#' the meaning and the action to take. Order matters: the first match wins, so
#' specific patterns precede general ones.
#' @return list of list(pattern, meaning, action)
#' @keywords internal
failure_catalogue <- function() {
  list(
    list(pattern = "conflicting values for the same subject",
         meaning = paste("Two rows of the sample sheet belonging to one subject",
                         "disagree about one of that subject's attributes. A",
                         "subject attribute repeats on every row of that subject,",
                         "so the rows have to agree; there is no defensible way",
                         "to pick one."),
         action = paste("The error names every conflict as",
                        "'subject / column: value vs value'. Correct them in the",
                        "sheet, then re-check with --check before re-running.")),
    list(pattern = "not in the sample (sheet|map)",
         meaning = paste("The run found FCS files that the sheet has no row for.",
                         "Identifiers are not guessed from filename order,",
                         "because a wrong guess silently mislabels a patient."),
         action = paste("Add a row for each file named in the error, or",
                        "regenerate the sheet with --write-samples, which covers",
                        "every file in the input directory. Note that the sheet",
                        "matches on basename, so a file moved between",
                        "subdirectories still matches.")),
    list(pattern = "duplicate file entries",
         meaning = paste("Two rows of the sheet name the same FCS file, so which",
                         "row describes it is undefined."),
         action = "Remove the duplicate row named in the error."),
    list(pattern = "must contain a 'file' column",
         meaning = paste("The sheet has no column naming the FCS files, which is",
                         "the one column it cannot do without."),
         action = paste("Add a 'file' column holding each acquisition's filename,",
                        "or start from --write-samples.")),
    list(pattern = "--samples supersedes",
         meaning = paste("The run was given both the unified sheet and one of the",
                         "separate tables it replaces, so one fact would have two",
                         "sources."),
         action = paste("Keep --samples and drop --sample-map, --patient-table",
                        "and --absolute-counts; the sheet carries all three.")),
    list(pattern = "none of the requested columns were found",
         meaning = paste("The patient table was read but none of its headers",
                         "matched a known column name in either language."),
         action = paste("Check the header row. Add the spellings in use to",
                        "metadata: column_map: in the config; the canonical",
                        "names are patient_id, cohort, sex, age_years,",
                        "date_of_birth, height_cm, weight_kg, wbc_per_ul.")),
    list(pattern = "--reference-date must be",
         meaning = "The reference date is not in YYYY-MM-DD form.",
         action = "Pass it as e.g. --reference-date 2025-01-01."),
    # Observed on a liposome-uptake panel whose specification named one channel.
    # Every table was written and the run then stopped at the embedding with a
    # stopifnot() from run_umap(), naming a matrix the user never supplied and
    # giving no hint of the cause.
    list(pattern = "ncol\\(mat\\) >= 2",
         meaning = paste("The shared UMAP needs at least two channels and was",
                         "given one. The transformed matrix holds only the",
                         "channels the specification names, so a specification",
                         "naming a single marker leaves the embedding with one",
                         "feature column. Everything before this point ran: the",
                         "thresholds, frequencies and uncertainty tables in the",
                         "results directory are complete and usable."),
         action = paste("Name more channels in the config. A population does not",
                        "have to be added: markers listed under",
                        "functional_blocks: are transformed too, so listing the",
                        "rest of the panel there is enough, and it also reports",
                        "each one's intensity inside the populations already",
                        "declared.")),
    # Same panel, earlier failure. The specification used the channel name
    # --list-channels prints, which is not the name the scorer matches when it
    # holds anything but letters, digits and dots.
    list(pattern = "marker\\(s\\) not in panel|no any_of marker available",
         meaning = paste("A population names a channel the scored matrix does not",
                         "contain. Where a channel name holds a hyphen, a slash or",
                         "a space, the reader makes it syntactically valid the way",
                         "R does -- PE-A becomes PE.A, 'Blue B 710/50-A' becomes",
                         "Blue.B.710.50.A -- while --list-channels prints the",
                         "original. A config written from what that flag printed",
                         "matches nothing, and no error is raised: every",
                         "population is simply reported UNAVAILABLE."),
         action = paste("Replace every hyphen, slash and space in the channel",
                        "names in the config with a dot, then re-run --check,",
                        "which reports the specification's markers as present or",
                        "absent before any file is analysed. Panels whose markers",
                        "are named CD3, CD4, CD14 are unaffected.")),
    list(pattern = "cannot open|No such file|does not exist|cannot find",
         meaning = paste("A path the run was given does not exist as seen from",
                         "inside the process. Under Docker this is almost always",
                         "a host path passed where a container path was needed:",
                         "every path in a flag is resolved inside the container,",
                         "not on the host."),
         action = paste("Check that each -v mount covers the path, and that the",
                        "flag names the container side of it. On Git Bash, set",
                        "MSYS_NO_PATHCONV=1, which stops the shell rewriting",
                        "/data into a Windows path.")),
    list(pattern = "cannot allocate|memory exhausted|vector of size",
         meaning = paste("The run exceeded available memory. Peak memory is set",
                         "by the number of events held at once, not by the number",
                         "of files."),
         action = paste("Lower --max-events-per-file (300000 is ample for gate",
                        "derivation), then --cells-per-sample and --max-cells,",
                        "which bound the embedding. Under Docker Desktop, also",
                        "raise the VM memory limit.")),
    list(pattern = "there is no package called|could not find function",
         meaning = paste("An optional component was requested but the package",
                         "providing it is not installed. FlowSOM backs --cluster",
                         "and readxl backs .xlsx inputs; both are Suggests, so",
                         "the pipeline runs without them until a flag needs one."),
         action = paste("Install the named package, or drop the flag that needs",
                        "it. The container image carries both.")),
    list(pattern = "no marker|0 markers|UNAVAILABLE",
         meaning = paste("Marker names did not resolve. cyRAVEN reads marker",
                         "symbols from the $PnS keyword, and a specification name",
                         "that does not match $PnS exactly matches nothing. This",
                         "is the most common cause of an empty frequency table."),
         action = paste("Run with --check, which lists every marker resolved from",
                        "the files and names the specification entries that match",
                        "none of them.")),
    list(pattern = "subscript out of bounds|undefined columns",
         meaning = paste("Something was addressed by a name the data does not",
                         "carry, which usually traces back to a marker or column",
                         "name in the config that this cohort does not have."),
         action = paste("Run with --check to compare the specification against",
                        "the panel, and confirm the config's column names against",
                        "the sheet's header row.")))
}

#' Interpret an error message against the catalogue
#' @param msg the error message.
#' @return list(meaning, action) or NULL when no entry matches.
#' @keywords internal
diagnose_failure <- function(msg) {
  for (e in failure_catalogue())
    if (grepl(e$pattern, msg, ignore.case = TRUE, perl = TRUE))
      return(list(meaning = e$meaning, action = e$action))
  NULL
}

#' Assemble the diagnosis block placed at the top of a failed run's report
#'
#' @param err a condition, or a character message.
#' @param outdir the results directory, used to say what survived.
#' @param opt parsed options.
#' @return a character string of HTML.
#' @keywords internal
failure_block <- function(err, outdir, opt = NULL) {
  msg <- if (inherits(err, "condition")) conditionMessage(err) else as.character(err)
  cl  <- if (inherits(err, "condition") && !is.null(conditionCall(err)))
    paste(deparse(conditionCall(err)), collapse = " ") else NA_character_
  stage <- .cyraven_log$stage
  d <- diagnose_failure(msg)
  logl <- .cyraven_log$lines
  tail_n <- min(length(logl), 60L)
  logtail <- if (tail_n) utils::tail(logl, tail_n) else "(no log recorded)"

  n_fig <- length(list.files(outdir, "[.]png$"))
  n_tab <- length(list.files(outdir, "[.]csv$"))

  paste0(
    "<div class='banner stop'><b>This run did not complete.</b> ",
    "Nothing below it is a result: the outputs that survive are those written ",
    "before the failure, and every stage after it is missing. They are embedded ",
    "anyway, because the evidence for what went wrong is usually in them.</div>",

    "<section class='sec' id='fail'><details open><summary>",
    "<span class='chev' aria-hidden='true'></span>",
    "<span class='sec-title'>0. Why this run failed</span>",
    "<span class='sec-count'>diagnosis</span></summary><div class='sec-body'>",

    "<h3>What stopped it</h3>",
    "<pre class='err'>", html_escape(msg), "</pre>",
    if (!is.na(cl))
      paste0("<p class='none'>Raised by: <span class='mono'>",
             html_escape(substr(cl, 1, 400)), "</span></p>") else "",

    "<h3>Where</h3><p>",
    if (is.na(stage))
      "Before the first stage, while reading options and inputs."
    else paste0("During <b>", html_escape(stage), "</b>. Every stage before it ",
                "completed; every stage after it did not run."),
    " ", n_fig, " figure(s) and ", n_tab,
    " table(s) had been written at that point.</p>",

    "<h3>What it means</h3><p>",
    if (!is.null(d)) html_escape(d$meaning)
    else paste0("This message is not one of the failure modes cyRAVEN recognises, ",
                "so it carries no interpretation beyond the text above. The log ",
                "below shows what the run was doing when it stopped."),
    "</p>",

    "<h3>What to do</h3><p>",
    if (!is.null(d)) html_escape(d$action)
    else paste0("Read the log below from the bottom up: the last line before the ",
                "failure names the file or the stage in question. Re-run with ",
                "<span class='mono'>--check</span>, which validates the inputs ",
                "from the FCS headers alone in seconds and reports marker, sheet ",
                "and group-column problems without analysing anything. If the ",
                "cause is still not evident, the log, run_manifest.txt and this ",
                "file together are what an issue report needs."),
    "</p>",

    "<h3>Evidence: the run log up to the failure</h3>",
    "<p class='none'>Last ", tail_n, " of ", length(logl), " line(s).</p>",
    "<pre class='log'>", html_escape(paste(logtail, collapse = "\n")), "</pre>",

    "<h3>Before re-running</h3><p>Validate the inputs without analysing them:</p>",
    "<pre class='cmd'>", html_escape(paste0(
      "docker run --rm -v \"$PWD/data:/data:ro\" -v \"$PWD/results:/results\" ",
      "cyraven:", tryCatch(as.character(utils::packageVersion("cyRAVEN")),
                           error = function(e) "latest"),
      " \\\n  --dir /data --samples /data/samples.csv --config /data/analysis.yaml ",
      "\\\n  --outdir /results --check")),
    "</pre>",
    "</div></details></section>")
}
