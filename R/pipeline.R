#' Fill an option list from the command-line defaults
#'
#' WHY THIS EXISTS. `run_cyraven()` is documented -- in the README, in the
#' vignette, in its own examples -- as taking a PARTIAL list: name the two or
#' three options that matter and let the rest default. That only worked when the
#' list arrived from `optparse`, which supplies every default itself. Called
#' directly with a short list, every unnamed option was NULL, and the failures
#' that produced were silent and remote from their cause: `singlet_mad_k`
#' missing made the singlet band `median +/- NULL * mad`, which selects no
#' events, so every marker threshold resolved over an empty vector to NA, and
#' the run reported that every population's markers were "not in panel" --
#' naming the panel for a fault in the caller's option list.
#'
#' Taking the defaults from `build_option_list()` rather than restating them
#' means there is one place a default is written down, and the programmatic and
#' command-line paths cannot drift apart.
#'
#' @param opt Named list, possibly partial.
#' @return `opt` with every unset option filled from its command-line default.
#' @keywords internal
fill_option_defaults <- function(opt) {
  ol <- build_option_list()
  dest <- vapply(ol, function(o) o@dest, character(1))
  defs <- lapply(ol, function(o) o@default)
  names(defs) <- dest
  defs <- defs[nzchar(dest) & !vapply(defs, is.null, logical(1))]
  # modifyList drops NULL entries on the right-hand side, which is what we want:
  # explicitly passing NULL should mean "unset", not "override the default".
  utils::modifyList(defs, opt %||% list())
}

#' run_cyraven
#'
#' Runs the pipeline. A run that stops on an error still writes `report.html`,
#' carrying the diagnosis and every output produced before the failure, and then
#' re-raises the original error unchanged so the exit status and the message a
#' caller sees are what they would have been.
#'
#' @param opt Named list of options; anything omitted takes its command-line
#'   default. See build_option_list() for the full set.
#' @keywords internal
#' @export
run_cyraven <- function(opt) {
  log_reset()
  # WHY THE REPORT IS WRITTEN HERE AND NOT IN on.exit(): an exit handler cannot
  # see the condition that ended the call, only that it ended, and "something
  # failed" is not a diagnosis. Catching it here gives the report the message,
  # the call and the stage; re-raising immediately afterwards leaves the caller's
  # view of the failure untouched.
  tryCatch(run_cyraven_impl(opt), error = function(e) {
    od <- tryCatch(fill_option_defaults(opt)$outdir, error = function(...) NULL)
    if (!is.null(od) && dir.exists(od) &&
        !isTRUE(fill_option_defaults(opt)[["check_only", exact = TRUE]]))
      try(write_run_report(od, opt = opt, failure = e), silent = TRUE)
    stop(e)
  })
}

#' The pipeline itself
#' @param opt Parsed options.
#' @keywords internal
run_cyraven_impl <- function(opt) {
  opt <- fill_option_defaults(opt)
  set.seed(opt$seed)
  fcs <- resolve_input_files(opt)
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  log_step("INPUT: ", length(fcs), " FCS file(s) -> ", opt$outdir)

  # ---- run manifest, written NOW and rewritten at the end -------------------
  # Written before the expensive steps on purpose: a run that dies at STEP 6 is
  # exactly when you want to know what it was invoked with, which files it had,
  # and which package versions were loaded. A manifest still reading
  # status: running is itself the diagnosis.
  .run_started <- Sys.time()
  .manifest <- file.path(opt$outdir, "run_manifest.txt")
  if (TRUE)
    try(write_run_manifest(.manifest, opt = opt, files = fcs,
                           status = "running", started = .run_started),
        silent = TRUE)
  # Rewrite on the way out however the run ends. on.exit fires on error and on
  # interrupt too, so a crashed run leaves "failed" rather than a stale
  # "running" that a later reader would mistake for a run still in progress.
  .finished_ok <- FALSE
  on.exit({
    if (TRUE)
      try(write_run_manifest(
            .manifest, opt = opt, files = fcs,
            extra = {
              .x <- list()
              if (exists("cofactors", inherits = FALSE))
                .x$cofactors <- unlist(cofactors)
              # Manual threshold corrections belong in the provenance record, not
              # only in the table they altered. A reader asking "was this run
              # touched by hand" should not have to grep a CSV to find out.
              if (exists("cfg_ovr", inherits = FALSE) && length(cfg_ovr))
                .x$manual_threshold_overrides <- unlist(lapply(names(cfg_ovr),
                  function(s) vapply(names(cfg_ovr[[s]]), function(m) {
                    e <- cfg_ovr[[s]][[m]]
                    if (!is.list(e)) e <- list(threshold = e)
                    paste0(s, " / ", m, " = ", e$threshold,
                           "  by=", e$set_by %||% "unattributed",
                           "  reason=", e$reason %||% "none given")
                  }, character(1))), use.names = FALSE)
              if (length(.x)) .x else NULL
            },
            status = if (.finished_ok) "completed" else "failed",
            started = .run_started), silent = TRUE)
  }, add = TRUE)

  if (!is.null(opt$write_sample_map)) {
    write_sample_map_template(fcs, opt$write_sample_map); return(invisible(NULL))
  }
  if (!is.null(opt[["write_samples", exact = TRUE]])) {
    write_samplesheet_template(fcs, opt$write_samples); return(invisible(NULL))
  }
  # Before the config is parsed, deliberately: --list-channels answers "what do I
  # write in the config", so requiring one to run it would be circular.
  if (isTRUE(opt$list_channels)) {
    list_channels(fcs); return(invisible(NULL))
  }

  cfg <- if (!is.null(opt$config)) yaml::read_yaml(opt$config) else list()
  spec   <- cfg$populations       %||% default_population_spec()
  blocks <- cfg$functional_blocks %||% default_functional_blocks()
  # Unlike `blocks`, no built-in default: a ratio hardcoded to population
  # names from one study's panel would be silently wrong (or absent) on
  # another, so this is opt-in per config (see compute_population_ratios()).
  ratios <- cfg$ratios %||% list()
  cfg_thr <- cfg$thresholds       %||% list()
  # Per-sample, per-marker manual corrections. Distinct from cfg_thr, which
  # applies one number to every sample. See resolve_threshold() for why the two
  # are not the same claim. `.ovr_applied` accumulates the keys that actually
  # matched, so an entry naming a sample this cohort does not contain is
  # reported rather than silently ignored.
  cfg_ovr <- cfg$sample_overrides %||% list()
  .ovr_applied <- character(0)
  # Fluorescence-minus-one controls, declared through the sample map rather than
  # the config because which tube is a control is a property of the acquisition.
  # Both are absent on a run that declares none, and nothing below changes.
  .fmo_map <- NULL; .fmo_group <- NULL; .fmo_used <- list()
  # Every fig_*() call below defaults to `colors = fcs_colors()`, resolved as a
  # default ARGUMENT at call time -- setting the package's palette here is what
  # makes a --config `colors:` override reach every figure with no other change.
  # See default_colors() for the full key list.
  set_fcs_colors(cfg$colors)
  # Anchors the study colours (green = reference) for every discrete scale.
  set_reference_group(opt$reference_group)
  min_cd45 <- cfg$gating$staining_qc$min_cd45_pct_of_live %||% opt$min_cd45_pct
  mad_k    <- cfg$gating$singlet_band$mad_k %||% opt$singlet_mad_k

  # ---- inputs: one sheet, or the three separate tables ----------------------
  # `.sheet` is non-NULL only on the unified route. Its three components are the
  # same structures the separate loaders return, so every consumer below reads
  # one of them without knowing which route supplied it.
  .sheet <- NULL
  if (!is.null(opt[["samples", exact = TRUE]])) {
    clash <- c("--sample-map"      = !is.null(opt$sample_map),
               "--patient-table"   = !is.null(opt$patient_table),
               "--absolute-counts" = !is.null(opt[["absolute_counts", exact = TRUE]]))
    if (any(clash))
      stop("--samples supersedes ", paste(names(clash)[clash], collapse = ", "),
           " and cannot be combined with ", if (sum(clash) > 1L) "them" else "it",
           ".\n  The sheet carries the same facts in one file; two sources of ",
           "truth for one fact is what this format removes.", call. = FALSE)
    log_step("STEP 0b - reading the sample sheet")
    .refd_sheet <- if (!is.null(opt$reference_date)) {
      rr <- as.Date(opt$reference_date, format = "%Y-%m-%d")
      if (is.na(rr)) stop("--reference-date must be YYYY-MM-DD, got: ",
                          opt$reference_date, call. = FALSE)
      rr
    } else Sys.Date()
    .sheet <- read_samplesheet(
      opt$samples, fcs,
      column_map = cfg$metadata$column_map %||% default_column_map(),
      value_map  = cfg$metadata$value_translations %||% default_value_map(),
      reference_date = .refd_sheet,
      count_unit = cfg$samples$count_unit %||% "cells/uL")
    # The YAML may name which study column plays which role, so that the pair of
    # files specifies a run without further flags. An explicit flag still wins:
    # a config is a study's standing choice, a flag is this run's.
    for (r in c("group_column", "batch_column")) {
      v <- cfg$samples[[r]]
      if (!is.null(v) && is.null(opt[[r]])) {
        opt[[r]] <- v
        log_msg("  ", sub("_", "-", paste0("--", r)), " ", v, " (from the config)")
      }
    }
    # clinical_columns is a LIST in the YAML and a comma-separated string on the
    # command line, so it is collapsed to one form here rather than handled
    # twice downstream.
    if (is.null(opt[["clinical_columns", exact = TRUE]]) &&
        !is.null(cfg$samples[["clinical_columns"]])) {
      opt$clinical_columns <- paste(unlist(cfg$samples[["clinical_columns"]]),
                                    collapse = ",")
      log_msg("  --clinical-columns ", opt$clinical_columns, " (from the config)")
    }
  }

  smap <- if (!is.null(.sheet)) .sheet$smap
          else if (!is.null(opt$sample_map)) load_sample_map(opt$sample_map, fcs)
          else NULL
  if (!is.null(smap)) {
    .fmo_map <- parse_fmo_map(smap)
    if (!is.null(.fmo_map)) {
      .sid <- smap$sample_id %||% smap$well %||% smap$file
      .fmo_group <- if ("control_group" %in% names(smap))
        setNames(trimws(as.character(smap$control_group)), as.character(.sid)) else NULL
      log_msg("fluorescence-minus-one control(s) declared for ",
              length(unique(.fmo_map$marker)), " marker(s) across ",
              length(unique(.fmo_map$sample_id)), " file(s): ",
              paste(utils::head(sort(unique(.fmo_map$marker)), 8), collapse = ", "),
              if (length(unique(.fmo_map$marker)) > 8) ", ..." else "")
    }
  }

  # ---- validate and stop, without analysing ---------------------------------
  # A run over a large cohort costs many minutes, and the errors this catches
  # (a filename the sheet does not cover, a marker name the specification does
  # not match, a subject whose rows disagree) are all knowable from the headers
  # alone. Reading only headers makes the check cost seconds.
  if (isTRUE(opt[["check_only", exact = TRUE]])) {
    log_step("CHECK - validating inputs, no analysis")
    report_input_check(fcs, smap, .sheet, spec, opt)
    return(invisible(NULL))
  }

  # ---- read -----------------------------------------------------------------
  log_step("STEP 1 - reading files and resolving markers")
  ignore_ch <- if (!is.null(opt$ignore_channels) && nzchar(opt$ignore_channels))
    trimws(strsplit(opt$ignore_channels, ",")[[1]]) else NULL
  if (length(ignore_ch))
    log_msg("ignoring channel(s) before panel fingerprinting: ",
            paste(ignore_ch, collapse = ", "))

  .sid_for <- function(f) {
    sid <- if (!is.null(smap)) {
      row <- smap[smap$file == basename(f), ]
      (row$sample_id %||% row$well %||% NA)[1]
    } else NA
    if (is.na(sid) || !nzchar(as.character(sid))) NULL else sid
  }
  .read_one <- function(f) {
    rd <- read_fcs_resolved(f, sample_id = .sid_for(f),
                            max_events = opt$max_events_per_file %||% 0L,
                            ignore_channels = ignore_ch)
    rd$exprs <- maybe_compensate(rd$exprs, rd$keywords)
    rd
  }

  # ---- reading, optionally in parallel --------------------------------------
  # Files are independent of one another, and reading is the slowest stage on a
  # large cohort, so this is the one place parallelism pays. Everything after it
  # stays sequential: cofactors, thresholds and the embedding are derived across
  # samples, and the seeded RNG that places them has to be consumed in one
  # order or a run stops reproducing itself.
  #
  # WHY FORKING AND NOT A SOCKET CLUSTER. Workers need the FCS data and give
  # back a matrix each. A fork shares memory copy-on-write and returns the
  # result directly; a socket cluster would serialise every matrix twice. Fork
  # is Unix-only, which is what the container is, and Windows falls back to
  # sequential reading with a note rather than silently ignoring the flag.
  n_threads <- suppressWarnings(as.integer(opt$read_threads %||% 1L))
  if (is.na(n_threads) || n_threads < 1L) n_threads <- 1L
  n_threads <- min(n_threads, length(fcs))
  if (n_threads > 1L && .Platform$OS.type != "unix") {
    log_msg("--read-threads ", n_threads, " ignored: forking needs a Unix host, ",
            "reading sequentially")
    n_threads <- 1L
  }

  if (n_threads > 1L) {
    log_msg("reading ", length(fcs), " file(s) on ", n_threads, " thread(s). ",
            "Peak memory scales with the thread count.")
    reads <- parallel::mclapply(fcs, .read_one, mc.cores = n_threads,
                                mc.preschedule = FALSE)
    bad <- vapply(reads, function(r) inherits(r, "try-error") || is.null(r$sample_id),
                  logical(1))
    if (any(bad))
      stop("parallel read failed for ", sum(bad), " file(s): ",
           paste(basename(fcs[bad]), collapse = ", "),
           ". Re-run with --threads 1 for the actual error.")
    names(reads) <- vapply(reads, `[[`, "", "sample_id")
    for (rd in reads)
      log_msg("  ", basename(rd$file), ": ", rd$n_events, " events",
              if (isTRUE(rd$subsampled))
                paste0(" (evenly sampled from ", rd$n_events_file, " in file)") else "",
              ", ", length(rd$marker_cols),
              " markers, sample_id '", rd$sample_id, "'")
  } else {
    reads <- list()
    for (f in fcs) {
      log_msg("reading ", basename(f), " (", round(file.size(f) / 1e6), " MB)")
      rd <- .read_one(f)
      log_msg("  ", rd$n_events, " events",
              if (isTRUE(rd$subsampled))
                paste0(" (evenly sampled from ", rd$n_events_file, " in file)") else "",
              ", ", length(rd$marker_cols),
              " markers, sample_id '", rd$sample_id, "'")
      reads[[rd$sample_id]] <- rd
    }
  }
  fpr <- fingerprint_panels(reads)

  # ---- bead calibration ------------------------------------------------------
  # Applied here, before any cofactor or threshold is derived, because the point
  # is to change the units everything downstream is expressed in. Applying it
  # later would leave the thresholds on one scale and the medians on another.
  if (!is.null(opt$calibration_beads)) {
    .cal <- tryCatch({
      if (is.null(opt$calibration_values))
        stop("--calibration-beads needs --calibration-values", call. = FALSE)
      bd <- read_fcs_resolved(opt$calibration_beads)
      av <- read.csv(opt$calibration_values, stringsAsFactors = FALSE,
                     check.names = FALSE)
      fit_bead_calibration(bd$exprs, bd$marker_cols, av,
                           min_r2 = opt$calibration_min_r2 %||% 0.98)
    }, error = function(e) {
      log_msg("WARNING bead calibration failed: ", conditionMessage(e),
              ". The run continues in instrument units")
      NULL
    })
    if (!is.null(.cal)) {
      write.csv(.cal, file.path(opt$outdir, "calibration.csv"), row.names = FALSE)
      .nok <- sum(grepl("^calibrated", .cal$verdict))
      log_msg("wrote calibration.csv (", .nok, " of ", nrow(.cal),
              " channel(s) converted to the assigned units)")
      if (.nok < nrow(.cal))
        log_msg("  NOTE ", nrow(.cal) - .nok, " channel(s) were NOT converted ",
                "and remain in instrument units. Intensities from calibrated ",
                "and uncalibrated channels are not on the same scale and must ",
                "not be compared. See the verdict column")
      if (.nok) {
        for (s in names(reads)) {
          ap <- apply_bead_calibration(reads[[s]]$exprs, reads[[s]]$marker_cols, .cal)
          reads[[s]]$exprs <- ap$exprs
        }
        log_msg("  applied to ", length(reads), " sample(s); every intensity ",
                "below is in the assigned units, not channel units")
      }
    }
  }

  # ---- acquisition-time stability -------------------------------------------
  # Placed here, before any threshold is derived, because this is the last point
  # at which excluding events is still a coherent operation: every cofactor,
  # threshold and frequency below is computed from whatever events survive this
  # step. Reporting is unconditional; removal is not.
  .aqc <- NULL
  if (!isTRUE(opt$no_acquisition_qc)) {
    log_step("STEP 1b - acquisition-time stability")
    .aqc <- tryCatch(
      run_acquisition_qc(reads, n_bins = opt$acquisition_bins %||% 40L,
                         mad_k = opt$acquisition_mad_k %||% 5),
      error = function(e) {
        log_msg("  WARNING acquisition QC failed: ", conditionMessage(e),
                ", every other output is unaffected")
        NULL
      })
    if (!is.null(.aqc)) {
      write.csv(.aqc$summary, file.path(opt$outdir, "acquisition_qc.csv"),
                row.names = FALSE)
      n_un <- sum(.aqc$summary$verdict %in% c("unstable", "unstable, minor"))
      log_msg("wrote acquisition_qc.csv (", n_un, " of ", nrow(.aqc$summary),
              " sample(s) carry a flagged interval)")
      if (!is.null(.aqc$bins))
        write.csv(.aqc$bins, file.path(opt$outdir, "acquisition_qc_bins.csv"),
                  row.names = FALSE)
      # The FIGURE is deliberately not drawn here. Drawing it at this point
      # changes the rendering of nine of the figures drawn later: the data
      # behind them is provably unaffected -- every table is byte-identical
      # either way -- but the pixels are not, through session state in the
      # ggplot2 and grid stack that could not be reproduced outside a full run.
      # A new diagnostic that silently re-renders every existing figure is
      # exactly the kind of change this package treats as unacceptable, so the
      # plot is drawn with the other figures instead, after them. See STEP 7c.
      if (any(.aqc$summary$verdict == "no Time channel"))
        log_msg("  NOTE ", sum(.aqc$summary$verdict == "no Time channel"),
                " file(s) carry no Time channel, so acquisition stability ",
                "cannot be assessed for them")

      # Removal, when asked for. Done by rewriting the event matrix, so every
      # stage below sees the cleaned file and nothing needs to know this
      # happened. The counts change, which is exactly why it is opt-in.
      if (isTRUE(opt$drop_unstable_events)) {
        n_drop <- 0L; n_files <- 0L
        for (s in names(.aqc$flagged)) {
          fl <- .aqc$flagged[[s]]
          if (is.null(fl) || !any(fl) || length(fl) != nrow(reads[[s]]$exprs)) next
          reads[[s]]$exprs <- reads[[s]]$exprs[!fl, , drop = FALSE]
          reads[[s]]$n_events <- nrow(reads[[s]]$exprs)
          n_drop <- n_drop + sum(fl); n_files <- n_files + 1L
        }
        log_msg("--drop-unstable-events removed ", n_drop, " event(s) from ",
                n_files, " file(s). Every count, threshold and frequency below ",
                "is computed on the remaining events, and the decision is ",
                "recorded in the run manifest")
      }
    }
  }

  # ---- explore-only ---------------------------------------------------------
  # Standalone unsupervised discovery: no specification, no gate hierarchy, no
  # declared populations. Only the transform is needed, and it is derived here
  # rather than by borrowing STEP 2, so this path cannot fail for a reason that
  # belongs to the declared analysis. Everything else below is skipped.
  if (isTRUE(opt$explore_only)) {
    log_step("EXPLORE ONLY - unsupervised discovery, declared analysis skipped")
    tr_only <- list()
    for (p in fpr$panels) {
      cf <- opt[["cofactor", exact = TRUE]] %||% derive_cofactor_pooled(reads, p$samples)
      tr_only[[p$name]] <- make_transform("arcsinh", cofactor = cf)
    }
    run_explore(reads, fpr, opt, opt$outdir, transforms = tr_only,
                file_paths = fcs)
    # No manifest write here: the on.exit handler above already owns that file,
    # writes it to the right path, and flips it to "completed" off .finished_ok.
    # An explicit call duplicated it and passed the DIRECTORY where the function
    # wants the file path, which failed with "cannot open file: Is a directory".
    log_step("DONE - outputs in ", normalizePath(opt$outdir))
    .finished_ok <- TRUE
    return(invisible(TRUE))
  }

  # ---- gate -----------------------------------------------------------------
  log_step("STEP 2 - deriving transform and gates")
  cofactors <- list(); gates <- list(); verdicts <- list(); recon <- list()
  transforms <- list()
  for (p in fpr$panels) {
    sids <- p$samples
    # Cofactor sources, in priority order: an explicit --cofactor, a value
    # recorded in the config, then derivation from the data. Derivation now pools
    # across samples by default (derive_cofactor_pooled) instead of trusting
    # whichever file sorted first -- see that function for why, and pass
    # --cofactor-from-first-sample to restore the previous behaviour exactly.
    # opt[["cofactor", exact = TRUE]], never opt$cofactor -- see the dest comment
    # on --cofactor-from-first-sample for what partial matching did here.
    .cf_opt <- opt[["cofactor", exact = TRUE]]
    cf <- .cf_opt %||% cfg$transform$cofactor$derived_from_batch[[p$name]] %||%
          (if (isTRUE(opt$first_sample_cofactor))
             derive_cofactor(reads[[sids[1]]]$exprs, reads[[sids[1]]]$marker_cols)
           else derive_cofactor_pooled(reads, sids))
    # A cofactor of zero or less is not a setting, it is a bug upstream: the
    # transform is asinh(x / cofactor), so it yields Inf for every marker and
    # every downstream threshold. Fail here, naming the value, rather than let it
    # propagate into 25 identical "no separable CD45+ mode" verdicts.
    if (!is.finite(cf) || cf <= 0)
      stop("cofactor for panel '", p$name, "' resolved to ", format(cf),
           ", which cannot be used: the arcsinh transform divides by it. ",
           "Pass a positive --cofactor, or omit it to derive one from the data.")
    cofactors[[p$name]] <- cf

    # The transform every later step applies. Built once per panel and threaded
    # through, so gating, scoring, embedding and the figures all use the same
    # one and none of them has to know which it is. Logicle parameters are
    # pooled over the panel for the reason given in transform-methods.R: a
    # per-file fit puts every sample on its own ruler and makes the
    # cross-sample comparisons this package exists for meaningless.
    tmethod <- opt$transform %||% "arcsinh"
    transforms[[p$name]] <- if (identical(tmethod, "logicle")) {
      lp <- derive_logicle_pooled(reads, sids, m = opt$logicle_m %||% 4.5,
                                  seed = opt$seed %||% 42L)
      if (is.null(lp))
        stop("panel '", p$name, "': no logicle parameters could be derived")
      log_msg(sprintf("%s: logicle transform, %d marker(s), pooled over %d sample(s)",
                      p$name, length(lp), attr(lp, "n_samples_pooled") %||% NA))
      make_transform("logicle", logicle = lp)
    } else {
      make_transform(tmethod, cofactor = cf)
    }
    log_msg(p$name, ": cofactor ", round(cf, 1),
            if (!is.null(.cf_opt)) " (supplied)"
            else if (isTRUE(opt$first_sample_cofactor))
              paste0(" (derived from ", sids[1], " only)")
            else paste0(" (median of ", attr(cf, "n_samples") %||% 1L,
                        " per-sample derivations)"))
    for (s in sids) {
      log_msg("gating ", s)
      rd <- reads[[s]]
      declared <- if (!is.null(smap) && "is_control" %in% names(smap))
        isTRUE(smap$is_control[smap$file == rd$file][1]) else NA
      g <- apply_gate_hierarchy(rd, cf, cfg_thr, control_ref = NULL,
                                singlet_k = mad_k,
                                viability_name = opt$viability_marker,
                                transform = transforms[[p$name]],
                                overrides = cfg_ovr[[s]])
      v <- staining_verdict(g, declared, min_cd45 * 1,
                            force_include = isTRUE(opt$include_qc_failed))
      log_msg("  ", v$verdict)
      gates[[s]] <- g; verdicts[[s]] <- v
    }
  }

  # unstained controls become the reference for unimodal markers, then re-gate
  ctrl_by_panel <- list()
  n_failed <- sum(vapply(verdicts, function(v)
    identical(v$qc_status, "failed"), logical(1)))
  if (n_failed)
    log_msg(n_failed, " sample(s) failed staining QC and were EXCLUDED but not ",
            "used as reference (declared samples, not controls)")
  if (isTRUE(opt$include_qc_failed)) {
    n_forced <- sum(vapply(verdicts, function(v)
      grepl("included anyway", v$verdict %||% "", fixed = TRUE), logical(1)))
    if (n_forced)
      log_msg("NOTE --include-qc-failed forced ", n_forced, " sample(s) that would ",
              "otherwise have failed staining QC into gating/UMAP/frequencies, see ",
              "staining_qc.csv for which, and treat their numbers as unreliable")
  }
  # WHY A DISCOVERED CONTROL IS NOT BELIEVED WHEN NONE WAS DECLARED.
  #
  # staining_verdict() reads an ABSENT is_control column as "might be a
  # control", because only `is_control = FALSE` is a positive assertion that a
  # file is a biological sample. On a sheet with no such column, any sample
  # whose CD45 gate found no density minimum was therefore promoted to the
  # unstained reference for its whole panel.
  #
  # The consequence is severe and silent. Every threshold in that panel then
  # comes from the 99.5th percentile of that one sample. If it is in fact a
  # stained biological sample rather than an unstained tube, every cut lands at
  # the top of a real distribution, and every population collapses to a
  # fraction of a percent. Observed on a real cohort: T cells at 0.034% of
  # CD45+ where the correct answer was over a hundred times that. No error, no
  # warning, just a table of near-zero frequencies that looks like a finding.
  #
  # So a control is believed only when the sheet says so. Absence of evidence
  # is not evidence of an unstained tube, and a cohort that contains no
  # controls is the ordinary case rather than the exception. --discover-controls
  # restores the old behaviour for anyone who really does ship unstained tubes
  # without labelling them.
  declared_control <- !is.null(smap) && "is_control" %in% names(smap) &&
    any(vapply(smap$is_control, isTRUE, logical(1)))
  discovered <- names(verdicts)[vapply(verdicts, function(v)
    isTRUE(v$is_reference), logical(1))]
  trust_discovered <- declared_control || isTRUE(opt$discover_controls)

  if (length(discovered) && !trust_discovered) {
    log_msg("NOTE ", length(discovered), " sample(s) failed staining QC and would ",
            "have been used as the unstained reference for every threshold: ",
            paste(discovered, collapse = ", "))
    log_msg("  the sheet declares no control, so they are treated as FAILED ",
            "samples instead and each threshold is derived from its own sample.")
    log_msg("  if one of them IS an unstained tube, declare it with ",
            "is_control=TRUE in the sheet, or pass --discover-controls.")
    for (s in discovered) {
      verdicts[[s]]$is_reference <- FALSE
      verdicts[[s]]$is_control   <- FALSE
      verdicts[[s]]$qc_status    <- if (isTRUE(opt$include_qc_failed)) "pass" else "failed"
      verdicts[[s]]$include      <- isTRUE(opt$include_qc_failed)
      verdicts[[s]]$verdict      <- paste0(
        verdicts[[s]]$verdict,
        " [no control declared in the sheet, so NOT used as a reference]")
    }
    discovered <- character(0)
  }

  for (s in discovered) {
    pn <- fpr$assignment[[s]]
    if (is.null(ctrl_by_panel[[pn]])) ctrl_by_panel[[pn]] <- s
  }
  if (length(ctrl_by_panel)) {
    log_msg("control reference per panel: ",
            paste(names(ctrl_by_panel), unlist(ctrl_by_panel), sep = "=", collapse = ", "))
    log_msg("  every threshold in these panel(s) is the 99.5th percentile of that ",
            "sample. Check thresholds_used.csv: a 'source' of control_q995 on ",
            "most markers, with frequencies far below what the panel should ",
            "yield, means the reference is not an unstained tube.")
  }

  # ---- QC figure FIRST ------------------------------------------------------
  log_step("STEP 3 - initial QC diagnostics")
  for (s in names(gates)) {
    g <- gates[[s]]; rd <- reads[[s]]
    lin <- intersect(c("CD3", "CD4", "CD8", "CD14", "CD19", "CD56", "HLA-DR"),
                     names(rd$marker_cols))
    .tr <- transforms[[fpr$assignment[[s]]]] %||%
             make_transform("arcsinh", cofactor = g$cofactor)
    dens <- lapply(lin, function(m)
      .tr$fn(rd$exprs[g$masks$cd45_pos, rd$marker_cols[[m]]], m))
    names(dens) <- lin
    recon[[s]] <- list(sample_id = s, panel = fpr$assignment[[s]],
                       verdict = verdicts[[s]]$verdict,
                       is_control = isTRUE(verdicts[[s]]$is_control),
                       qc_status = verdicts[[s]]$qc_status %||% "pass",
                       cd45_source = g$cd45_source %||% NA_character_,
                       fsc_log10 = g$scatter_gate$fsc_log10,
                       ssc_log10 = g$scatter_gate$ssc_log10,
                       cd45 = g$cd45_x,
                       gate = g$scatter_gate[c("fsc_lo","fsc_hi","ssc_lo","ssc_hi")],
                       cd45_threshold = g$cd45_threshold,
                       marker_densities = dens)
  }
  fig_recon_diagnostics(recon, file.path(opt$outdir, "recon_diagnostics.png"))

  # ---- populations ----------------------------------------------------------
  log_step("STEP 4 - scoring populations")
  pops <- list(); thr_rows <- list(); qc_recon <- list()
  for (s in names(gates)) {
    g <- gates[[s]]; rd <- reads[[s]]; cf <- g$cofactor
    pn <- fpr$assignment[[s]]
    ctrl_s <- ctrl_by_panel[[pn]]
    # Channels a population may reference: fluorescence markers, plus any
    # SCATTER channel the spec names.
    #
    # WHY scatter must be gateable: granulocytes are conventionally defined by
    # high side scatter (dense cytoplasmic granules), not by a lineage marker --
    # there is no antibody in the panel that identifies them. A DSL that can only
    # threshold fluorescence cannot express that population at all. Scatter is
    # log10-transformed rather than asinh-transformed because it is strictly
    # positive and spans decades, so the asinh cofactor derived for fluorescence
    # is meaningless for it; the threshold machinery only needs a scale on which
    # the modes separate.
    needed <- unique(unlist(lapply(spec, function(d)
      c(setdiff(names(d), "any_of"), names(d[["any_of"]])))))
    # Functional-block markers must be transformed too, even when NO gate
    # references them.
    #
    # WHY: tmat is the only transformed matrix downstream, and both MFI tables read
    # from it. A marker that is in the panel but appears in no population gate
    # (CD11c here -- a dendritic/myeloid marker reported per population but never
    # used to define one) is therefore absent from tmat, and its functional-block
    # rows vanish. The failure is silent: the table is written, it simply has no
    # rows for that marker, and the "absent from all panels" note blames the panel
    # for a marker the panel actually contains.
    needed <- unique(c(needed, unlist(lapply(blocks, `[[`, "markers"))))
    avail   <- intersect(needed, names(rd$marker_cols))
    sc_need <- intersect(needed, names(rd$scatter_cols))
    tr <- transforms[[pn]] %||% make_transform("arcsinh", cofactor = cf)
    tmat <- vapply(avail, function(m) tr$fn(rd$exprs[, rd$marker_cols[[m]]], m),
                   numeric(nrow(rd$exprs)))
    colnames(tmat) <- avail
    if (length(sc_need)) {
      smat <- vapply(sc_need, function(k)
        log10(pmax(rd$exprs[, rd$scatter_cols[[k]]], 1)), numeric(nrow(rd$exprs)))
      colnames(smat) <- sc_need
      tmat <- cbind(tmat, smat)
      avail <- c(avail, sc_need)
      log_msg("  scatter channel(s) used as gates: ", paste(sc_need, collapse = ", "),
              " (log10 scale)")
    }
    thr <- setNames(rep(NA_real_, length(avail)), avail)
    tdet <- list(); qdens <- list()
    for (m in avail) {
      cx <- if (!is.null(ctrl_s) && ctrl_s != s && m %in% names(reads[[ctrl_s]]$marker_cols))
        tr$fn(reads[[ctrl_s]]$exprs[gates[[ctrl_s]]$masks$single_cells,
                                    reads[[ctrl_s]]$marker_cols[[m]]], m) else NULL
      .ov <- sample_override(cfg_ovr, s, m)
      if (!is.null(.ov)) .ovr_applied <- c(.ovr_applied, paste0(s, "\r", m))
      # An FMO, where one is declared for this marker and applies to this
      # sample, is the better reference: it is the same panel with one reagent
      # left out, so its distribution in this channel IS the negative population
      # under the spreading the real samples experience. An unstained tube
      # cannot show that, and a cut anchored to it sits too low. The arithmetic
      # is the same either way; only the reference and the recorded source
      # differ. See fmo.R.
      .kind <- "control_q995"
      .fs <- fmo_for_sample(.fmo_map, s, m, .fmo_group)
      if (!is.na(.fs) && !is.null(reads[[.fs]]) &&
          m %in% names(reads[[.fs]]$marker_cols)) {
        cx <- tr$fn(reads[[.fs]]$exprs[gates[[.fs]]$masks$single_cells,
                                       reads[[.fs]]$marker_cols[[m]]], m)
        .kind <- "fmo_q995"
        .fmo_used[[length(.fmo_used) + 1L]] <- data.frame(
          sample_id = s, marker = m,
          fmo_threshold = as.numeric(stats::quantile(cx, 0.995, na.rm = TRUE)),
          fmo_sample = .fs, stringsAsFactors = FALSE)
      }
      rr <- resolve_threshold(m, tmat[g$masks$cd45_pos, m], cfg_thr[[m]]$threshold,
                              cx, override = .ov, control_kind = .kind)
      thr[[m]] <- rr$threshold; tdet[[m]] <- rr
      qdens[[m]] <- tmat[g$masks$cd45_pos, m]
      thr_rows[[length(thr_rows) + 1L]] <- data.frame(
        sample_id = s, panel = pn, marker = m, threshold = rr$threshold,
        source = rr$source, needs_review = rr$needs_review, cofactor = cf,
        # Carried on every row, but dropped again before the table is written
        # unless the run actually declares an override. Two all-NA columns are
        # still a change to a published file, and a run with no overrides must
        # produce the table it always produced.
        override_reason = rr$override_reason %||% NA_character_,
        override_by = rr$override_by %||% NA_character_,
        stringsAsFactors = FALSE)
    }
    hi <- derive_intermediate_bounds(tmat, thr, g$masks$cd45_pos, spec)
    sp <- score_populations(tmat, thr, g$masks$cd45_pos, spec, hi_thr = hi)
    pops[[s]] <- list(scored = sp, thresholds = thr, hi_thresholds = hi,
                      details = tdet, tmat = tmat)
    qc_recon[[s]] <- list(sample_id = s, thresholds = tdet, marker_densities = qdens)

    # ---- release the raw expression matrix -----------------------------------
    #
    # THE PEAK THIS FLATTENS. Until here, every sample held TWO full matrices at
    # once: the raw `exprs` read from the file, and the transformed `tmat` just
    # derived from it. Both are events x channels doubles, so at 300,000 events
    # the pair costs roughly 60-70 MB per sample, and the loop accumulates them
    # for all 25 -- which is precisely where an uncapped run gets OOM-killed
    # (exit 137 at STEP 4, no message, because the kernel does not send one).
    #
    # `tmat` is genuinely needed downstream: it is what the embedding, the MFI
    # tables and the marker-state tests read. `exprs` is not -- after this line
    # nothing reads it, which save_session()'s existing keep_exprs = FALSE
    # default already asserted. Dropping it here rather than at the end of the
    # loop is the whole point: releasing it at the end would be after the peak,
    # not before it.
    #
    # TWO EXCEPTIONS, both required for correctness rather than tidiness:
    #   - an unstained control is still read by LATER samples in this panel, as
    #     the reference for unimodal markers (see resolve_threshold above), so
    #     its matrix must outlive its own iteration;
    #   - --keep-exprs asks for these matrices in session_state.RData.
    if (!isTRUE(opt$keep_exprs) && !identical(s, ctrl_s)) {
      reads[[s]]$exprs <- NULL
      rd <- NULL
    }
  }
  # The controls held back above are released now that every sample that could
  # reference them has been scored.
  if (!isTRUE(opt$keep_exprs))
    for (.cs in unlist(ctrl_by_panel)) if (!is.null(.cs)) reads[[.cs]]$exprs <- NULL
  fig_gating_qc(qc_recon, file.path(opt$outdir, "gating_qc.png"))

  if (!is.null(opt$write_config)) {
    derived <- list()
    for (s in names(pops)) for (m in names(pops[[s]]$details)) {
      d <- pops[[s]]$details[[m]]
      if (is.null(derived[[m]]))
        derived[[m]] <- list(threshold = round(d$threshold, 4), source = d$source,
                             needs_review = d$needs_review)
    }
    write_config(opt$write_config, derived, cofactors, spec, blocks, colors = fcs_colors())
    return(invisible(NULL))
  }

  # ---- gate placement uncertainty -------------------------------------------
  # Placed here because it needs the transformed matrices and the parent masks,
  # which are at their most available immediately after scoring.
  #
  # IT READS AND WRITES NOTHING BACK. `pops`, `thr_rows` and every threshold
  # already derived are untouched, so a run with this enabled and a run without
  # produce identical frequencies, MFIs and p-values. What it adds is a second
  # number beside each of them.
  #
  # ITS RNG IS BORROWED, NOT SPENT. run_gate_uncertainty() saves .Random.seed and
  # puts it back, for the reason spelled out at run_unsupervised_clusters(): this
  # function seeds once at the top and STEP 6's cell selection draws from that one
  # stream, so a step that consumed draws here would change which cells are
  # embedded and quietly redraw every UMAP in the run.
  # ---- what excluding the flagged intervals would cost ----------------------
  # The quantity that turns a QC flag into a decision. Computed only when the
  # events are still present: with --drop-unstable-events they are already gone
  # and the difference is zero by construction.
  if (!is.null(.aqc) && !is.null(.aqc$flagged) &&
      !isTRUE(opt$drop_unstable_events)) {
    .imp <- list()
    for (s in names(.aqc$flagged)) {
      fl <- .aqc$flagged[[s]]
      P <- pops[[s]]; g <- gates[[s]]
      if (is.null(P$tmat) || is.null(g) || is.null(fl)) next
      if (length(fl) != nrow(P$tmat) || !any(fl)) next
      d <- tryCatch(frequency_delta_if_cleaned(P$tmat, P$thresholds,
                                               g$masks$cd45_pos, spec, !fl),
                    error = function(e) NULL)
      if (is.null(d)) next
      d$sample_id <- s
      .imp[[length(.imp) + 1L]] <- d
    }
    if (length(.imp)) {
      .imp <- do.call(rbind, .imp)
      write.csv(.imp, file.path(opt$outdir, "acquisition_qc_impact.csv"),
                row.names = FALSE)
      .worst <- max(abs(.imp$pct_delta_if_cleaned), na.rm = TRUE)
      log_msg("wrote acquisition_qc_impact.csv (largest movement if the flagged ",
              "intervals were excluded: ", round(.worst, 3), " percentage ",
              "points). Compare against that population's u_pct_points before ",
              "deciding the file needs re-acquiring")
    }
  }

  unc <- NULL
  if (!isTRUE(opt$no_uncertainty)) {
    log_step("STEP 4b - gate placement uncertainty")
    .u_t0 <- Sys.time()
    unc <- tryCatch(
      run_gate_uncertainty(pops, gates, verdicts, fpr$assignment, spec,
                           B = opt$uncertainty_boot %||% 100L,
                           seed = opt$seed %||% 42L,
                           max_events = opt$uncertainty_max_events %||% 20000L,
                           lod_events = opt$lod_events %||% 20L,
                           loq_events = opt$loq_events %||% 50L),
      error = function(e) {
        log_msg("  WARNING uncertainty analysis failed: ", conditionMessage(e),
                ", every other output is unaffected")
        NULL
      })
    if (!is.null(unc) && !is.null(unc$thresholds)) {
      write.csv(unc$thresholds,
                file.path(opt$outdir, "threshold_uncertainty.csv"),
                row.names = FALSE)
      ut <- unc$thresholds
      log_msg("wrote threshold_uncertainty.csv (", nrow(ut), " threshold(s), ",
              round(as.numeric(difftime(Sys.time(), .u_t0, units = "secs"))),
              "s)")
      # A threshold whose bootstrap rarely finds the valley at all is not
      # imprecise, it is unresolved, and that is worth a line in the log rather
      # than a wide interval nobody reads.
      shaky <- ut[is.finite(ut$bootstrap_valley_rate) &
                    ut$bootstrap_valley_rate < 0.8, , drop = FALSE]
      if (nrow(shaky))
        log_msg("  NOTE ", nrow(shaky), " threshold(s) were found in fewer than ",
                "80% of resamples, so the cut is barely determined by the data. ",
                "Worst: ",
                paste(utils::head(paste0(shaky$sample_id[order(shaky$bootstrap_valley_rate)],
                                         "/",
                                         shaky$marker[order(shaky$bootstrap_valley_rate)]),
                                  4L), collapse = ", "))
    }
    if (!is.null(unc) && !is.null(unc$budget)) {
      write.csv(unc$budget, file.path(opt$outdir, "uncertainty_budget.csv"),
                row.names = FALSE)
      bb <- stats::aggregate(u_pct_points ~ term, unc$budget, median)
      bb <- bb[order(-bb$u_pct_points), , drop = FALSE]
      log_msg("wrote uncertainty_budget.csv (largest median contribution: ",
              bb$term[1], ", ", round(bb$u_pct_points[1], 3),
              " percentage points)")
    }
    # Counting is not a gate term and is deliberately absent from the budget
    # above, which answers "which threshold". It is reported here and on the
    # frequency table instead, because a population can sit behind a perfectly
    # placed cut and still rest on too few events to mean anything.
    if (!is.null(unc) && !is.null(unc$frequencies) &&
        "detection" %in% names(unc$frequencies)) {
      .df <- qc_pass_rows(unc$frequencies)
      if (!is.null(.df) && nrow(.df)) {
        .nq <- sum(.df$detection %in% c("below LOD", "detected, below LOQ"))
        if (.nq)
          log_msg("NOTE ", .nq, " of ", nrow(.df), " population-sample values ",
                  "rest on fewer than ", opt$loq_events %||% 50L, " events and ",
                  "are not quantifiable at this acquisition depth (",
                  sum(.df$detection == "below LOD"), " below the limit of ",
                  "detection). See the detection column of ",
                  "population_frequencies.csv; raising --max-events-per-file ",
                  "lowers the limit, changing the gating does not")
      }
    }
  }

  # ---- metadata -------------------------------------------------------------
  patients <- NULL
  if (!is.null(.sheet)) {
    patients <- .sheet$patients
    if (!is.null(patients)) {
      log_step("STEP 5 - patient metadata")
      write.csv(patients, file.path(opt$outdir, "patient_metadata_english.csv"),
                row.names = FALSE, na = "")
      log_msg("wrote patient_metadata_english.csv (from the sample sheet)")
    }
  } else if (!is.null(opt$patient_table)) {
    log_step("STEP 5 - patient metadata")
    cm <- cfg$metadata$column_map %||% default_column_map()
    vm <- cfg$metadata$value_translations %||% default_value_map()
    refd <- if (!is.null(opt$reference_date)) {
      rr <- as.Date(opt$reference_date, format = "%Y-%m-%d")
      if (is.na(rr)) stop("--reference-date must be YYYY-MM-DD, got: ",
                          opt$reference_date)
      rr
    } else Sys.Date()
    log_msg("  reference date for age derivation: ", format(refd),
            if (is.null(opt$reference_date))
              " (today - pass --reference-date for reproducible ages)" else "")
    patients <- load_patient_table(opt$patient_table, cm, vm,
                                   reference_date = refd)
    write.csv(patients, file.path(opt$outdir, "patient_metadata_english.csv"),
              row.names = FALSE, na = "")
    log_msg("wrote patient_metadata_english.csv")
  }

  # ---- embed per panel ------------------------------------------------------
  log_step("STEP 6 - embedding")
  embeddings <- list(); all_cells <- list()
  for (p in fpr$panels) {
    inc <- p$samples[vapply(p$samples, function(s) isTRUE(verdicts[[s]]$include), TRUE)]
    if (!length(inc)) {
      log_msg(p$name, ": no sample passed staining QC, embedding skipped")
      next
    }
    log_msg(p$name, ": embedding ", length(inc), " sample(s)")
    # Built-in default is a panel-agnostic lineage-marker list, not this
    # study's panel specifically -- select_umap_features() already falls back
    # to every eligible marker when fewer than 2 of these are present. On a
    # foreign panel that happens to share 2+ names with this list (common
    # ones like CD3/CD4) but is otherwise unrelated, that fallback never
    # fires and the embedding would silently run on just those markers.
    # --umap-markers overrides the list outright; --umap-markers-all skips
    # lineage preference entirely and uses every eligible marker.
    # opt[["umap_markers", exact = TRUE]], not opt$umap_markers: "umap_markers"
    # is an exact prefix of the sibling flag's dest "umap_markers_all", and $
    # partial-matches onto it (returning that flag's logical value) whenever
    # --umap-markers is absent, since optparse drops NULL-default options
    # from the list instead of keying them to NULL.
    um_opt <- opt[["umap_markers", exact = TRUE]]
    umap_prefer <- if (!is.null(um_opt)) trimws(strsplit(um_opt, ",")[[1]])
    else c("CD3","CD4","CD8","CD14","CD19","CD56","CD16","HLA-DR",
          "TCR-Vd1","TCR-Vd2","CD25","CD127","CD57")
    feats <- select_umap_features(
      setNames(seq_along(p$markers), p$markers),
      prefer = intersect(umap_prefer, p$markers),
      exclude = na.omit(c(gates[[inc[1]]]$viability_marker)),
      lineage_only = !isTRUE(opt$umap_markers_all))
    navail <- vapply(inc, function(s) sum(gates[[s]]$masks$cd45_pos), integer(1))
    planned <- plan_subsample(setNames(navail, inc), cap = opt$cells_per_sample,
                              total_cap = opt$max_cells)
    .rare <- identical(opt$subsample %||% "uniform", "rare")
    rows <- list()
    for (s in inc) {
      idx <- which(gates[[s]]$masks$cd45_pos)
      tm <- pops[[s]]$tmat
      keep <- intersect(feats, colnames(tm))
      if (.rare) {
        # Inverse-density draw. Changes WHICH cells are embedded and therefore
        # every embedding figure, which is why it is opt-in. sampling_weight
        # carries the reciprocal inclusion probability so anything computed from
        # the embedded cells can be weighted back to the true composition.
        dr <- draw_subsample_rare(idx, tm[idx, keep, drop = FALSE],
                                  min(length(idx), planned[[s]]),
                                  seed = opt$seed %||% 42L)
        take <- dr$idx; wts <- dr$weight
      } else {
        take <- sort(sample(idx, min(length(idx), planned[[s]])))
        wts <- NULL
      }
      d <- as.data.frame(tm[take, keep, drop = FALSE])
      d$sample_id <- s
      d$population_label <- pops[[s]]$scored$labels[take]
      d$panel <- p$name
      d$event_index <- take
      if (!is.null(wts)) d$sampling_weight <- round(wts, 5)
      rows[[s]] <- d
    }
    cells <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
    fcols <- intersect(feats, colnames(cells))
    .fmat <- as.matrix(cells[, fcols, with = FALSE])

    # ---- optional batch correction ------------------------------------------
    # Applied to the assembled marker matrix, after gating and scoring and
    # before the embedding. That position is deliberate: the gate hierarchy and
    # every threshold are derived PER SAMPLE and are already batch-local, so
    # correcting before them would only fight the per-sample derivation, while
    # correcting after the embedding would fix the picture and leave the
    # statistics untouched.
    #
    # correct_batch() refuses when batch and group are confounded. That refusal
    # is the feature, not an obstacle to it -- see batch-correct.R.
    batch_info <- NULL
    if (isTRUE(opt$correct_batch)) {
      bcol <- opt$batch_column
      if (is.null(bcol)) {
        log_msg("  --correct-batch needs --batch-column; skipping correction")
      } else {
        # Both vectors come from the sample map rather than from group_of(),
        # which is not resolved until the statistics stage several hundred lines
        # below. The confounding check only needs the labels, and taking them
        # from the map keeps this step independent of how groups are later
        # merged with the patient table.
        .col <- function(nm) if (!is.null(smap) && !is.null(nm) && nm %in% names(smap))
          smap[[nm]][match(cells$sample_id, smap$sample_id)] else NULL
        bvec <- .col(bcol)
        gvec <- .col(opt$group_column)
        if (is.null(bvec) || all(is.na(bvec))) {
          log_msg("  batch column '", bcol, "' not found in the sample map; ",
                  "skipping correction")
        } else {
          .bm <- tolower(trimws(opt$batch_method %||% "quantile"))
          if (!.bm %in% c("quantile", "cluster", "cytonorm")) {
            log_msg("  --batch-method '", .bm, "' not recognised; using ",
                    "'quantile'. Valid: quantile, cluster, cytonorm")
            .bm <- "quantile"
          }
          bres <- correct_batch(.fmat, bvec, gvec,
                                max_cramers_v = opt$batch_max_cramers_v %||% 0.6,
                                force = isTRUE(opt$force_batch_correction),
                                method = .bm,
                                cluster_k = opt$batch_cluster_k %||% 10L,
                                seed = opt$seed %||% 42L)
          .fmat <- bres$tmat
          for (cc in colnames(.fmat)) data.table::set(cells, j = cc, value = .fmat[, cc])
          batch_info <- bres
          # The decision was previously held in memory and never written, so a
          # reader could not tell from the outputs whether a correction had run,
          # which method fitted it, or why it was refused. One row, always
          # written once --correct-batch is set, including on a refusal.
          write.csv(data.frame(
            method        = bres$method,
            fitted_as     = bres$how %||% NA_character_,
            corrected     = bres$corrected,
            cramers_v     = bres$cramers_v,
            max_cramers_v = opt$batch_max_cramers_v %||% 0.6,
            forced        = isTRUE(opt$force_batch_correction),
            n_markers     = length(bres$markers),
            markers       = paste(bres$markers, collapse = ";"),
            reason        = bres$reason,
            stringsAsFactors = FALSE),
            file.path(opt$outdir, "batch_correction.csv"), row.names = FALSE)
          log_msg("wrote batch_correction.csv")
        }
      }
    }

    # ---- project into a saved embedding, or train a new one -----------------
    # Projection keeps coordinates comparable between runs, which is the whole
    # point of --umap-model: cluster 4 stays where cluster 4 was. It falls back
    # to training rather than failing, because every reason a projection is
    # refused (missing model, missing marker, incompatible uwot) is a reason to
    # produce a correct fresh embedding, not no embedding.
    emb <- NULL
    if (!is.null(opt$umap_model) && TRUE) {
      .saved <- load_umap_model(paste0(opt$umap_model,
                                       if (length(fpr$panels) > 1L)
                                         paste0(".", p$name) else ""))
      .pc <- if (!is.null(.saved))
        project_umap(.saved, .fmat,
                     n_threads = if (opt$threads > 0) opt$threads else 1L) else NULL
      if (!is.null(.pc)) {
        log_msg("  projected ", nrow(.pc), " cells into the saved embedding ",
                "(coordinates are comparable with the run that trained it)")
        emb <- list(coords = .pc, model = NULL, scale_params = .saved$scale_params,
                    params = c(.saved$meta, list(projected = TRUE,
                                                 n_cells = nrow(.pc))))
      }
    }
    if (is.null(emb))
      emb <- run_umap(.fmat,
                      n_neighbors = opt$n_neighbors, min_dist = opt$min_dist,
                      n_threads = if (opt$threads > 0) opt$threads else NULL,
                      seed = opt$seed,
                      ret_model = !is.null(opt$save_umap_model))
    if (!is.null(opt$save_umap_model) && !is.null(emb$model) &&
        TRUE) {
      # One model per panel when there is more than one: populations scored under
      # different marker sets do not share an embedding, so a single file would
      # be silently wrong for all but the first.
      .mp <- paste0(opt$save_umap_model,
                    if (length(fpr$panels) > 1L) paste0(".", p$name) else "")
      try(save_umap_model(emb$model, .mp, emb$scale_params, fcols,
                          meta = emb$params), silent = TRUE)
    }
    # names must match the figure layer's aesthetics and the documented
    # cells_umap.csv schema; run_umap() returns them as umap_1/umap_2
    cells$umap_1 <- emb$coords[, 1]; cells$umap_2 <- emb$coords[, 2]
    # attach patient covariates
    # The join needs patient_id on BOTH sides. Check the patient table too, not
    # just the sample map: if the table's identifier header matched no entry in the
    # column map, `patients` arrives with every other column intact and no
    # patient_id, and the join below then fails deep inside setNames() with
    # "attempt to set an attribute on NULL" -- an error that says nothing about the
    # actual cause. Name the cause and list what the table does carry.
    if (!is.null(patients) && !"patient_id" %in% names(patients)) {
      stop("the patient table has no usable patient identifier column.\n",
           "  Columns recognised: ",
           paste(names(patients), collapse = ", "), "\n",
           "  Add the header spelling to column_map$patient_id in your config, ",
           "or rename the column to 'patient_id'.")
    }
    if (!is.null(patients) && !is.null(smap) && "patient_id" %in% names(smap)) {
      key <- smap[, intersect(c("sample_id","well","patient_id","timepoint"), names(smap)),
                  drop = FALSE]
      key$sample_id <- key$sample_id %||% key$well
      cells <- merge(cells, key, by = "sample_id", all.x = TRUE, sort = FALSE)

      # Match patient_id case- and whitespace-insensitively, then adopt the
      # PATIENT TABLE's spelling before merging.
      #
      # WHY: the sample map is typed by hand at the bench, the patient table is
      # exported from a clinical system, and "HS5_004" vs "hs5_004" is the normal
      # result. A plain merge() on the raw strings is case-sensitive, so such a
      # pair silently produces all-NA covariates -- every patient panel in the
      # overview figure comes out blank with no error anywhere. Normalising the
      # JOIN KEY while keeping the clinical table's canonical spelling in the
      # output means the figures label patients the way the clinical record does.
      pmap <- setNames(patients$patient_id, norm_id(patients$patient_id))
      canon <- unname(pmap[norm_id(cells$patient_id)])
      # Report before overwriting: a mismatch the user should know about.
      raw_ok  <- !is.na(cells$patient_id) & cells$patient_id %in% patients$patient_id
      fuzzy   <- !is.na(canon) & !raw_ok
      unmatched <- unique(cells$patient_id[!is.na(cells$patient_id) & is.na(canon)])
      if (any(fuzzy))
        log_msg("  patient_id matched ignoring case/whitespace for ",
                length(unique(cells$patient_id[fuzzy])),
                " sample-map value(s); using the patient table's spelling")
      if (length(unmatched))
        log_msg("  WARNING: no patient-table row for patient_id: ",
                paste(utils::head(unmatched, 5), collapse = ", "),
                if (length(unmatched) > 5) sprintf(" (+%d more)", length(unmatched) - 5) else "",
                ", covariates will be NA for these samples")
      cells$patient_id[!is.na(canon)] <- canon[!is.na(canon)]
      cells <- merge(cells, patients, by = "patient_id", all.x = TRUE, sort = FALSE)
    }

    # The group column, from the sample map, when the patient table did not
    # supply it. Everything downstream that splits cells by study group reads it
    # off THIS table: umap_density_by_group.png, umap_overview_by_group.png and
    # the per-group FlowJo files. Without this a run whose grouping lives in the
    # sample map -- which is the whole design when there is no clinical metadata
    # -- produced neither figure and no per-group export, and said only that a
    # 'cohort' column was missing, naming a table the user may never have had.
    .gcol <- opt$group_column %||% "cohort"
    if (!.gcol %in% names(cells) && !is.null(smap) &&
        .gcol %in% names(smap) && "sample_id" %in% names(smap)) {
      cells[[.gcol]] <- as.character(smap[[.gcol]])[
        match(cells$sample_id, smap$sample_id)]
    }

    embeddings[[p$name]] <- list(cells = as.data.frame(cells), features = fcols,
                                 params = emb$params, samples = inc)
    all_cells[[p$name]] <- as.data.frame(cells)

    sfx <- if (length(fpr$panels) > 1L) paste0("_", p$name) else ""
    # Drop the "Other CD45+" leftover from the FIGURE data only (see --other).
    # embeddings/all_cells keep every cell, so cells_umap.csv, the FlowJo export
    # and every table below are untouched by this.
    if (!isTRUE(opt$include_other)) {
      .ec <- embeddings[[p$name]]$cells
      .keep <- !grepl("^Other|unclassified", .ec$population_label, ignore.case = TRUE)
      if (any(.keep) && !all(.keep)) {
        log_msg("  hiding '", paste(unique(.ec$population_label[!.keep]), collapse = "', '"),
                "' from the UMAP figures (", sum(!.keep), " of ", nrow(.ec),
                " cells), pass --other to draw it; tables are unaffected")
        embeddings[[p$name]]$cells <- .ec[.keep, , drop = FALSE]
      }
    }
    fig_umap_overview(embeddings[[p$name]]$cells,
                      file.path(opt$outdir, paste0("umap_overview", sfx, ".png")),
                      panel_label = p$name, feature_cols = fcols)
    fig_marker_grid(embeddings[[p$name]]$cells, fcols,
                    file.path(opt$outdir, paste0("umap_markers", sfx, ".png")),
                    panel_label = p$name)
    fig_density_by_sample(embeddings[[p$name]]$cells,
                          file.path(opt$outdir, paste0("umap_density", sfx, ".png")),
                          panel_label = p$name)

    gcol <- opt$group_column %||% "cohort"
    .has_group <- gcol %in% names(embeddings[[p$name]]$cells) &&
      length(unique(stats::na.omit(embeddings[[p$name]]$cells[[gcol]]))) > 1L

    # One UMAP per marker, at full size, in its own folder. Written on EVERY
    # run, not only grouped ones. Where a group column resolves to 2+ groups
    # each marker is faceted by it; where it does not, the same folder holds one
    # unfaceted panel per marker.
    #
    # WHY IT IS NOT COVERED BY THE FIGURES ABOVE. umap_markers.png colours by
    # intensity but pools the groups and shrinks every marker into a grid cell.
    # umap_overview_by_group.png splits by group but colours by population or
    # covariate, never by intensity. The question "is this marker brighter, or
    # somewhere else, in one group" needs both at once, and on an ungrouped run
    # the same file is still the only full-size view of one marker.
    #
    # On by default, following the rule for anything that only ADDS an output:
    # an analysis you have to know to ask for is one nobody runs.
    # --no-marker-group-umaps turns it off, for a panel wide enough that one PNG
    # per marker is a folder you do not want.
    #
    # try(): this runs before .ext_ok() is defined, and a figure failure here
    # would otherwise take down a run whose tables are already written.
    if (!isTRUE(opt$no_marker_group_umaps)) {
      .mg <- try(fig_marker_umaps_by_group(
        embeddings[[p$name]]$cells, fcols,
        file.path(opt$outdir, paste0("marker_umaps_by_group", sfx)),
        group_col = if (.has_group) gcol else NULL,
        panel_label = p$name), silent = TRUE)
      if (inherits(.mg, "try-error"))
        log_msg("  NOTE marker UMAPs failed, continuing: ",
                conditionMessage(attr(.mg, "condition")))
    }

    # Same comparison, one panel per GROUP instead of per sample -- the combined
    # embedding (umap_overview.png) already shows every group overlaid; this is
    # the side-by-side view of the same shared axes, so a group-specific shift
    # or dropout is visible without hunting through overlaid colours.
    if (gcol %in% names(embeddings[[p$name]]$cells)) {
      n_grp <- length(unique(stats::na.omit(embeddings[[p$name]]$cells[[gcol]])))
      if (n_grp > 1L) {
        fig_density_by_sample(
          embeddings[[p$name]]$cells,
          file.path(opt$outdir, paste0("umap_density_by_group", sfx, ".png")),
          panel_label = p$name, facet_by = gcol)
        # Combined + one column per group, for every colouring umap_overview.png
        # draws (population, sample, covariates) -- not a replacement for either
        # figure above, an additional side-by-side view across them.
        fig_umap_overview_by_group(
          embeddings[[p$name]]$cells,
          file.path(opt$outdir, paste0("umap_overview_by_group", sfx, ".png")),
          group_col = gcol, panel_label = p$name, feature_cols = fcols)
        # FlowJo-style multigraph overlay. The two figures above answer "where is
        # this population and how big is it"; this one answers "what does it
        # express", per cluster per subcluster per marker, one curve per group.
        #
        # Subclusters are fitted on the REFERENCE group only (see
        # subcluster_by_reference) so each one is a compartment of normal
        # biology that the patient studies are then projected onto -- the
        # comparison the study exists to make. With no --reference-group there
        # is nothing to be normal relative to, so the clusters stay whole and
        # the figure is exactly what it was before.
        ec  <- embeddings[[p$name]]$cells
        # --auto-subcluster-k replaces the fixed k with a per-population choice
        # by mean silhouette on the reference cells. Opt-in, because it CHANGES
        # the subcluster lettering -- "4b" under auto-k is not "4b" under fixed k
        # -- and every published overlay and shift table is indexed by that
        # lettering. The full score curve is written so the choice is auditable
        # rather than implicit.
        .k_sub <- 3L
        if (isTRUE(opt$auto_subcluster_k) && !is.null(opt$reference_group) &&
            TRUE) {
          .ck <- choose_subcluster_k(ec, markers = fcols, group_col = gcol,
                                     reference = opt$reference_group,
                                     seed = opt$seed)
          if (!is.null(.ck)) {
            .k_sub <- .ck$k
            write.csv(.ck$curve, file.path(opt$outdir,
                      paste0("subcluster_k_selection", sfx, ".csv")), row.names = FALSE)
            log_msg("  auto-k by silhouette: ",
                    paste(sprintf("%s=%d", names(.ck$k), .ck$k), collapse = ", "))
          } else {
            log_msg("  NOTE --auto-subcluster-k: no population had enough reference ",
                    "cells to score; falling back to k = 3")
          }
        }
        sub <- if (!is.null(opt$reference_group))
          subcluster_by_reference(ec, markers = fcols, group_col = gcol,
                                  reference = opt$reference_group, k = .k_sub,
                                  seed = opt$seed)
        else NULL
        if (!is.null(sub))
          log_msg("  subclustered on reference '", opt$reference_group, "': ",
                  length(unique(paste(ec$population_label, sub))),
                  " cluster-subcluster compartment(s)")
        fig_multigraph_overlay(
          ec, file.path(opt$outdir, paste0("umap_multigraph_overlay", sfx, ".png")),
          markers = fcols, group_col = gcol, panel_label = p$name, subcluster = sub,
          reference = opt$reference_group)
        if (!is.null(sub)) {
          sh <- stats_subcluster_shifts(ec, markers = fcols, subcluster = sub,
                                        group_col = gcol, reference = opt$reference_group)
          if (!is.null(sh)) {
            write.csv(sh, file.path(opt$outdir,
                      paste0("subcluster_marker_shifts", sfx, ".csv")), row.names = FALSE)
            log_msg("wrote subcluster_marker_shifts", sfx, ".csv (", nrow(sh),
                    " cluster x subcluster x marker x study comparison(s), ",
                    "sorted by |Cliff's delta|)")
          }
        }
      }
    }
  }

  # ---- tables ---------------------------------------------------------------
  log_step("STEP 7 - writing tables")
  gate_counts <- do.call(rbind, lapply(gates, `[[`, "counts"))
  write.csv(gate_counts, file.path(opt$outdir, "gate_counts.csv"), row.names = FALSE)

  qc <- do.call(rbind, lapply(names(verdicts), function(s) data.frame(
    sample_id = s, panel = fpr$assignment[[s]], file = reads[[s]]$file,
    n_events = reads[[s]]$n_events,
    pct_cd45_of_live = round(verdicts[[s]]$pct_cd45_of_live, 3),
    verdict = verdicts[[s]]$verdict, qc_status = verdicts[[s]]$qc_status %||% NA,
    included = verdicts[[s]]$include,
    stringsAsFactors = FALSE)))
  write.csv(qc, file.path(opt$outdir, "staining_qc.csv"), row.names = FALSE)

  # ---- cross-sample threshold consistency -----------------------------------
  # WHY THIS CHECK EXISTS: needs_review is a per-sample flag about how ONE
  # threshold was derived (it is TRUE only for the quantile fallback). It cannot
  # see the failure mode that actually corrupts a cohort comparison: a threshold
  # that was found by a legitimate valley, looks entirely plausible on its own
  # sample, and yet sits nowhere near the same marker's threshold in every other
  # sample of the same panel.
  #
  # That situation is not a subtle bias. When a marker's threshold collapses
  # toward (or below) zero in one tube while its peers sit in the normal range,
  # essentially every event in that tube is called positive, so a population is
  # reported near 0% in some samples and high in others -- an all-or-nothing
  # pattern across the cohort that reads as a dramatic biological effect and is
  # purely an artefact of gate placement. Per-sample gating is correct and
  # necessary (gain and staining drift between acquisitions), but the resulting
  # thresholds are only comparable if they are actually comparable, and nothing
  # was checking that.
  #
  # Compare each threshold against the MEDIAN OF THE OTHER SAMPLES for the same
  # marker within the same panel (leave-one-out, so one bad tube cannot drag the
  # reference it is being judged against), scaled by the MAD of those others so
  # the tolerance adapts to how variable that marker legitimately is.
  thr_all <- do.call(rbind, thr_rows)
  # The override columns exist on every row but only earn their place in the
  # written table when the run declares at least one. A run with none must write
  # the file it has always written.
  if (!is.null(thr_all) && !length(cfg_ovr))
    thr_all <- thr_all[, setdiff(names(thr_all), c("override_reason", "override_by")),
                       drop = FALSE]
  if (length(cfg_ovr)) {
    report_unused_overrides(cfg_ovr, .ovr_applied)
    if (length(.ovr_applied))
      log_msg("applied ", length(.ovr_applied), " manual threshold override(s); ",
              "they are recorded with their reason in thresholds_used.csv and ",
              "listed in the run manifest")
  }
  thr_all$scale_outlier <- FALSE
  thr_all$peer_median   <- NA_real_
  thr_all$robust_z      <- NA_real_
  if (!is.null(thr_all) && nrow(thr_all)) {
    if (!"panel" %in% names(thr_all))
      thr_all$panel <- unlist(fpr$assignment[thr_all$sample_id], use.names = FALSE)
    for (key in unique(paste(thr_all$panel, thr_all$marker, sep = "\r"))) {
      idx <- which(paste(thr_all$panel, thr_all$marker, sep = "\r") == key)
      v <- thr_all$threshold[idx]
      # Need enough peers for "the others" to mean anything.
      if (sum(is.finite(v)) < 4L) next
      for (j in seq_along(idx)) {
        others <- v[-j]; others <- others[is.finite(others)]
        if (length(others) < 3L || !is.finite(v[j])) next
        med <- median(others)
        # MAD floored so a marker with near-identical peer thresholds does not
        # divide by ~0 and flag every trivial difference.
        s <- max(mad(others), 0.25)
        z <- abs(v[j] - med) / s
        thr_all$peer_median[idx[j]] <- round(med, 4)
        thr_all$robust_z[idx[j]]    <- round(z, 2)
        thr_all$scale_outlier[idx[j]] <- z > 3.5
      }
    }
    n_out <- sum(thr_all$scale_outlier, na.rm = TRUE)
    if (n_out) {
      bad <- thr_all[which(thr_all$scale_outlier), ]
      log_msg("WARNING ", n_out, " threshold(s) inconsistent with the same marker ",
              "in other samples of the same panel, per-sample frequencies for the ",
              "affected populations are NOT comparable. See scale_outlier in ",
              "thresholds_used.csv:")
      for (i in seq_len(min(nrow(bad), 12L)))
        log_msg("    ", bad$sample_id[i], " / ", bad$marker[i], ": ",
                round(bad$threshold[i], 2), " vs peer median ",
                bad$peer_median[i], " (robust z ", bad$robust_z[i], ")")
      if (nrow(bad) > 12L) log_msg("    ... and ", nrow(bad) - 12L, " more")
      # Per-marker summary: which markers are unstable, and over what range.
      sm <- do.call(rbind, lapply(split(thr_all, paste(thr_all$panel, thr_all$marker)),
        function(d) data.frame(panel = d$panel[1], marker = d$marker[1],
          n = nrow(d), median = round(median(d$threshold, na.rm = TRUE), 3),
          min = round(min(d$threshold, na.rm = TRUE), 3),
          max = round(max(d$threshold, na.rm = TRUE), 3),
          n_outlier = sum(d$scale_outlier, na.rm = TRUE),
          n_fallback = sum(d$source == "quantile_fallback"),
          stringsAsFactors = FALSE)))
      write.csv(sm[order(-sm$n_outlier), ],
                file.path(opt$outdir, "threshold_scale_qc.csv"), row.names = FALSE)
    }
  }
  write.csv(thr_all, file.path(opt$outdir, "thresholds_used.csv"),
            row.names = FALSE)

  # ---- spillover spreading ---------------------------------------------------
  # Diagnostic only, and the explanation for a class of quantile_fallback this
  # package could otherwise only report without a cause.
  if (!isTRUE(opt$no_spreading)) {
    .sp <- tryCatch(run_spreading_report(pops, gates, thr_all, fpr$assignment,
                                         max_samples = opt$spreading_max_samples %||% 8L),
                    error = function(e) NULL)
    if (!is.null(.sp)) {
      write.csv(.sp$pairs, file.path(opt$outdir, "spreading_pairs.csv"),
                row.names = FALSE)
      write.csv(.sp$receivers, file.path(opt$outdir, "spreading_receivers.csv"),
                row.names = FALSE)
      .nb <- sum(grepl("^unresolved and heavily spread", .sp$receivers$verdict))
      log_msg("wrote spreading_pairs.csv and spreading_receivers.csv (",
              sum(.sp$pairs$substantial, na.rm = TRUE), " of ", nrow(.sp$pairs),
              " channel pair(s) widen the receiver's negative population by 25% ",
              "or more)")
      if (.nb)
        log_msg("  NOTE ", .nb, " marker(s) both fail to resolve a density ",
                "minimum in most samples AND receive substantial spreading. ",
                "For those the cut is unresolved because of the panel, and no ",
                "gating strategy recovers it. See spreading_receivers.csv")
    }
  }

  # ---- derived cut against its FMO-anchored equivalent -----------------------
  # The diagnostic the FMO feature exists for. Scaled by the threshold's own
  # uncertainty where that is available, because a gap of 0.3 units means nothing
  # until you know whether the cut moves by 0.05 or by 0.5 under resampling. Runs
  # after thr_all is complete and after STEP 4b, so `unc` is populated.
  if (length(.fmo_used)) {
    .fa <- tryCatch(fmo_agreement(thr_all, do.call(rbind, .fmo_used),
                                  unc = unc$thresholds),
                    error = function(e) NULL)
    if (!is.null(.fa) && nrow(.fa)) {
      write.csv(.fa, file.path(opt$outdir, "fmo_agreement.csv"), row.names = FALSE)
      .nd <- sum(grepl("^derived cut well", .fa$verdict))
      log_msg("wrote fmo_agreement.csv (", nrow(.fa), " comparison(s), ",
              .nd, " where the derived cut and the FMO disagree by more than ",
              "the cut's own uncertainty can explain)")
      if (.nd)
        log_msg("  NOTE a derived cut far ABOVE its FMO is discarding real ",
                "signal; far BELOW, it is calling spillover positive. See the ",
                "verdict column")
    }
  }

  # Two DIFFERENT exclusion reasons are recorded separately, because conflating
  # them produces a false statement in every figure caption:
  #   is_control = TRUE  -- an unstained control tube. Its CD45+ parent is a
  #     quantile-fallback slice, so every percentage is a fraction of an arbitrary
  #     top-N% cut.
  #   qc_status = "failed" -- a BIOLOGICAL sample whose staining or CD45 gate did
  #     not work. Also not comparable, but it is a patient, not a control, and a
  #     caption that reports it as a control misstates the cohort.
  # Rows are retained either way (useful for spotting reagent carry-over and for
  # auditing exclusions) but neither may be compared against stained samples,
  # which is why the figures plot only `qc_status == "pass"`.
  freq <- do.call(rbind, lapply(names(pops), function(s) {
    m <- pops[[s]]$scored$masks; par <- sum(gates[[s]]$masks$cd45_pos)
    if (!length(m)) return(NULL)
    data.frame(sample_id = s, panel = fpr$assignment[[s]], population = names(m),
               count = vapply(m, sum, integer(1)),
               pct_of_cd45_pos = 100 * vapply(m, sum, integer(1)) / max(1L, par),
               is_control = isTRUE(verdicts[[s]]$is_control),
               qc_status = verdicts[[s]]$qc_status %||% "pass",
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  # ---- absolute concentrations, when an external blood count is supplied ------
  # `count` in this table is an EVENT count, not a cell number: it counts events
  # the instrument recorded, which depends on how long the tube was run and on
  # --max-events-per-file, so it is not comparable between samples and must never
  # be plotted as if it were a biological quantity.
  #
  # A frequency is comparable but compositional -- populations are constrained to
  # sum to 100%, so one lineage expanding forces every other to fall. Deciding
  # whether a cohort has FEWER NK cells or simply MORE of something else requires
  # an absolute scale, which a cytometer cannot supply because it does not measure
  # the volume it sampled. wbc_per_ul from a haemogram supplies it:
  #
  #     cells/uL(population) = pct_of_cd45_pos/100 * wbc_per_ul
  #
  # CD45+ is the denominator and CD45 marks leukocytes, so the total leukocyte
  # concentration is the correct scale factor for this hierarchy.
  # Built from `patients` + `smap`, NOT from the per-cell table: `cells` is scoped
  # to the embedding loop and holds only the LAST panel, so reading it here would
  # silently drop every sample in every other panel.
  freq$cells_per_ul <- NA_real_
  wbc_lookup <- NULL
  if (!is.null(patients) && "wbc_per_ul" %in% names(patients) &&
      !is.null(smap) && "patient_id" %in% names(smap)) {
    k <- smap[, intersect(c("sample_id", "well", "patient_id"), names(smap)),
              drop = FALSE]
    k$sample_id <- k$sample_id %||% k$well
    pw <- patients[, c("patient_id", "wbc_per_ul")]
    pw <- pw[!is.na(pw$wbc_per_ul), , drop = FALSE]
    if (nrow(pw)) {
      # Same case/whitespace-insensitive matching the covariate join uses, so a
      # bench-typed sample map and a clinically-exported table still line up.
      k$.j  <- norm_id(k$patient_id)
      pw$.j <- norm_id(pw$patient_id)
      m <- merge(k, pw[, c(".j", "wbc_per_ul")], by = ".j", all.x = TRUE)
      m <- m[!is.na(m$wbc_per_ul) & !is.na(m$sample_id), , drop = FALSE]
      if (nrow(m))
        wbc_lookup <- setNames(as.numeric(m$wbc_per_ul), m$sample_id)
    }
  }
  if (!is.null(wbc_lookup)) {
    hit <- freq$sample_id %in% names(wbc_lookup)
    freq$cells_per_ul[hit] <- round(
      freq$pct_of_cd45_pos[hit] / 100 * wbc_lookup[freq$sample_id[hit]], 2)
    n_s <- length(unique(freq$sample_id[hit]))
    log_msg("absolute concentrations derived for ", n_s, "/",
            length(unique(freq$sample_id)), " sample(s) from wbc_per_ul")
    if (n_s < length(unique(freq$sample_id)))
      log_msg("  NOTE samples without wbc_per_ul have cells_per_ul = NA; ",
              "absolute-count comparisons are restricted to the ", n_s,
              " sample(s) that have it")
  } else {
    # Stated once, plainly, because the absence changes which conclusions the
    # output can support -- it is not a missing nicety.
    log_msg("NOTE no wbc_per_ul column in the patient table, so no absolute cell ",
            "counts. `count` is an EVENT count (acquisition-dependent, capped by ",
            "--max-events-per-file) and is NOT a cell number: report ",
            "pct_of_cd45_pos, not count. Add wbc_per_ul (cells/uL, from a ",
            "haemogram) to the patient table to obtain cells/uL.")
  }
  # Uncertainty columns are APPENDED, never inserted, and only when the analysis
  # ran. Every column this table has always had keeps its name, its meaning and
  # its position, so a script reading it by name or by position is unaffected and
  # the values themselves are bit-identical to a run without --uncertainty.
  if (!is.null(freq) && !is.null(unc) && !is.null(unc$frequencies)) {
    .uf <- unc$frequencies
    .k <- match(paste(freq$sample_id, freq$population, sep = "\r"),
                paste(.uf$sample_id, .uf$population, sep = "\r"))
    freq$u_pct_points     <- .uf$u_pct_points[.k]
    freq$pct_lo           <- round(pmax(0, freq$pct_of_cd45_pos -
                                          .uf$u_pct_points[.k]), 4)
    freq$pct_hi           <- round(freq$pct_of_cd45_pos +
                                     .uf$u_pct_points[.k], 4)
    freq$u_n_terms        <- .uf$n_terms[.k]
    freq$u_n_terms_missing <- .uf$n_terms_missing[.k]
    # Counting uncertainty and the detection limits follow the same rule: append
    # only. u_pct_points above still means gate placement alone and still holds
    # the value it always did, and pct_lo/pct_hi are still built from it.
    if ("u_total_pct_points" %in% names(.uf)) {
      freq$n_parent_events        <- .uf$n_parent_events[.k]
      freq$u_counting_pct_points  <- .uf$u_counting_pct_points[.k]
      freq$u_total_pct_points     <- .uf$u_total_pct_points[.k]
      freq$pct_lo_total <- round(pmax(0, freq$pct_of_cd45_pos -
                                        .uf$u_total_pct_points[.k]), 4)
      freq$pct_hi_total <- round(freq$pct_of_cd45_pos +
                                   .uf$u_total_pct_points[.k], 4)
      freq$lod_pct    <- .uf$lod_pct[.k]
      freq$loq_pct    <- .uf$loq_pct[.k]
      freq$detection  <- .uf$detection[.k]
    }
  }
  write.csv(freq, file.path(opt$outdir, "population_frequencies.csv"), row.names = FALSE)

  unav <- do.call(rbind, lapply(names(pops), function(s) {
    u <- pops[[s]]$scored$unavailable
    if (!length(u)) return(NULL)
    data.frame(sample_id = s, population = names(u), reason = unlist(u),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  if (!is.null(unav))
    write.csv(unav, file.path(opt$outdir, "populations_unavailable.csv"), row.names = FALSE)

  mfi <- do.call(rbind, lapply(names(pops), function(s) {
    P <- pops[[s]]; tm <- P$tmat; m <- P$scored$masks
    if (!length(m)) return(NULL)
    do.call(rbind, lapply(names(m), function(nm) {
      sel <- m[[nm]]
      if (sum(sel) < 10L) return(NULL)
      # is_control / qc_status added here so qc_pass_rows() can do its job on this
      # table too. Without them a control tube's marker medians -- real numbers
      # describing an UNSTAINED sample -- enter the differential-state tests and
      # the phenotype heatmap silently, because qc_pass_rows() passes a table
      # through unchanged when the columns it filters on are absent. The
      # functional_markers table already carried both; this brings the MFI table
      # in line with it. Strictly additive: two new columns, no row or value
      # changes.
      data.frame(sample_id = s, panel = fpr$assignment[[s]], population = nm,
                 n_cells = sum(sel), marker = colnames(tm),
                 median_asinh = apply(tm[sel, , drop = FALSE], 2, median),
                 pct_positive = vapply(colnames(tm), function(mk)
                   100 * mean(tm[sel, mk] > P$thresholds[[mk]]), numeric(1)),
                 is_control = isTRUE(verdicts[[s]]$is_control),
                 qc_status = verdicts[[s]]$qc_status %||% "pass",
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }))
  write.csv(mfi, file.path(opt$outdir, "population_marker_mfi.csv"), row.names = FALSE)

  # ---- functional marker blocks, scoped as the gating strategy specifies -----
  # population_marker_mfi.csv reports EVERY marker in EVERY population, which is
  # useful but unscoped. This table applies the document's four scoping rules
  # (exhaustion outside monocytes, homing and activation in CD3+ subsets,
  # HLA-DR/BTN in the monocyte subsets) so each functional marker is reported
  # only where it is interpretable.
  fx <- do.call(rbind, lapply(names(pops), function(s) {
    P <- pops[[s]]; tm <- P$tmat; m <- P$scored$masks
    if (!length(m)) return(NULL)
    do.call(rbind, lapply(names(blocks), function(bn) {
      b <- blocks[[bn]]
      mk_use <- intersect(b$markers, colnames(tm))
      pop_use <- resolve_block_populations(b, spec, names(m))
      if (!length(mk_use) || !length(pop_use)) return(NULL)
      do.call(rbind, lapply(pop_use, function(nm) {
        sel <- m[[nm]]
        if (sum(sel) < 10L) return(NULL)
        data.frame(sample_id = s, panel = fpr$assignment[[s]], block = bn,
                   population = nm, n_cells = sum(sel), marker = mk_use,
                   median_asinh = vapply(mk_use, function(k)
                     median(tm[sel, k]), numeric(1)),
                   pct_positive = vapply(mk_use, function(k)
                     100 * mean(tm[sel, k] > P$thresholds[[k]]), numeric(1)),
                   is_control = isTRUE(verdicts[[s]]$is_control),
                   qc_status = verdicts[[s]]$qc_status %||% "pass",
                   row.names = NULL, stringsAsFactors = FALSE)
      }))
    }))
  }))
  if (!is.null(fx)) {
    write.csv(fx, file.path(opt$outdir, "functional_markers.csv"), row.names = FALSE)
    log_msg("wrote functional_markers.csv (", nrow(fx), " rows across ",
            length(unique(fx$block)), " blocks)")
    # Distinguish the two reasons a block marker has no rows. Reporting both as
    # "absent from all panels" blames the panel for a marker it contains, and
    # sends the reader to the acquisition template instead of to the gate spec.
    want_b   <- unique(unlist(lapply(blocks, `[[`, "markers")))
    in_panel <- unique(unlist(lapply(pops, function(P) colnames(P$tmat))))
    miss_b   <- setdiff(want_b, in_panel)
    empty_b  <- setdiff(intersect(want_b, in_panel), unique(fx$marker))
    if (length(miss_b))
      log_msg("NOTE functional markers absent from all panels: ",
              paste(sort(miss_b), collapse = ", "))
    if (length(empty_b))
      log_msg("NOTE functional markers present in the panel but with no reportable ",
              "population (fewer than 10 cells in every scoped population): ",
              paste(sort(empty_b), collapse = ", "))
  }

  # ---- derived abundance ratios (e.g. CD4:CD8), from the config's `ratios:` -
  rt <- compute_population_ratios(freq, ratios)
  if (!is.null(rt)) {
    write.csv(rt, file.path(opt$outdir, "population_ratios.csv"), row.names = FALSE)
    log_msg("wrote population_ratios.csv (", length(unique(rt$ratio)), " ratio(s), ",
            length(unique(rt$sample_id)), " sample(s))")
  }

  if (length(all_cells)) {
    cu <- data.table::rbindlist(all_cells, use.names = TRUE, fill = TRUE)
    data.table::fwrite(cu, file.path(opt$outdir, "cells_umap.csv"))
    log_msg("wrote cells_umap.csv (", nrow(cu), " cells)")

    # ---- optional FlowJo hand-off ------------------------------------------
    # Strictly ADDITIVE: everything above has already been written, and the run
    # continues to the group comparison below either way.
    #
    # Sourced rather than duplicated, and called with `cu` already in memory, so
    # the exporter cannot drift from the schema that just produced it and the
    # 13 MB CSV is not written and re-read for no reason.
    #
    # tryCatch, not a bare call: a failure here (unwritable --flowjo-outdir, a
    # flowCore complaint about one odd sample) must not discard a 40-minute run
    # whose real outputs are already on disk. It degrades to a warning.
    if (isTRUE(opt$flowjo_export)) {
      fj_dir <- opt$flowjo_outdir %||% file.path(opt$outdir, "flowjo")
      tryCatch({
        export_flowjo_fcs(cu, fj_dir,
                          concat = !isTRUE(opt$flowjo_no_concat),
                          groups = !isTRUE(opt$flowjo_no_groups),
                          log    = log_msg)
      }, error = function(e) {
        log_msg("WARNING FlowJo export failed: ", conditionMessage(e))
        log_msg("        every default output above is unaffected")
      })
    }
  }

  # ---- between-group abundance comparison -----------------------------------
  # Groups come from the patient table via the sample map, so the grouping is
  # declared in the study metadata rather than parsed out of filenames. Absent a
  # group column the same figure still renders with all samples pooled, which is
  # the frequency overview.
  #
  # Computed here (not inside `if (!is.null(freq))` below) because
  # --absolute-counts wants the same group_of/p_src/grouped and is independent
  # of this run's own gating having produced a freq table at all.
  group_of <- NULL
  gcol <- opt$group_column %||% "cohort"
  if (!is.null(patients) && gcol %in% names(patients) &&
      !is.null(smap) && "patient_id" %in% names(smap)) {
    group_of <- resolve_group_of(patients, smap, gcol)
  } else if (!is.null(smap) && gcol %in% names(smap) &&
             "sample_id" %in% names(smap)) {
    # FALL BACK TO THE SAMPLE MAP. The patient table is optional -- a design
    # with no clinical metadata still has study groups -- and a sample map that
    # names the group column is a complete specification of the grouping on its
    # own. Requiring the patient table meant a run with `cohort` sitting in the
    # sample map reported it "not found in the patient table" and silently
    # skipped the grouped comparison, the compositional tests and the paired
    # test, which is most of the point of the run.
    #
    # The batch-correction step above already reads the group straight from the
    # sample map, so without this the pipeline disagreed with itself about where
    # the group lives.
    g <- as.character(smap[[gcol]])
    keep <- !is.na(g) & nzchar(g) & !is.na(smap$sample_id)
    if (any(keep)) {
      group_of <- setNames(g[keep], smap$sample_id[keep])
      log_msg("  group column '", gcol, "' taken from the sample map (",
              length(unique(group_of)), " group(s))")
    }
  }
  if (is.null(group_of) && !is.null(opt$group_column)) {
    log_msg("NOTE --group-column '", gcol, "' found in neither the patient ",
            "table nor the sample map; grouped comparison skipped. Available ",
            "in the sample map: ",
            paste(setdiff(names(smap), c("file", "sample_id")), collapse = ", "),
            if (!is.null(patients))
              paste0("; in the patient table: ",
                     paste(setdiff(names(patients), "patient_id"), collapse = ", ")))
  }
  # Defined here rather than at STEP 7c, where it used to be, because the design
  # feasibility table and the parametric tests below are the first things to use
  # it and a function is not visible before the line that creates it.
  #
  # Everything it wraps is an addition to a run that has already produced its
  # baseline output, so a failure inside one is logged and stepped over. An
  # addition that can delete the thing it was added to is not an addition.
  .ext_ok <- function(label, expr) {
    tryCatch(expr, error = function(e)
      log_msg("  WARNING ", label, " failed: ", conditionMessage(e),
              ", every other output is unaffected"))
  }

  p_src <- if (identical(opt$p_adjust_display, "BH")) "BH" else "raw"
  grouped <- !is.null(group_of) && length(unique(group_of)) > 1L

  # WHY GROUPING AND TESTING ARE SEPARATE DECISIONS. `grouped` says a group
  # column resolved, and it goes on driving how figures are coloured and split.
  # `group_tests` says the between-group null hypothesis should actually be
  # tested. A cohort can legitimately want the first without the second: with one
  # donor in a group there is no comparison to make, yet the per-group figures
  # are still the point of the run. --no-group-tests separates the two.
  #
  # Diagnostics are deliberately NOT covered by it. Confounding, threshold drift
  # and batch mixing keep their own --no-* flags, because on a design too small
  # to test they become more informative, not less.
  min_group_n <- suppressWarnings(as.integer(opt$min_group_n %||% 3L))
  if (is.na(min_group_n) || min_group_n < 2L) min_group_n <- 3L
  group_tests <- grouped && !isTRUE(opt$no_group_tests)
  if (grouped) {
    .donor_of <- if (!is.null(smap) &&
                     all(c("sample_id", "patient_id") %in% names(smap)))
      setNames(as.character(smap$patient_id), as.character(smap$sample_id)) else NULL
    .ext_ok("design feasibility", {
      write_design_feasibility(group_of, gcol, min_group_n, opt$reference_group,
                               opt$outdir, donor_of = .donor_of,
                               tested = group_tests)
    })
  }
  if (grouped && !group_tests)
    log_msg("--no-group-tests: figures stay grouped by '", gcol,
            "', and no between-group test is run. Diagnostics are unaffected.")

  if (!is.null(freq)) {
    # One pass per marker panel when more than one is present -- the same
    # split fig_umap_overview()/fig_umap_markers() already use for the
    # embedding (sfx). Pooling populations scored under DIFFERENT marker sets
    # and thresholds into one bar chart, with no indication in the figure of
    # which panel a given point came from, is misleading in exactly the way
    # this pipeline's per-sample gating is careful to avoid everywhere else.
    # With a single panel -- the common case -- this is one iteration and every
    # filename below is unchanged from before.
    multi_panel <- length(fpr$panels) > 1L
    panel_names <- if (multi_panel) vapply(fpr$panels, `[[`, "", "name") else NA_character_
    for (pn in panel_names) {
      sfx    <- if (multi_panel) paste0("_", pn) else ""
      plabel <- if (multi_panel) pn else ""
      freq_p <- if (multi_panel) freq[freq$panel == pn, , drop = FALSE] else freq
      fx_p   <- if (multi_panel && !is.null(fx)) fx[fx$panel == pn, , drop = FALSE] else fx
      rt_p   <- if (multi_panel && !is.null(rt))
        rt[!is.na(rt$panel) & rt$panel == pn, , drop = FALSE] else rt

      gstats <- NULL
      if (group_tests) {
        gstats <- stats_group_comparison(freq_p, group_of,
                                         reference = opt$reference_group,
                                         min_n = min_group_n)
        if (!is.null(gstats)) {
          # Two columns appended after every existing one: the typical gate
          # uncertainty behind this population, and how many multiples of it the
          # observed difference is. A ratio below 1 says the groups differ by less
          # than the distance the cut itself moves, which no p-value discloses.
          #
          # Guarded on `unc` rather than left to the function's own NULL handling,
          # so --no-uncertainty writes the table it always wrote. Two columns of
          # NA are not information, and a file whose shape depends on a flag that
          # produced nothing is a worse contract than one that does not change.
          if (!is.null(unc))
            gstats <- tryCatch(annotate_gate_uncertainty(gstats, unc$frequencies),
                               error = function(e) gstats)
          write.csv(gstats, file.path(opt$outdir,
                                      paste0("group_comparison_stats", sfx, ".csv")),
                    row.names = FALSE)
          if (!is.null(unc)) {
            .nres <- sum(is.finite(gstats$difference_over_gate_u) &
                           gstats$difference_over_gate_u < 1 &
                           gstats$significant_raw, na.rm = TRUE)
            if (.nres > 0)
              log_msg("  NOTE ", .nres, " nominally significant result(s) differ ",
                      "by LESS than the gate placement uncertainty behind them. ",
                      "See difference_over_gate_u and threshold_uncertainty.csv.")
          }
          nsig <- sum(gstats$significant_raw, na.rm = TRUE)
          nbh  <- sum(gstats$significant_BH,  na.rm = TRUE)
          log_msg("wrote group_comparison_stats", sfx, ".csv (", nrow(gstats),
                  " tests; ", nsig, " with raw p<0.05, ", nbh,
                  " surviving BH correction)")
          if (nsig > 0 && nbh == 0)
            log_msg("NOTE no comparison survives multiple-testing correction, treat ",
                    "the raw-p hits as hypotheses to confirm, not findings.")
        }

        # The parametric equivalents, in their own files. Written beside the
        # rank tests rather than in place of them, so the primary result does
        # not change and the t-test or ANOVA a journal asks for exists with its
        # assumption checks attached.
        if (!isTRUE(opt$no_parametric)) {
          .ext_ok("parametric tests", {
            pt <- parametric_group_tests(freq_p, group_of,
                                         reference = opt$reference_group,
                                         min_n = min_group_n)
            if (!is.null(pt)) {
              write.csv(pt, file.path(opt$outdir,
                                      paste0("parametric_tests", sfx, ".csv")),
                        row.names = FALSE)
              nok <- sum(pt$assumptions_met, na.rm = TRUE)
              log_msg("wrote parametric_tests", sfx, ".csv (", nrow(pt),
                      " test(s); ", nok, " of ", nrow(pt),
                      " met both normality and equal-variance assumptions)")
              if (nok < nrow(pt))
                log_msg("  where assumptions_met is FALSE the rank test in ",
                        "group_comparison_stats", sfx, ".csv is the defensible one")
            }
            ph <- posthoc_group_tests(freq_p, group_of, min_n = min_group_n)
            if (!is.null(ph)) {
              write.csv(ph, file.path(opt$outdir,
                                      paste0("posthoc_tests", sfx, ".csv")),
                        row.names = FALSE)
              log_msg("wrote posthoc_tests", sfx, ".csv (", nrow(ph),
                      " pairwise comparison(s) across ",
                      length(unique(ph$test)), " method(s); read the one the ",
                      "assumption columns of parametric_tests", sfx,
                      ".csv point to)")
            }
          })
        }
        fig_group_comparison(
          freq_p, file.path(opt$outdir, paste0("group_comparison", sfx, ".png")),
          group_of = group_of, stats = gstats, reference = opt$reference_group,
          p_source = p_src, panel_label = plabel)
        # The same tests as one figure: effect against evidence, every population
        # at once. group_comparison.png is for reading one population, this is for
        # finding which one to read.
        if (!is.null(gstats))
          .ext_ok("group volcano", {
            fig_group_volcano(
              gstats, file.path(opt$outdir,
                                paste0("group_differences", sfx, ".png")),
              title_noun = if (nzchar(plabel))
                paste0("Population abundance, ", plabel) else "Population abundance")
          })
      }
      # Always rendered, pooled and in per-cent -- the cohort-composition view.
      fig_population_frequencies(
        freq_p, file.path(opt$outdir, paste0("population_frequencies", sfx, ".png")),
        panel_label = plabel)

      # ---- functional-marker positivity, same grouping/stats as abundance ---
      if (!is.null(fx_p) && nrow(fx_p)) {
        fxstats <- NULL
        if (group_tests) {
          fxstats <- stats_group_comparison(
            fx_panel_key(fx_p), group_of, reference = opt$reference_group,
            value_col = "pct_positive", min_n = min_group_n)
          if (!is.null(fxstats)) {
            write.csv(fxstats, file.path(opt$outdir,
                                         paste0("functional_markers_stats", sfx, ".csv")),
                      row.names = FALSE)
            log_msg("wrote functional_markers_stats", sfx, ".csv (",
                    nrow(fxstats), " tests)")
          }
        }
        fig_functional_markers(
          fx_p, file.path(opt$outdir, paste0("functional_markers", sfx, ".png")),
          group_of = if (grouped) group_of else NULL, stats = fxstats,
          reference = opt$reference_group, p_source = p_src, panel_label = plabel)
      }

      # ---- derived-ratio comparison, same grouping/stats as abundance ------
      if (!is.null(rt_p) && nrow(rt_p)) {
        rtstats <- NULL
        if (group_tests) {
          rtstats <- stats_group_comparison(rt_p, group_of,
                                            reference = opt$reference_group,
                                            value_col = "value",
                                            min_n = min_group_n)
          if (!is.null(rtstats)) {
            write.csv(rtstats, file.path(opt$outdir,
                                         paste0("population_ratios_stats", sfx, ".csv")),
                      row.names = FALSE)
            log_msg("wrote population_ratios_stats", sfx, ".csv (",
                    nrow(rtstats), " tests)")
          }
        }
        fig_population_ratios(
          rt_p, file.path(opt$outdir, paste0("population_ratios", sfx, ".png")),
          group_of = if (grouped) group_of else NULL, stats = rtstats,
          reference = opt$reference_group, p_source = p_src, panel_label = plabel)
      }
    }
  }

  # ---- external absolute cell counts ----------------------------------------
  # Either from --absolute-counts, or from the sheet's count.<population>
  # columns. Both arrive here as the same long table, so everything below is
  # blind to which supplied it.
  ac_path <- opt[["absolute_counts", exact = TRUE]]
  if (!is.null(ac_path) || !is.null(.sheet$counts)) {
    log_step("STEP 7b - absolute cell counts")
    ac <- if (!is.null(.sheet$counts)) .sheet$counts
          else load_absolute_counts(ac_path, smap, opt$outdir)
    if (!is.null(ac)) {
      run_samples <- names(verdicts)
      not_in_run <- setdiff(unique(ac$sample_id), run_samples)
      if (length(not_in_run)) {
        log_msg("  NOTE --absolute-counts: ", length(not_in_run), " matched sample(s) are",
                " not part of this run (different --dir/--files subset?), dropped: ",
                paste(not_in_run, collapse = ", "))
        ac <- ac[ac$sample_id %in% run_samples, , drop = FALSE]
      }
      missing_from_ac <- setdiff(run_samples, unique(ac$sample_id))
      if (length(missing_from_ac))
        log_msg("  NOTE --absolute-counts: ", length(missing_from_ac), " sample(s) in this",
                " run have no matching row in the counts file: ",
                paste(missing_from_ac, collapse = ", "))

      if (!nrow(ac)) {
        log_msg("NOTE --absolute-counts: no matched rows overlap this run's samples; skipped")
      } else {
        ac$is_control <- vapply(ac$sample_id,
                                function(s) isTRUE(verdicts[[s]]$is_control), logical(1))
        ac$qc_status  <- vapply(ac$sample_id,
                                function(s) verdicts[[s]]$qc_status %||% "pass", character(1))

        write.csv(ac[, c("sample_id", "population", "cells_per_ul", "is_control", "qc_status")],
                  file.path(opt$outdir, "absolute_counts.csv"), row.names = FALSE)
        log_msg("wrote absolute_counts.csv (", length(unique(ac$sample_id)), " sample(s), ",
                length(unique(ac$population)), " population(s))")

        # QC figure FIRST and unconditionally: this data was measured outside this
        # pipeline, so a mismatch or a unit error needs to be visible before the
        # (much more inviting) group-comparison figure gets read as a finding.
        fig_absolute_counts_qc(ac, file.path(opt$outdir, "absolute_counts_qc.png"),
                               group_of = group_of)

        acstats <- NULL
        if (group_tests) {
          acstats <- stats_group_comparison(ac, group_of, reference = opt$reference_group,
                                            value_col = "cells_per_ul",
                                            min_n = min_group_n)
          if (!is.null(acstats)) {
            write.csv(acstats, file.path(opt$outdir, "absolute_counts_stats.csv"),
                      row.names = FALSE)
            log_msg("wrote absolute_counts_stats.csv (", nrow(acstats), " tests)")
          }
        }
        fig_group_comparison(
          ac, file.path(opt$outdir, "absolute_counts.png"),
          group_of = if (grouped) group_of else NULL, stats = acstats,
          reference = opt$reference_group, p_source = p_src,
          value_col = "cells_per_ul", value_label = "cells / \u00b5L",
          value_caveat = paste(
            "directly measured (--absolute-counts), independent of this pipeline's",
            "own gating, see absolute_counts_qc.png first."),
          title_noun = "Absolute cell count (external)")
      }
    }
  }

  # ===========================================================================
  # STEP 7c -- EXTENDED ANALYSES
  # ===========================================================================
  # Placed after every baseline table and figure is on disk, and each item
  # wrapped in tryCatch, so nothing here can cost a completed run its results.
  # That is the same contract the FlowJo export above operates under, and for the
  # same reason: these are additions, and an addition that can delete the thing
  # it was added to is not an addition.
  log_step("STEP 7c - extension analyses")
  # .ext_ok is defined earlier in this function, before its first use in the
  # group-comparison section.

  .cells_all <- if (length(all_cells))
    as.data.frame(data.table::rbindlist(all_cells, use.names = TRUE, fill = TRUE))
  else NULL

  # ---- gate uncertainty figures ---------------------------------------------
  # Drawn here rather than at STEP 4b because the first of them wants the study
  # groups, which are not resolved until the comparison above.
  if (!is.null(unc)) {
    .ext_ok("gate uncertainty figures", {
      fig_frequency_uncertainty(unc$frequencies,
                                file.path(opt$outdir, "frequency_uncertainty.png"),
                                group_of = if (grouped) group_of else NULL)
      fig_uncertainty_budget(unc$budget,
                             file.path(opt$outdir, "uncertainty_budget.png"))
      fig_detection_limits(unc$frequencies,
                           file.path(opt$outdir, "detection_limits.png"))
    })
  }
  # Drawn here rather than at STEP 1b, where it is computed, so that adding this
  # diagnostic leaves every pre-existing figure byte-identical. See the note at
  # STEP 1b.
  if (!is.null(.aqc) && !is.null(.aqc$bins)) {
    .ext_ok("acquisition QC figure", {
      fig_acquisition_qc(.aqc$bins, file.path(opt$outdir, "acquisition_qc.png"),
                         summary = .aqc$summary)
    })
  }

  # ---- differential state: population x marker, tested on SAMPLES -----------
  if (!isTRUE(opt$no_differential_state) && !is.null(mfi) && TRUE) {
    .ext_ok("differential state", {
      ms <- if (group_tests)
        stats_marker_state(mfi, group_of, reference = opt$reference_group,
                           min_n = min_group_n) else NULL
      if (!is.null(ms)) {
        write.csv(ms, file.path(opt$outdir, "marker_state_stats.csv"), row.names = FALSE)
        log_msg("wrote marker_state_stats.csv (", nrow(ms), " tests; ",
                sum(ms$significant_raw, na.rm = TRUE), " with raw p<0.05, ",
                sum(ms$significant_BH, na.rm = TRUE), " surviving BH)")
        log_msg("  NOTE these are SAMPLE-level tests (n = donors). ",
                "subcluster_marker_shifts.csv ranks the same kind of shift over ",
                "POOLED CELLS and carries no p-value, the two answer different ",
                "questions and should not be read as one.")
      }
      fig_marker_state(mfi, file.path(opt$outdir, "marker_state.png"),
                       group_of = if (grouped) group_of else NULL, stats = ms,
                       reference = opt$reference_group, p_source = p_src)
    })
  }

  # ---- phenotype and composition heatmaps ----------------------------------
  if (!isTRUE(opt$no_heatmaps) && TRUE) {
    .ext_ok("phenotype heatmap", {
      if (!is.null(mfi))
        fig_population_marker_heatmap(
          mfi, file.path(opt$outdir, "population_marker_heatmap.png"),
          annotate_expected = expected_positive_markers(spec))
    })
    .ext_ok("cohort composition heatmap", {
      if (!is.null(.cells_all) && grouped)
        fig_cohort_confusion(.cells_all,
                             file.path(opt$outdir, "cohort_composition_heatmap.png"),
                             group_col = gcol)
    })
  }

  # ---- compositional (CLR) re-test of the abundance comparison -------------
  if (!isTRUE(opt$no_compositional) && group_tests && !is.null(freq) &&
      TRUE) {
    .ext_ok("compositional analysis", {
      # The CLR closure is per SAMPLE, and a sample belongs to exactly one panel,
      # so the geometric mean is never taken across panels even in a multi-panel
      # run. What IS pooled is the resulting stats table -- the same convention
      # population_frequencies.csv already follows. Said out loud, because the
      # baseline splits FIGURES by panel for good reason.
      if (length(fpr$panels) > 1L)
        log_msg("  NOTE multi-panel run: CLR is computed within each sample (so ",
                "never across panels), but the stats table below pools panels, ",
                "as population_frequencies.csv does.")
      fclr <- clr_frequencies(freq)
      if (!is.null(fclr)) {
        cstats <- stats_group_comparison(fclr, group_of,
                                         reference = opt$reference_group,
                                         value_col = "clr")
        if (!is.null(cstats)) {
          write.csv(cstats, file.path(opt$outdir, "compositional_clr_stats.csv"),
                    row.names = FALSE)
          log_msg("wrote compositional_clr_stats.csv (", nrow(cstats),
                  " tests on centred log-ratios; ", attr(fclr, "clr_n_zeros"),
                  " zero value(s) replaced at delta = ",
                  signif(attr(fclr, "clr_delta"), 3), ")")
          raw <- try(stats_group_comparison(freq, group_of,
                                            reference = opt$reference_group),
                     silent = TRUE)
          if (!inherits(raw, "try-error")) {
            conc <- compositional_concordance(raw, cstats)
            if (!is.null(conc)) {
              write.csv(conc, file.path(opt$outdir, "compositional_concordance.csv"),
                        row.names = FALSE)
              n_art <- sum(grepl("^raw_only", conc$verdict))
              n_msk <- sum(grepl("^clr_only", conc$verdict))
              log_msg("wrote compositional_concordance.csv (",
                      sum(conc$verdict == "robust_to_composition"),
                      " robust, ", n_art, " raw-only, ", n_msk, " CLR-only)")
              if (n_art > 0)
                log_msg("  NOTE ", n_art, " comparison(s) are significant on raw ",
                        "percentages but NOT on CLR. Percentages are forced to sum ",
                        "to 100, so those may be a consequence of another ",
                        "population moving rather than their own change.")
            }
          }
        }
      }
    })
  }

  # ---- confounding: is a cohort difference really an age or sex difference? -
  if (!isTRUE(opt$no_confounding) && grouped && !is.null(freq) &&
      TRUE) {
    .ext_ok("confounding diagnostic", {
      cvs <- trimws(strsplit(opt$covariates %||% "age,sex", ",")[[1]])
      cvs <- cvs[nzchar(cvs)]
      cf <- stats_confounding(freq, patients, smap, group_of, covariates = cvs)
      if (!is.null(cf)) {
        write.csv(cf, file.path(opt$outdir, "confounding_diagnostics.csv"),
                  row.names = FALSE)
        hi <- sum(grepl("^HIGH", cf$confounder_risk))
        log_msg("wrote confounding_diagnostics.csv (", nrow(cf), " checks; ",
                hi, " flagged HIGH)")
        if (hi > 0)
          log_msg("  NOTE ", hi, " population(s) have a covariate that is BOTH ",
                  "unbalanced across cohorts AND associated with abundance. For ",
                  "those, a cohort difference and a covariate difference are not ",
                  "separable by this design.")
      }
      if (isTRUE(opt$rank_ancova)) {
        ra <- stats_rank_ancova(freq, patients, smap, group_of, covariates = cvs)
        if (!is.null(ra)) {
          write.csv(ra, file.path(opt$outdir, "covariate_adjusted_stats.csv"),
                    row.names = FALSE)
          log_msg("wrote covariate_adjusted_stats.csv (", nrow(ra),
                  " EXPLORATORY rank-ANCOVA fits; ",
                  sum(grepl("^NOT FITTED", ra$status)),
                  " refused for too few residual df)")
        }
      }
    })
  }

  # ---- clinical variables against populations and markers -------------------
  # Read from the SHEET, not from `patients`: the sheet recognises a fixed set of
  # subject attributes and carries everything else as a study column, which is
  # exactly where a severity score or an outcome flag ends up. Taking them from
  # smap is what makes those columns reachable at all.
  .clin_opt <- opt[["clinical_columns", exact = TRUE]]
  if (!is.null(.clin_opt) && nzchar(.clin_opt) && !is.null(smap)) {
    .ext_ok("clinical association", {
      cl <- trimws(strsplit(.clin_opt, ",")[[1]]); cl <- cl[nzchar(cl)]
      .sid <- as.character(smap$sample_id %||% smap$file)
      # A CLINICAL COLUMN CAN LIVE IN EITHER TABLE, so both are searched.
      # read_samplesheet() moves the columns it recognises as subject attributes
      # -- age_years, sex, infection_focus, cohort and the rest -- out of the
      # sample map and into `patients`, keyed by patient_id. Everything it does
      # not recognise stays in the map. Looking in only one of the two made half
      # the sheet's own columns invisible: age_years and infection_focus were
      # reported "not in the sheet" for a sheet that plainly contains them.
      .clin_lookup <- function(cv) {
        if (cv %in% names(smap))
          return(stats::setNames(smap[[cv]], .sid))
        if (!is.null(patients) && cv %in% names(patients) &&
            "patient_id" %in% names(smap)) {
          row <- match(norm_id(smap$patient_id), norm_id(patients$patient_id))
          return(stats::setNames(patients[[cv]][row], .sid))
        }
        NULL
      }
      found <- stats::setNames(lapply(cl, .clin_lookup), cl)
      have <- cl[!vapply(found, is.null, logical(1))]
      miss <- setdiff(cl, have)
      if (length(miss))
        log_msg("  NOTE --clinical-columns in neither the sample map nor the ",
                "patient attributes, skipped: ", paste(miss, collapse = ", "))
      if (length(have)) {
        clin <- found[have]
        # BEFORE the association, and outside the branch that needs it to have
        # produced something. This figure is about the sheet's own columns: it
        # needs no frequency table and no estimable test, and on a cohort too
        # small for the association it is the one thing still worth drawing --
        # whether the variables the study collected say the same thing twice.
        fig_clinical_correlogram(
          clin, file.path(opt$outdir, "clinical_variables_correlation.png"))
        # THE PATIENT MAP MATTERS TO THE RESULT, not just to the labels. A column
        # that is constant within a patient is tested on one value per patient;
        # without this the survivor arm of a cohort sampled three times over is
        # one patient counted three times. See clin_variable_unit().
        .pat_of <- if ("patient_id" %in% names(smap))
          stats::setNames(as.character(smap$patient_id), .sid) else NULL
        ca <- stats_clinical_association(freq, mfi, clin, patient_of = .pat_of)
        if (!is.null(ca)) {
          if (!is.null(ca$populations)) {
            write.csv(ca$populations,
                      file.path(opt$outdir, "clinical_association.csv"),
                      row.names = FALSE)
            log_msg("wrote clinical_association.csv (", nrow(ca$populations),
                    " population x variable test(s), ",
                    sum(ca$populations$significant_BH, na.rm = TRUE),
                    " significant after BH within variable; ",
                    sum(ca$populations$underpowered, na.rm = TRUE),
                    " on fewer than 10 samples)")
            fig_clinical_heatmap(ca$populations, "population",
                                 file.path(opt$outdir, "clinical_association.png"),
                                 title = "Clinical variables against population abundance")
            for (cv in have) {
              .cvsafe <- gsub("[^A-Za-z0-9]+", "_", cv)
              fig_clinical_detail(
                freq, clin[[cv]], cv,
                file.path(opt$outdir, paste0("clinical_", .cvsafe, ".png")))
              # The scatter shows the data, the forest shows how well the effect is
              # pinned down. Both per variable, because an interval is only
              # readable against one axis at a time.
              fig_clinical_forest(
                ca$populations, "population", cv,
                file.path(opt$outdir,
                          paste0("clinical_effects_", .cvsafe, ".png")))
            }
            # One picture of the whole cohort: severity, outcome, timepoint and
            # infection focus as strips over the population matrix.
            fig_clinical_landscape(
              freq, clin, file.path(opt$outdir, "clinical_landscape.png"))
          }
          if (!is.null(ca$markers)) {
            write.csv(ca$markers,
                      file.path(opt$outdir, "clinical_association_markers.csv"),
                      row.names = FALSE)
            log_msg("wrote clinical_association_markers.csv (",
                    nrow(ca$markers), " marker x variable test(s))")
            fig_clinical_heatmap(ca$markers, "marker",
                                 file.path(opt$outdir,
                                           "clinical_association_markers.png"),
                                 title = "Clinical variables against marker intensity")
          }
        } else {
          log_msg("  NOTE clinical association produced no testable pair; every ",
                  "variable was constant, all-missing, or below the minimum ",
                  "sample count")
        }
      }
    })
  }

  # ---- population abundance per acquisition batch ---------------------------
  # Separate from the batch diagnostics, which ask whether the EMBEDDING
  # separates by batch. This asks whether the reported numbers do.
  # opt$batch_column, not `bcol`: that name is assigned inside an earlier block
  # and is not in scope here.
  .bc <- opt[["batch_column", exact = TRUE]]
  if (!is.null(.bc) && !is.null(smap) && .bc %in% names(smap)) {
    .ext_ok("populations by batch", {
      .sidb <- as.character(smap$sample_id %||% smap$file)
      fig_populations_by_batch(
        freq, stats::setNames(as.character(smap[[.bc]]), .sidb),
        file.path(opt$outdir, "populations_by_batch.png"))
    })
  }

  # ---- paired / repeated-measures design ------------------------------------
  if (!is.null(opt$paired_column) && !is.null(freq) && TRUE) {
    .ext_ok("paired comparison", {
      # The pairing and condition columns may live in either metadata table; look
      # in both rather than making the user remember which.
      .col_of <- function(cn) {
        if (!is.null(smap) && cn %in% names(smap) && "sample_id" %in% names(smap))
          return(setNames(as.character(smap[[cn]]), smap$sample_id))
        if (!is.null(patients) && cn %in% names(patients) && !is.null(smap) &&
            all(c("sample_id", "patient_id") %in% names(smap))) {
          r <- match(norm_id(smap$patient_id), norm_id(patients$patient_id))
          return(setNames(as.character(patients[[cn]][r]), smap$sample_id))
        }
        NULL
      }
      pr <- .col_of(opt$paired_column)
      cd <- if (!is.null(opt$condition_column)) .col_of(opt$condition_column) else NULL
      if (is.null(pr)) {
        log_msg("  NOTE --paired-column '", opt$paired_column,
                "' not found in the sample map or patient table; paired test skipped")
      } else if (is.null(cd)) {
        log_msg("  NOTE --paired-column needs --condition-column (which condition ",
                "each member of a pair is); paired test skipped")
      } else {
        ps <- stats_paired_comparison(freq, pr, cd)
        if (!is.null(ps)) {
          write.csv(ps, file.path(opt$outdir, "paired_comparison_stats.csv"),
                    row.names = FALSE)
          log_msg("wrote paired_comparison_stats.csv (", nrow(ps), " tests, ",
                  max(ps$n_pairs), " complete pairs; ",
                  sum(ps$n_incomplete_pairs_dropped), " incomplete pair-slots dropped)")
        } else {
          log_msg("  NOTE paired test produced nothing, check that every pairing ",
                  "unit has at least two conditions and at least 3 complete pairs exist")
        }
        # ---- trajectories, one line per patient across the conditions --------
        # Drawn whether or not the test above produced anything: the test needs
        # complete pairs and refuses without them, and a cohort with incomplete
        # follow-up is exactly the one where seeing the lines matters most.
        #
        # The colour is the first TWO-LEVEL clinical column, which is what an
        # outcome flag looks like. A numeric score cannot colour a line usefully
        # (twelve unique values, twelve colours) and a three-level variable makes
        # three medians out of a dozen patients, so neither is used here.
        .oc <- NULL; .ocn <- "outcome"
        .clin_names <- if (!is.null(.clin_opt) && nzchar(.clin_opt))
          trimws(strsplit(.clin_opt, ",")[[1]]) else character(0)
        for (cn in .clin_names) {
          v <- .col_of(cn)
          if (is.null(v)) next
          u <- unique(stats::na.omit(trimws(as.character(v))))
          if (length(u) == 2L) { .oc <- v; .ocn <- cn; break }
        }
        fig_clinical_trajectory(
          freq, cd, pr, file.path(opt$outdir, "population_trajectories.png"),
          outcome_of = .oc, outcome_name = .ocn)
      }
    })
  }

  # ---- threshold drift: does the GATE move with cohort? ---------------------
  if (!isTRUE(opt$no_threshold_drift) && grouped && !is.null(thr_all) &&
      TRUE) {
    .ext_ok("threshold drift", {
      td <- stats_threshold_drift(thr_all, group_of, spec = spec)
      if (!is.null(td)) {
        write.csv(td, file.path(opt$outdir, "threshold_drift_stats.csv"),
                  row.names = FALSE)
        nfl <- sum(grepl("^FLAGGED", td$drift_flag))
        log_msg("wrote threshold_drift_stats.csv (", nrow(td), " markers; ",
                nfl, " flagged)")
        if (nfl > 0)
          log_msg("  NOTE ", nfl, " marker(s) have per-sample thresholds that ",
                  "differ systematically by cohort. For the populations they ",
                  "define, part of any abundance difference is a difference in ",
                  "the DEFINITION, not in the biology.")
      }
      fig_threshold_drift(thr_all, file.path(opt$outdir, "threshold_drift.png"),
                          group_of = group_of, stats = td,
                          reference = opt$reference_group)
    })
  }

  # ---- batch-effect diagnostic ---------------------------------------------
  if (!is.null(.cells_all) && TRUE) {
    .ext_ok("batch diagnostic", {
      bcol <- opt$batch_column
      bcells <- .cells_all
      if (!is.null(bcol) && !bcol %in% names(bcells)) {
        # Not already merged onto the cells table -- pull it from the metadata.
        v <- NULL
        if (!is.null(smap) && bcol %in% names(smap) && "sample_id" %in% names(smap))
          v <- setNames(as.character(smap[[bcol]]), smap$sample_id)
        else if (!is.null(patients) && bcol %in% names(patients) && !is.null(smap) &&
                 all(c("sample_id", "patient_id") %in% names(smap)))
          v <- setNames(as.character(patients[[bcol]][
            match(norm_id(smap$patient_id), norm_id(patients$patient_id))]),
            smap$sample_id)
        if (!is.null(v)) bcells[[bcol]] <- unname(v[bcells$sample_id])
        else {
          log_msg("  NOTE --batch-column '", bcol, "' not found anywhere; ",
                  "falling back to the FCS acquisition date")
          bcol <- NULL
        }
      }
      if (is.null(bcol)) {
        # Free batch variable: the $DATE keyword every acquisition writes. Used
        # only when it actually varies -- one date is one batch, and a diagnostic
        # of a single batch is a figure that can only say "yes, they match".
        ad <- vapply(names(reads), function(s) {
          k <- reads[[s]]$keywords
          v <- k[["$DATE"]] %||% k[["DATE"]] %||% NA_character_
          as.character(v)[1]
        }, character(1))
        if (length(unique(stats::na.omit(ad))) > 1L) {
          bcells$acquisition_date <- unname(ad[bcells$sample_id])
          bcol <- "acquisition_date"
          log_msg("  batch variable: FCS $DATE keyword (",
                  length(unique(stats::na.omit(ad))), " distinct dates)")
        }
      }
      if (!is.null(bcol)) {
        br <- batch_mixing_report(bcells, bcol, group_col = gcol, seed = opt$seed)
        if (!is.null(br)) {
          write.csv(br$summary, file.path(opt$outdir, "batch_mixing_stats.csv"),
                    row.names = FALSE)
          if (!is.null(br$confounding))
            write.csv(br$confounding,
                      file.path(opt$outdir, "batch_group_confounding.csv"),
                      row.names = FALSE)
          log_msg("wrote batch_mixing_stats.csv, ", br$summary$verdict)
          if (!is.null(br$confounding) && br$confounding$cramers_v >= 0.5)
            log_msg("  NOTE batch and ", gcol, " overlap strongly (Cramer's V = ",
                    br$confounding$cramers_v, "). No batch correction can ",
                    "separate them, which is why this pipeline reports the ",
                    "effect rather than removing it.")
          fig_batch_diagnostic(bcells, br,
                               file.path(opt$outdir, "batch_diagnostic.png"), bcol)
        }

        # WHICH channel moved, not merely whether the batches mix. iLISI above is
        # a property of the embedding and names nothing anyone can act on; a
        # flagged marker names a reagent lot or a detector.
        #
        # Two views, because they fail differently. The threshold test reuses
        # stats_threshold_drift() with batch as the grouping variable -- it was
        # always generic in that argument -- and sees drift that moves the CUT.
        # The distributional test sees a marker that changed its spread or lost
        # its separation while the density minimum between the modes stayed put,
        # which the first cannot detect at all.
        mbd <- marker_batch_drift(bcells, bcol, seed = opt$seed %||% 42L)
        if (!is.null(mbd)) {
          write.csv(mbd, file.path(opt$outdir, "marker_batch_drift.csv"),
                    row.names = FALSE)
          nfl <- sum(mbd$verdict == "differs between batches")
          log_msg("wrote marker_batch_drift.csv (", nfl, " of ", nrow(mbd),
                  " marker(s) differ between batches by at least half their own ",
                  "spread", if (nfl) paste0("; worst: ", mbd$marker[1], ", ",
                                            mbd$emd_over_mad[1], " MAD") else "", ")")
        }
        bmap <- setNames(as.character(bcells[[bcol]]), bcells$sample_id)
        bmap <- bmap[!duplicated(names(bmap))]
        btd <- stats_threshold_drift(thr_all, bmap, spec)
        if (!is.null(btd)) {
          write.csv(btd, file.path(opt$outdir, "threshold_batch_drift.csv"),
                    row.names = FALSE)
          log_msg("wrote threshold_batch_drift.csv (the same test as ",
                  "threshold_drift_stats.csv, grouped by ", bcol,
                  " instead of by study group)")
        }
      }
    })
  }

  # ---- unsupervised clustering vs the gate spec -----------------------------
  if (isTRUE(opt$unsupervised) && TRUE) {
    for (pn in names(embeddings)) {
      .ext_ok(paste0("unsupervised clustering (", pn, ")"), {
        sfx2 <- if (length(embeddings) > 1L) paste0("_", pn) else ""
        ec <- embeddings[[pn]]$cells
        uc <- run_unsupervised_clusters(ec, embeddings[[pn]]$features,
                                        n_clusters = opt$cluster_k,
                                        grid = opt$cluster_grid, seed = opt$seed)
        if (!is.null(uc)) {
          out <- data.frame(sample_id = ec$sample_id,
                            event_index = ec$event_index,
                            population_label = ec$population_label,
                            cluster = uc$cluster, som_node = uc$node,
                            stringsAsFactors = FALSE)
          write.csv(out, file.path(opt$outdir,
                    paste0("unsupervised_clusters", sfx2, ".csv")), row.names = FALSE)
          ag <- cluster_gate_agreement(ec, uc$cluster)
          if (!is.null(ag)) {
            write.csv(ag$per_cluster, file.path(opt$outdir,
                      paste0("cluster_gate_agreement_clusters", sfx2, ".csv")),
                      row.names = FALSE)
            write.csv(ag$per_population, file.path(opt$outdir,
                      paste0("cluster_gate_agreement_populations", sfx2, ".csv")),
                      row.names = FALSE)
            nun <- sum(grepl("^UNDESCRIBED", ag$per_cluster$interpretation))
            nsus <- sum(grepl("^SUSPECT", ag$per_population$interpretation))
            log_msg("wrote cluster_gate_agreement_*", sfx2, ".csv (", nun,
                    " undescribed cluster(s), ", nsus,
                    " population(s) with a suspect threshold)")
            if (nsus > 0)
              log_msg("  NOTE a population flagged SUSPECT THRESHOLD has cells ",
                      "that cluster together phenotypically but fail its Boolean ",
                      "gate. That is a gating problem, not a biological result.")
          }
          fig_unsupervised_clusters(ec, uc$cluster, file.path(opt$outdir,
                                    paste0("unsupervised_clusters", sfx2, ".png")),
                                    agreement = ag, panel_label = pn)

          # ---- learned gating strategies for undescribed clusters ------------
          # Runs only on request. It answers the question the agreement table
          # RAISES but cannot settle -- what IS that cluster -- by fitting a
          # gate to it. Nothing it produces re-enters the analysis: the scored
          # populations, frequencies and every tested result are untouched
          # whether this runs or not. It proposes; a human disposes.
          if (isTRUE(opt$explain_clusters) && !is.null(ag)) {
            .ext_ok(paste0("learned gate proposals (", pn, ")"), {
              ex <- explain_unmatched_clusters(
                ec, uc$cluster, embeddings[[pn]]$features, ag,
                max_clusters = opt$explain_max_clusters %||% 4L,
                max_depth    = opt$explain_max_depth %||% 4L)
              if (!is.null(ex)) {
                write.csv(ex$summary, file.path(opt$outdir,
                          paste0("cluster_gate_proposals", sfx2, ".csv")),
                          row.names = FALSE)
                write.csv(ex$polygons, file.path(opt$outdir,
                          paste0("cluster_gate_polygons", sfx2, ".csv")),
                          row.names = FALSE)
                if (isTRUE(opt$export_gates))
                  .ext_ok(paste0("gate export (", pn, ")"), {
                    .tr <- transforms[[pn]]
                    write_gating_ml(ex$polygons,
                      file.path(opt$outdir,
                                paste0("cluster_gates", sfx2, ".gatingml.xml")),
                      transform = .tr, id_col = "cluster",
                      x_col = "x_asinh", y_col = "y_asinh")
                    .lin <- polygons_linear_table(ex$polygons, .tr,
                              id_col = "cluster", x_col = "x_asinh",
                              y_col = "y_asinh")
                    if (!is.null(.lin))
                      write.csv(.lin, file.path(opt$outdir,
                        paste0("cluster_gate_polygons_linear", sfx2, ".csv")),
                        row.names = FALSE)
                  })
                log_msg("wrote cluster_gate_proposals", sfx2, ".csv and ",
                        "cluster_gate_polygons", sfx2, ".csv (",
                        length(ex$strategies), " strategy/strategies)")
                log_msg("  NOTE a proposed gate is DESCRIPTIVE. It says where ",
                        "these cells sit in marker space, not that they are a ",
                        "real population. Metrics are held-out, but the cluster ",
                        "itself came from the same cells.")
                fx <- as.matrix(ec[, intersect(embeddings[[pn]]$features,
                                               names(ec)), drop = FALSE])
                for (kk in names(ex$strategies)) {
                  fig_gate_strategy(
                    fx, as.integer(uc$cluster == as.integer(kk)),
                    ex$strategies[[kk]],
                    file.path(opt$outdir, paste0("cluster_gate_strategy_",
                                                 kk, sfx2, ".png")),
                    label = paste0("cluster ", kk))
                }
              }
            })
          }
        }
      })
    }
  }

  # ---- labels supplied by another tool --------------------------------------
  # Runs on the embedded cells, so it sits after the clustering block and needs
  # nothing from it. What it produces is a gate for a population this package did
  # not define, measured on donors it was not fitted to. Descriptive throughout:
  # no scored population, frequency or test is touched.
  if (!is.null(opt$external_labels)) {
    .ext_ok("external label gates", {
      .lab <- read_external_labels(opt$external_labels)
      for (pn in names(embeddings)) {
        sfx3 <- if (length(embeddings) > 1L) paste0("_", pn) else ""
        .ec <- join_external_labels(embeddings[[pn]]$cells, .lab)
        if (is.null(.ec)) next
        .xl <- explain_external_labels(
          .ec, embeddings[[pn]]$features, seed = opt$seed %||% 42L,
          max_labels = opt$external_max_labels %||% 6L,
          max_donors = opt$transfer_max_donors %||% 8L,
          transfer_max_cells = opt$transfer_max_cells %||% 20000L,
          max_depth = opt$explain_max_depth %||% 4L)
        if (is.null(.xl)) next
        write.csv(.xl$summary, file.path(opt$outdir,
                  paste0("external_label_gates", sfx3, ".csv")), row.names = FALSE)
        if (!is.null(.xl$polygons))
          write.csv(.xl$polygons, file.path(opt$outdir,
                    paste0("external_label_polygons", sfx3, ".csv")),
                    row.names = FALSE)
        if (!is.null(.xl$transfer)) {
          write.csv(.xl$transfer, file.path(opt$outdir,
                    paste0("gate_transferability", sfx3, ".csv")), row.names = FALSE)
          write.csv(.xl$transfer_summary, file.path(opt$outdir,
                    paste0("gate_transferability_summary", sfx3, ".csv")),
                    row.names = FALSE)
          .weak <- .xl$transfer_summary[
            is.finite(.xl$transfer_summary$f1_min) &
              .xl$transfer_summary$f1_min < 0.5, , drop = FALSE]
          if (nrow(.weak))
            log_msg("  NOTE ", nrow(.weak), " gate(s) score below F1 0.5 on at ",
                    "least one held-out donor, so they describe the donors they ",
                    "were fitted to rather than the population: ",
                    paste(.weak$label, collapse = ", "))
        }
        for (kk in names(.xl$strategies))
          fig_gate_strategy(
            as.matrix(.ec[, intersect(embeddings[[pn]]$features, names(.ec)),
                          drop = FALSE]),
            as.integer(!is.na(.ec$external_label) & .ec$external_label == kk),
            .xl$strategies[[kk]],
            file.path(opt$outdir, paste0("external_label_strategy_",
                                         gsub("[^A-Za-z0-9_.-]", "_", kk),
                                         sfx3, ".png")),
            label = kk)
        if (isTRUE(opt$export_gates) && !is.null(.xl$polygons))
          .ext_ok(paste0("gate export (", pn, ", external labels)"), {
            .tr <- transforms[[pn]]
            write_gating_ml(.xl$polygons, file.path(opt$outdir,
                            paste0("external_label_gates", sfx3, ".gatingml.xml")),
                            transform = .tr)
            .lin <- polygons_linear_table(.xl$polygons, .tr)
            if (!is.null(.lin))
              write.csv(.lin, file.path(opt$outdir,
                        paste0("external_label_polygons_linear", sfx3, ".csv")),
                        row.names = FALSE)
          })
      }
    })
  }

  # ---- conformance against an accepted baseline -----------------------------
  # The second reference. thresholds_used.csv compares each sample against its
  # peers in this run and catches one bad tube; this compares the run against a
  # run that was accepted earlier and catches the case where everything moved
  # together, which the peer check cannot see because the peers moved with it.
  .drift_failed <- FALSE
  if (!is.null(opt$baseline)) {
    # The verdict is taken from what the block RETURNS rather than assigned from
    # inside it, so it does not depend on which environment .ext_ok's promise
    # happens to evaluate in.
    .drift_failed <- isTRUE(.ext_ok("specification conformance", {
      .bl <- read_spec_baseline(opt$baseline)
      .cf <- specification_conformance(thr_all, freq, spec, .bl,
                                       transform = opt$transform %||% "arcsinh")
      if (!is.null(.cf$markers))
        write.csv(.cf$markers, file.path(opt$outdir,
                  "specification_conformance.csv"), row.names = FALSE)
      if (!is.null(.cf$populations))
        write.csv(.cf$populations, file.path(opt$outdir,
                  "specification_conformance_populations.csv"), row.names = FALSE)
      if (!is.null(.cf$spec_changes))
        write.csv(.cf$spec_changes, file.path(opt$outdir,
                  "specification_changes.csv"), row.names = FALSE)
      log_msg("conformance against baseline of ", .cf$summary$baseline_created,
              ": ", .cf$summary$n_fail, " failure(s), ",
              .cf$summary$n_qualify, " qualified")
      if (isTRUE(.cf$summary$transform_changed))
        log_msg("  NOTE the transform differs from the baseline's, so thresholds ",
                "are on a different scale and no marker comparison is meaningful. ",
                "Re-baseline, or re-run with the baseline's transform.")
      if (isTRUE(.cf$summary$spec_changed))
        log_msg("  NOTE the population specification changed since the baseline; ",
                "see specification_changes.csv. A redefined population is a ",
                "different measurement, not a drifting one.")
      if (.cf$summary$n_fail > 0)
        log_msg("  WARNING this run no longer places its cuts where the baseline ",
                "did. Frequencies from the two runs are not the same measurement ",
                "and should not be pooled until someone has looked.")
      .cf$summary$n_fail > 0
    }))
  }
  if (!is.null(opt$write_baseline)) {
    .ext_ok("baseline", {
      write_spec_baseline(opt$write_baseline, thr_all, freq, spec, opt, cofactors)
      log_msg("wrote baseline ", opt$write_baseline,
              " (compare a later run against it with --baseline)")
    })
  }

  # ---- MIFlowCyt report -----------------------------------------------------
  # Written by default, because a reporting standard nobody remembers to ask for
  # is one nobody files. Costs nothing: the instrument section is keywords this
  # run already holds and discards, and the analysis section is the run itself.
  if (!isTRUE(opt$no_miflowcyt)) {
    .ext_ok("MIFlowCyt report", {
      write_miflowcyt(file.path(opt$outdir, "miflowcyt.md"), reads, opt = opt,
                      spec = spec, transforms = transforms, fpr = fpr,
                      thresholds = thr_all)
    })
  }

  # ---- run report ------------------------------------------------------------
  # Written last so it can index everything the run actually produced, and
  # written by default because the reading order is the part of this package a
  # directory listing cannot convey.
  if (!isTRUE(opt$no_report)) {
    # ---- the methods catalogue, and the evidence for it ---------------------
    # Immunophenotyping papers report t-tests and ANOVA; this run reports rank
    # tests. statistical_methods.csv states every commonly reported method and
    # what was done about it, and normality_tests.csv carries the Shapiro-Wilk
    # and Brown-Forsythe results that justify the choice, so a reader can check
    # it rather than take it on trust. Both are written on every run, including
    # ungrouped ones, where the catalogue records why no test ran at all.
    .ext_ok("statistical methods catalogue", {
      ngrp <- if (!is.null(group_of))
        length(unique(stats::na.omit(unname(group_of)))) else NA_integer_
      ntest <- if (exists("gstats", inherits = FALSE) && !is.null(gstats))
        nrow(gstats) else NA_integer_
      smt <- statistical_methods_table(n_groups = ngrp, paired = FALSE,
                                       n_tests = ntest)
      write.csv(smt, file.path(opt$outdir, "statistical_methods.csv"),
                row.names = FALSE)
      log_msg("wrote statistical_methods.csv (", nrow(smt),
              " methods; ", sum(smt$role == "used"), " used)")
    })
    .ext_ok("normality and variance diagnostics", {
      if (!is.null(freq) && !is.null(group_of)) {
        nt <- normality_report(freq, group_of)
        if (!is.null(nt) && nrow(nt)) {
          write.csv(nt, file.path(opt$outdir, "normality_tests.csv"),
                    row.names = FALSE)
          nnorm <- sum(grepl("^departs", nt$interpretation))
          log_msg("wrote normality_tests.csv (", nrow(nt), " tests; ", nnorm,
                  " departing from normal). Shapiro-Wilk at these group sizes ",
                  "has little power, so a non-significant result is not ",
                  "evidence of normality -- which is the argument for the ",
                  "rank tests used above.")
        }
      }
    })

    # ---- explore mode -------------------------------------------------------
    # Everything below writes into <outdir>/explore/ and touches no existing
    # output. Wrapped like every other extension, so a failure here cannot take
    # down a run whose declared deliverables are already on disk.
    if (isTRUE(opt$explore)) {
      .ext_ok("explore mode", {
        # --maybe-learn is the only thing that lets the two analyses see each
        # other. Without it explore gets the transform and the group labels and
        # nothing else: the transform because it is a property of the panel
        # derived from the data rather than anything the specification decided,
        # and the group labels because they come from the sample sheet, not from
        # the analysis. Thresholds, QC verdicts, scored populations and the
        # confounding verdict are all products of the declared run, so they are
        # withheld and explore names its clusters from pooled medians instead --
        # which is what a blind unsupervised run has.
        .learn <- isTRUE(opt$maybe_learn)
        .conf <- NULL
        if (.learn) {
          .cpath <- file.path(opt$outdir, "batch_group_confounding.csv")
          if (file.exists(.cpath))
            .conf <- utils::read.csv(.cpath, stringsAsFactors = FALSE)
        }
        log_msg("explore: ", if (.learn)
          "--maybe-learn is set, the declared run lends its thresholds, QC verdicts and confounding verdict"
          else "isolated from the declared analysis; pass --maybe-learn to link them")
        .ex <- run_explore(reads, fpr, opt, opt$outdir, transforms = transforms,
                           gates = if (.learn) gates else NULL,
                           verdicts = if (.learn) verdicts else NULL,
                           pops = if (.learn) pops else NULL,
                           group_of = group_of, confounding = .conf,
                           file_paths = fcs)

        # The other direction, and the only thing explore is allowed to add to
        # the declared deliverables. One file, written only under the flag, so a
        # run without it produces byte-identical output to a run with no
        # --explore at all.
        if (.learn && !is.null(.ex$gaps) && nrow(.ex$gaps)) {
          write.csv(.ex$gaps, file.path(opt$outdir, "spec_gaps.csv"),
                    row.names = FALSE)
          nlump <- sum(.ex$gaps$issue == "population spans several clusters")
          nmiss <- sum(.ex$gaps$issue == "cluster no declared population covers")
          log_msg("wrote spec_gaps.csv (", nlump, " population(s) spanning ",
                  "several clusters, ", nmiss, " cluster(s) nothing covers). ",
                  "This is what unsupervised clustering says about the ",
                  "specification, and it exists only because --maybe-learn is set.")
        }
      })
    }

    .ext_ok("run report", {
      write_run_report(opt$outdir, opt = opt, verdicts = verdicts)
    })
  }

  # ---- session state --------------------------------------------------------
  if (!isTRUE(opt$no_session)) {
    log_step("STEP 8 - saving session state")
    state <- list(reads = reads, gates = gates, verdicts = verdicts,
                  panels = fpr, cofactors = cofactors, populations = pops,
                  embeddings = embeddings, patients = patients,
                  tables = list(gate_counts = gate_counts, staining_qc = qc,
                                frequencies = freq, mfi = mfi,
                                thresholds = do.call(rbind, thr_rows),
                                unavailable = unav),
                  recon = recon, options = opt, spec = spec, blocks = blocks)
    save_session(file.path(opt$outdir, "session_state.RData"), state,
                 keep_exprs = isTRUE(opt$keep_exprs))
  }

  log_step("DONE - outputs in ", normalizePath(opt$outdir))
  log_msg("INSPECT recon_diagnostics.png AND gating_qc.png BEFORE using any numbers.")
  # Flips the on.exit manifest writer from "failed" to "completed". Set here, at
  # the last statement, so it can only be TRUE if everything above actually ran.
  .finished_ok <- TRUE
  # --fail-on-drift raises AFTER the run is marked completed, because it is not a
  # run failure: every output was produced and every one of them is on disk. It
  # is a verdict about whether this cohort is comparable to the baseline, raised
  # as an error only so a scheduled job can act on it.
  if (isTRUE(opt$fail_on_drift) && isTRUE(.drift_failed))
    stop("specification conformance failed against ", opt$baseline,
         ": see specification_conformance.csv. Every output was written; ",
         "this is --fail-on-drift reporting the verdict as an exit code.",
         call. = FALSE)
  invisible(TRUE)
}

