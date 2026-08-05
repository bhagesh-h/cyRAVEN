# =============================================================================
# cyRAVEN -- turn the embedding into FlowJo-native parameters
# =============================================================================
#
# DUAL-USE: exported as a function and also runnable as a stand-alone script.
#
#   export_flowjo_fcs(cells, "output/flowjo")
#   Rscript -e 'cyRAVEN::export_flowjo_fcs(...)'
#
# The CLI block at the bottom is guarded by `sys.nframe() == 0L`, true only when
# Rscript runs this file directly and false when the package loads it. That is
# what lets run_cyraven(flowjo_export = TRUE) call the function in-process, with
# no second container, no re-read of a CSV it already has in memory, and no
# second copy of this logic to keep in step.
#
# WHY IT EXISTS
#   The pipeline produces a STATIC umap_overview.png. To interrogate the
#   embedding interactively -- lasso a cluster, back-gate it, overlay a marker --
#   the coordinates have to live inside FlowJo as real parameters.
#
#   FlowJo's command line CANNOT compute a UMAP. Its documented operations are
#   loading, compensation, gating, statistics, batching layouts/tables, and
#   export; UMAP is a plugin, and plugins are not among them. So the embedding
#   must be computed outside FlowJo and carried in.
#
# WHY EXTRA FCS PARAMETERS RATHER THAN CLR
#   FlowJo's CLR import (`-importCLR`) is the other candidate route, and it is
#   the wrong one here for two documented reasons:
#
#     1. CLR values are PROBABILITIES -- "a value between 0 and 1" -- and are
#        renormalised to a 4096 range on import. UMAP coordinates are unbounded
#        reals including negatives. Squeezing them into \code{[0,1]} costs precision
#        and makes the axes lie about their own units.
#     2. "The number of data rows in the CLR/CSV file must be equal to the
#        number of events in the sample file." The pipeline subsamples
#        (--cells-per-sample), so a CLR built from cells_umap.csv would have
#        far fewer rows than the source FCS and would be rejected.
#
#   Writing UMAP-1 / UMAP-2 as ordinary FCS parameters sidesteps both. FlowJo
#   treats them exactly like acquired channels: plottable, gateable, and
#   available to statistics -- which is the whole point.
#
# WHAT IT WRITES (into outdir)
#   <sample>.fcs           one per sample: markers + UMAP-1 + UMAP-2 +
#                          Population + SampleID + CohortID
#   _ALL_SAMPLES.fcs       every cell concatenated -- open THIS one first; it is
#                          the whole embedding in a single plot
#   _GROUP_<cohort>.fcs    one per study group, that group's cells concatenated
#   population_codes.csv   Population / CohortID / SampleID code -> label legend
#   samples.txt            newline-delimited sample list, the input FlowJo's CLI
#                          accepts as a first argument to build a workspace with
#                          no template
#
# WHY A PER-GROUP FILE AND NOT JUST A CohortID GATE: the cohort comparison is the
# question this study exists to ask, and in the GUI it is made by overlaying the
# groups on one axis. Doing that from _ALL_SAMPLES means re-deriving a range gate
# on the CohortID code channel every session, in every layout, before any overlay
# can be built. Shipping each group as its own openable sample makes the overlay a
# drag-and-drop: the three cohorts land in the Layout Editor as three series with
# no gating at all. The per-sample files stay for within-group work, and
# _ALL_SAMPLES stays for the pooled embedding -- these are additive, not a
# replacement for either.
# =============================================================================

# Columns that are metadata rather than measured signal. Everything NOT in this
# list is treated as a marker channel, so a panel change needs no edit here.
FLOWJO_META_COLS <- c("patient_id", "sample_id", "population_label", "panel",
                      "event_index", "umap_1", "umap_2", "timepoint", "sex",
                      "age_years", "cohort")

#' Write UMAP-annotated FCS files for FlowJo
#'
#' @param cells  data.frame/data.table of per-cell rows, or a path to
#'   cells_umap.csv. Requires sample_id, umap_1, umap_2, population_label.
#' @param outdir directory to create and write into.
#' @param concat also write the concatenated _ALL_SAMPLES.fcs.
#' @param groups also write one concatenated _GROUP_<level>.fcs per study group.
#' @param group_col column holding the study group. Defaults to `cohort`, the
#'   name the pipeline gives the grouping variable in cells_umap.csv; absent that
#'   column the group files are skipped with a note rather than failing, since
#'   an ungrouped run is legitimate.
#' @param log    message sink; the main pipeline passes its own log_msg so the
#'   export's progress lands in the same run log as everything else.
#' @return character vector of written FCS paths, invisibly.
#' @export
export_flowjo_fcs <- function(cells, outdir, concat = TRUE, groups = TRUE,
                              group_col = "cohort", log = message) {

  for (p in c("flowCore")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("required package not installed: ", p,
           "\n  run: BiocManager::install(\"", p, "\")", call. = FALSE)
    }
  }

  if (is.character(cells) && length(cells) == 1L) {
    if (!file.exists(cells)) {
      stop("cells file not found: ", cells,
           "\nRun the main pipeline first \u2014 this consumes its output.", call. = FALSE)
    }
    cells <- utils::read.csv(cells, check.names = FALSE, stringsAsFactors = FALSE)
  }
  cells <- as.data.frame(cells, stringsAsFactors = FALSE)

  if (!nrow(cells)) stop("no cells to export", call. = FALSE)

  required <- c("sample_id", "umap_1", "umap_2", "population_label")
  missing  <- setdiff(required, names(cells))
  if (length(missing)) {
    stop("cells is missing required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  markers <- setdiff(names(cells), FLOWJO_META_COLS)
  # A non-numeric column here would silently become NA in the matrix, so drop
  # anything unmeasurable rather than exporting a channel of NaNs.
  markers <- markers[vapply(cells[markers], is.numeric, logical(1))]
  if (!length(markers)) stop("no numeric marker columns found", call. = FALSE)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # FCS carries numbers, not strings. Every categorical becomes an integer code
  # plus a legend CSV, because a factor level silently renumbering between runs
  # would make two exports look comparable when they are not. Levels are sorted
  # so the mapping is deterministic.
  # NA/blank is folded into an explicit "unassigned" LEVEL rather than left to
  # become an NA code. sort() drops NA, so match() against its output returns NA
  # for those rows; that NA reaches the FCS matrix, makes the channel's $PnR
  # unwritable, and flowCore dies with "missing value where TRUE/FALSE needed" --
  # an error that names nothing relevant. This is reachable in a normal run: the
  # pipeline deliberately tolerates a sample whose patient_id has no patient-table
  # row ("covariates will be NA for these samples"), and every such cell arrives
  # here with cohort = NA. Giving them a real level keeps them in the export, in
  # the legend, and visible as their own group instead of crashing the write.
  code_of <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(trimws(x))] <- "unassigned"
    lv <- sort(unique(x))
    list(codes = match(x, lv), levels = lv)
  }

  pop    <- code_of(cells$population_label)
  samp   <- code_of(cells$sample_id)
  cohort <- if ("cohort" %in% names(cells)) code_of(cells$cohort)
            else list(codes = rep(1L, nrow(cells)), levels = "unknown")

  legend <- rbind(
    data.frame(channel = "Population", code = seq_along(pop$levels),    label = pop$levels),
    data.frame(channel = "SampleID",   code = seq_along(samp$levels),   label = samp$levels),
    data.frame(channel = "CohortID",   code = seq_along(cohort$levels), label = cohort$levels),
    stringsAsFactors = FALSE
  )
  utils::write.csv(legend, file.path(outdir, "population_codes.csv"), row.names = FALSE)

  mat <- cbind(
    as.matrix(cells[, markers, drop = FALSE]),
    `UMAP-1`   = cells$umap_1,
    `UMAP-2`   = cells$umap_2,
    Population = pop$codes,
    SampleID   = samp$codes,
    CohortID   = cohort$codes
  )
  storage.mode(mat) <- "numeric"

  # FCS has no representation for NA/Inf, and leaving one in place fails the same
  # way an NA code does -- an unwritable $PnR and a flowCore error that names
  # nothing useful. This is reachable whenever a run spans more than one panel:
  # the pipeline concatenates panels with rbindlist(fill = TRUE), so a marker
  # present in one panel and absent from another arrives as NA for the other
  # panel's cells. Name the affected channels and the count, then zero-fill, so
  # one absent marker costs a channel's honesty rather than the entire export.
  bad <- !is.finite(mat)
  if (any(bad)) {
    hit <- colnames(mat)[apply(bad, 2L, any)]
    log("  WARNING ", sum(bad), " non-finite value(s) written as 0, in channel(s): ",
        paste(hit, collapse = ", "))
    log("          (typically a marker absent from one panel of a multi-panel run;",
        " those cells are NOT ", paste(hit, collapse = "/"), "-negative)")
    mat[bad] <- 0
  }

  # $PnR is the channel's declared range. flowCore defaults it to the data max,
  # which for a categorical code channel would make FlowJo autoscale the axis to
  # 13 and render every population on top of the axis label. Declaring ranges
  # explicitly keeps the code channels legible and gives the UMAP axes a little
  # padding so points do not sit on the plot border.
  build_frame <- function(sub_mat, sample_name) {
    mins <- apply(sub_mat, 2, min, na.rm = TRUE)
    maxs <- apply(sub_mat, 2, max, na.rm = TRUE)
    # apply(min) over an all-NA or zero-row column returns +/-Inf, which is as
    # unwritable as NA. Cannot happen once the non-finite sweep above has run,
    # but $PnR is not worth betting a whole run on.
    mins[!is.finite(mins)] <- 0
    maxs[!is.finite(maxs)] <- 1
    pad  <- (maxs - mins) * 0.02
    pad[!is.finite(pad) | pad == 0] <- 1

    params <- data.frame(
      name     = colnames(sub_mat),
      desc     = colnames(sub_mat),
      range    = maxs - mins + 2 * pad,
      minRange = mins - pad,
      maxRange = maxs + pad,
      stringsAsFactors = FALSE
    )
    rownames(params) <- sprintf("$P%d", seq_len(ncol(sub_mat)))

    flowCore::flowFrame(
      exprs      = sub_mat,
      parameters = methods::new("AnnotatedDataFrame", data = params),
      description = list(
        `$FIL`     = paste0(sample_name, ".fcs"),
        `GUID`     = sample_name,
        `$SRC`     = sample_name,
        `EXPORTED` = "export_flowjo_fcs.R"
      )
    )
  }

  safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

  sample_paths <- character(0)
  for (s in samp$levels) {
    idx <- cells$sample_id == s
    if (!any(idx)) next
    path <- file.path(outdir, paste0(safe_name(s), ".fcs"))
    suppressWarnings(flowCore::write.FCS(build_frame(mat[idx, , drop = FALSE], s), path))
    sample_paths <- c(sample_paths, path)
  }

  group_paths <- character(0)
  if (isTRUE(groups)) {
    if (!group_col %in% names(cells)) {
      log("  NOTE no '", group_col, "' column \u2014 per-group FCS files skipped")
    } else {
      gv <- as.character(cells[[group_col]])
      # An unmatched patient carries NA here. Those cells are still real and
      # still in _ALL_SAMPLES, so they get their own file rather than being
      # dropped silently from the per-group view.
      gv[is.na(gv) | !nzchar(trimws(gv))] <- "unassigned"
      lv <- sort(unique(gv))
      # One level means the group file would be a byte-for-byte copy of
      # _ALL_SAMPLES under a different name.
      if (length(lv) < 2L) {
        log("  NOTE only one '", group_col, "' level (", lv,
            ") \u2014 per-group FCS files skipped, _ALL_SAMPLES.fcs already is it")
      } else {
        for (g in lv) {
          idx <- gv == g
          if (!any(idx)) next
          path <- file.path(outdir, paste0("_GROUP_", safe_name(g), ".fcs"))
          suppressWarnings(flowCore::write.FCS(
            build_frame(mat[idx, , drop = FALSE], paste0("GROUP ", g)), path))
          group_paths <- c(group_paths, path)
          log("  group '", g, "': ", sum(idx), " cells -> ", basename(path))
        }
      }
    }
  }

  all_path <- character(0)
  if (isTRUE(concat)) {
    all_path <- file.path(outdir, "_ALL_SAMPLES.fcs")
    suppressWarnings(flowCore::write.FCS(build_frame(mat, "ALL_SAMPLES"), all_path))
  }

  # Order is deliberate -- FlowJo opens samples in list order, so the pooled
  # embedding lands first, the groups next (they are what the overlay is built
  # from), and the individual samples last.
  written <- c(all_path, group_paths, sample_paths)

  writeLines(normalizePath(written, winslash = "/", mustWork = FALSE),
             file.path(outdir, "samples.txt"))

  log("wrote ", length(written), " FlowJo FCS file(s) to ", outdir,
      " (", nrow(mat), " cells, ", length(markers) + 5L, " parameters): ",
      length(all_path), " pooled, ", length(group_paths), " group, ",
      length(sample_paths), " per-sample")
  log("  open _ALL_SAMPLES.fcs, plot UMAP-1 vs UMAP-2; ",
      "decode Population/SampleID/CohortID with population_codes.csv")
  if (length(group_paths))
    log("  overlay the _GROUP_*.fcs files in the Layout Editor to compare cohorts")

  invisible(written)
}

# --- command-line entry point ------------------------------------------------
# sys.nframe() is 0 only under `Rscript this-file.R`; source() runs the file
# inside a call frame, so this block stays inert when the pipeline imports it.
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(optparse))

  opts <- parse_args(OptionParser(option_list = list(
    make_option("--cells",  type = "character", default = "results/cells_umap.csv",
                help = "cells_umap.csv produced by run_cyraven() [%default]"),
    make_option("--outdir", type = "character", default = "output/flowjo",
                help = "directory to write FCS files into [%default]"),
    # dest is set explicitly: optparse's derived name for a flag containing a
    # dash is not worth guessing, and a silently-NULL option would make
    # --no-concat a no-op rather than an error.
    make_option("--no-concat", action = "store_true", default = FALSE,
                dest = "no_concat",
                help = "skip the concatenated _ALL_SAMPLES.fcs"),
    make_option("--no-groups", action = "store_true", default = FALSE,
                dest = "no_groups",
                help = "skip the per-group _GROUP_<cohort>.fcs files"),
    make_option("--group-col", type = "character", default = "cohort",
                dest = "group_col",
                help = "column defining the study group [%default]")
  )))

  export_flowjo_fcs(opts$cells, opts$outdir,
                    concat    = !isTRUE(opts$no_concat),
                    groups    = !isTRUE(opts$no_groups),
                    group_col = opts$group_col)
}
