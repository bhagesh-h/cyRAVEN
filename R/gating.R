# SECTION 4 -- AUTO-DERIVED GATE HIERARCHY
# =============================================================================

#' Derive the leukocyte scatter gate from the FSC-A density
#'
#' WHAT: lower FSC bound = the deepest valley in log10 FSC-A; upper bound = a
#' high percentile; SSC bounds = robust percentiles within the FSC window.
#'
#' WHY THIS IS COUNTER-INTUITIVE AND IMPORTANT: on this data a large majority of
#' recorded events (roughly two thirds in the test batch) are sub-cellular debris
#' sitting BELOW the FSC valley, and that debris is CD45-negative. A conventional
#' "draw a box low on FSC/SSC for lymphocytes" gate therefore selects debris, not
#' cells. Deriving the boundary from the valley puts the gate above the debris
#' mode wherever that mode happens to fall for a given instrument and sample.
#' @param ex Numeric expression matrix, events x channels.
#' @param sc The sc.
#' @param fsc_hi_q The fsc hi q. Default `0.999`.
#' @param ssc_q The ssc q. Default `c(0.005, 0.995)`.
#' @keywords internal
derive_scatter_gate <- function(ex, sc, fsc_hi_q = 0.999,
                               ssc_q = c(0.005, 0.995)) {
  if (!all(c("FSC-A", "SSC-A") %in% names(sc)))
    stop("FSC-A and SSC-A are required for the scatter gate but were not found")
  lf <- log10(pmax(ex[, sc[["FSC-A"]]], 1))
  ls <- log10(pmax(ex[, sc[["SSC-A"]]], 1))
  v <- density_valley(lf, bins = 240L, smooth = 5, min_gap_frac = 0.10)
  if (!is.finite(v)) {
    v <- as.numeric(quantile(lf, 0.60, na.rm = TRUE))
    warning("no FSC-A density valley found; using the 60th percentile (",
            round(v, 2), "), CHECK THE QC FIGURE")
  }
  fsc_hi <- as.numeric(quantile(lf, fsc_hi_q, na.rm = TRUE))
  inw <- lf > v & lf < fsc_hi
  sb <- as.numeric(quantile(ls[inw], ssc_q, na.rm = TRUE))
  list(fsc_lo = v, fsc_hi = fsc_hi, ssc_lo = sb[1], ssc_hi = sb[2],
       mask = inw & ls > sb[1] & ls < sb[2],
       fsc_log10 = lf, ssc_log10 = ls)
}

#' Derive the singlet band from the FSC-H/FSC-A ridge
#'
#' WHY: doublet exclusion assumes height and area track together for single
#' cells. The RATIO at which they track is instrument- and gain-dependent -- on
#' this cytometer the ridge sits near 0.56, not the 1.0 that a hardcoded gate
#' would assume, and a band centred on 1.0 discards essentially every real cell.
#' Centering on the measured median with a MAD-scaled width adapts to any
#' instrument. MAD is used rather than SD because doublets are a heavy tail.
#' @param ex Numeric expression matrix, events x channels.
#' @param sc The sc.
#' @param parent_mask The parent mask.
#' @param k Neighbourhood size, or number of clusters, depending on the function. Default `3`.
#' @keywords internal
derive_singlet_band <- function(ex, sc, parent_mask, k = 3) {
  # The two must be present AND be different columns. On a height-only
  # acquisition the reading stage points FSC-A at the height channel so the
  # scatter gate can run, which would otherwise make this ratio identically 1:
  # a MAD of zero, a band of zero width, and every event discarded.
  have_both <- all(c("FSC-H", "FSC-A") %in% names(sc)) &&
               !identical(unname(sc[["FSC-H"]]), unname(sc[["FSC-A"]]))
  if (!have_both) {
    log_msg("  no distinct FSC-H and FSC-A, singlet gate skipped ",
            "(all parent events retained)")
    return(list(lo = NA_real_, hi = NA_real_, ratio = NULL,
                mask = parent_mask, skipped = TRUE))
  }
  r <- ex[, sc[["FSC-H"]]] / pmax(ex[, sc[["FSC-A"]]], 1)
  rp <- r[parent_mask]
  med <- median(rp, na.rm = TRUE)
  md  <- mad(rp, na.rm = TRUE)
  if (!is.finite(md) || md == 0) md <- sd(rp, na.rm = TRUE)
  lo <- med - k * md; hi <- med + k * md
  list(lo = lo, hi = hi, ratio = r, median = med, mad = md,
       mask = parent_mask & r > lo & r < hi, skipped = FALSE)
}

#' Locate the viability dye among the available markers
#' WHY: the dye differs per panel (Zombie-NIR, Zombie-Violet, L/D, 7-AAD, DAPI,
#' PI, Live/Dead...). Detect by pattern so no panel needs code changes, but let
#' the user name it explicitly when the pattern is ambiguous.
#' @param markers Character vector of marker names to use.
#' @param explicit The explicit.
#' @keywords internal
detect_viability_marker <- function(markers, explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit)) {
    if (explicit %in% markers) return(explicit)
    warning("named viability marker '", explicit, "' not in panel; auto-detecting")
  }
  pat <- "zombie|live.?dead|^l/d$|^ld$|7.?aad|dapi|viability|propidium|^pi$|sytox|dracula"
  hit <- markers[grepl(pat, markers, ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

#' Apply the full gate hierarchy to one file
#'
#' Chain: all events -> leukocytes (scatter) -> single cells -> live -> CD45+.
#' Each step degrades gracefully: a missing viability dye or missing CD45 skips
#' that step with a warning rather than aborting, because a panel that lacks one
#' is still analysable.
#'
#' @param rd The rd.
#' @param cofactor Arcsinh cofactor. See [derive_cofactor()].
#' @param cfg The cfg. Default `list()`.
#' @param control_ref The control ref.
#' @param singlet_k The singlet k. Default `3`.
#' @param viability_name The viability name.
#' @param cd45_name The cd45 name. Default `"CD45"`.
#' @param transform Intensity transform from [make_transform()]. Defaults to
#'   arcsinh with `cofactor`, which is what every caller got before the
#'   transform became selectable.
#' @param overrides Optional per-marker override entries for THIS sample, as
#'   returned by indexing the config's `sample_overrides` block by sample id.
#'   See [sample_override()]. Absent by default, in which case every threshold is
#'   derived exactly as before.
#' @return list of masks, derived geometry, thresholds, and a tidy counts table
#' @export
apply_gate_hierarchy <- function(rd, cofactor, cfg = list(), control_ref = NULL,
                                 singlet_k = 3, viability_name = NULL,
                                 cd45_name = "CD45", transform = NULL,
                                 overrides = NULL) {
  ex <- rd$exprs; sc <- rd$scatter_cols; mc <- rd$marker_cols
  sg <- derive_scatter_gate(ex, sc)
  sb <- derive_singlet_band(ex, sc, sg$mask, k = singlet_k)

  # Default to arcsinh with the supplied cofactor so every existing caller keeps
  # its behaviour unchanged; a caller that wants logicle passes the object built
  # by make_transform() instead. The gate code below never learns which is in
  # force, which is the point of routing both through one closure.
  tr <- transform %||% make_transform("arcsinh", cofactor = cofactor)
  tf <- function(m) tr$fn(ex[, mc[[m]]], m)           # transform one marker

  # A --config thresholds: entry is {threshold, source, needs_review} (what
  # write_config() emits and what the population-marker loop below reads via
  # cfg_thr[[m]]$threshold) -- but resolve_threshold()'s cfg_value must be a
  # bare number. Unwrap here too, for the two markers (viability, CD45) gated
  # in THIS function rather than the spec-driven loop: passing the list
  # through unwrapped previously crashed resolve_threshold's is.finite()
  # check with "default method not implemented for type 'list'" the moment
  # anyone tried to override either from --config, rather than applying it.
  cfg_threshold <- function(v) if (is.list(v)) v$threshold else v

  # --- live gate -------------------------------------------------------------
  vmk <- detect_viability_marker(names(mc), viability_name)
  live <- sb$mask; v_thr <- NA_real_; v_src <- "skipped"
  if (!is.na(vmk)) {
    vx <- tf(vmk)
    rr <- resolve_threshold(vmk, vx[sb$mask], cfg_threshold(cfg[[vmk]]),
                            control_ref[[vmk]], fallback_q = 0.90,
                            override = sample_override(list(x = overrides), "x", vmk))
    v_thr <- rr$threshold; v_src <- rr$source
    live <- sb$mask & vx < v_thr        # dead cells are the BRIGHT tail
    log_msg("  viability '", vmk, "' threshold ", round(v_thr, 2), " (", v_src, ")")
  } else {
    log_msg("  no viability dye detected, live gate skipped")
  }

  # --- CD45 gate -------------------------------------------------------------
  cd45 <- live; c_thr <- NA_real_; c_src <- "skipped"; cd45_x <- NULL
  if (cd45_name %in% names(mc)) {
    cd45_x <- tf(cd45_name)
    rr <- resolve_threshold(cd45_name, cd45_x[live], cfg_threshold(cfg[[cd45_name]]),
                            control_ref[[cd45_name]],
                            override = sample_override(list(x = overrides), "x", cd45_name))
    c_thr <- rr$threshold; c_src <- rr$source
    cd45 <- live & cd45_x > c_thr
    log_msg("  CD45 threshold ", round(c_thr, 2), " (", c_src, ")")
  } else {
    warning("CD45 absent from this panel, all live cells treated as the parent")
  }

  masks <- list(all_events = rep(TRUE, nrow(ex)), leukocytes = sg$mask,
                single_cells = sb$mask, live_cells = live, cd45_pos = cd45)
  parent_of <- c(all_events = NA, leukocytes = "all_events",
                 single_cells = "leukocytes", live_cells = "single_cells",
                 cd45_pos = "live_cells")
  counts <- data.frame(
    sample_id = rd$sample_id, gate = names(masks),
    count = vapply(masks, sum, integer(1)),
    pct_of_parent = vapply(names(masks), function(g) {
      p <- parent_of[[g]]
      if (is.na(p)) 100 else 100 * sum(masks[[g]]) / max(1L, sum(masks[[p]]))
    }, numeric(1)),
    pct_of_all = 100 * vapply(masks, sum, integer(1)) / nrow(ex),
    row.names = NULL)

  list(masks = masks, counts = counts, scatter_gate = sg, singlet = sb,
       viability_marker = vmk, viability_threshold = v_thr, viability_source = v_src,
       cd45_threshold = c_thr, cd45_source = c_src, cd45_x = cd45_x,
       cofactor = cofactor)
}

#' Staining QC verdict per file
#'
#' WHY: an unstained or failed-staining file has essentially no CD45+ events.
#' Embedding it silently places a cloud of background-only cells in the shared
#' space and corrupts every frequency table. Flag and exclude, but still report.
#' @param gate The gate.
#' @param declared_control The declared control. Default `NA`.
#' @param min_cd45_pct The min cd45 pct. Default `5`.
#' @return list(..., qc_status, is_control, is_reference). The three differ and
#'   every difference is load-bearing:
#'     qc_status    -- "pass" | "control" | "failed". The REPORTING category.
#'                    "control" is an unstained reference tube; "failed" is a
#'                    biological sample whose staining or CD45 gate did not work.
#'                    Both are excluded, for different reasons, and a figure that
#'                    calls a failed patient sample a "control" is simply wrong.
#'     is_control   -- TRUE only for an actual control tube. NOT an exclusion flag:
#'                    exclusion is `!include`, which covers failures too.
#'     is_reference -- additionally usable as the UNSTAINED REFERENCE from which
#'                    thresholds for unimodal markers are taken.
#'   Only a sample the operator DECLARED as a control earns is_reference. A
#'   biological sample that merely failed staining is excluded, never promoted to
#'   reference: its distributions are a broken assay, and adopting them as the
#'   negative reference propagates one bad tube's failure into every threshold in
#'   the panel. When no sample map declares controls at all (declared = NA) the
#'   inference may still promote, because then it is the only signal available.
#' @param force_include when TRUE (--include-qc-failed), a DECLARED sample (not
#'   a control) that would otherwise fail is included anyway, using whatever
#'   threshold was already computed. It is still excluded when the sample map
#'   itself says this is a control tube -- that case is an unstained reference,
#'   and embedding it puts a cloud of background-only events in the shared
#'   space regardless of this flag. Only meant for the "too few samples to
#'   afford losing any" situation; the verdict text still records the forcing
#'   so it stays auditable in staining_qc.csv.
#' @export
staining_verdict <- function(gate, declared_control = NA, min_cd45_pct = 5,
                             force_include = FALSE) {
  pct <- gate$counts$pct_of_parent[gate$counts$gate == "cd45_pos"]
  # A sample map that says FALSE is a positive assertion: this is a biological
  # sample. Failures below are then failures, not discovered controls.
  may_reference <- !identical(declared_control, FALSE)
  if (isTRUE(declared_control))
    return(list(pct_cd45_of_live = pct, verdict = "declared control - reference only",
                qc_status = "control",
                include = FALSE, is_control = TRUE, is_reference = TRUE))
  force <- force_include && !may_reference

  # A CD45 threshold that came from the quantile fallback carries NO evidence of
  # a positive population: the fallback exists precisely because CD45 was
  # unimodal, and the fraction above a q-quantile cut is (1-q) BY CONSTRUCTION.
  # Testing that fraction against a floor is circular -- it would pass an
  # unstained file every time. Absence of a CD45 mode is itself the diagnosis.
  if (identical(gate$cd45_source, "quantile_fallback")) {
    if (force)
      return(list(pct_cd45_of_live = pct,
                  verdict = paste("no separable CD45+ mode (threshold from quantile",
                                  "fallback) - included anyway (--include-qc-failed);",
                                  "pct_of_cd45_pos for this sample carries NO evidence,",
                                  "treat it with caution"),
                  qc_status = "pass",
                  include = TRUE, is_control = FALSE, is_reference = FALSE))
    return(list(pct_cd45_of_live = pct,
                verdict = paste0("no separable CD45+ mode (threshold from quantile ",
                                 "fallback) - ",
                                 if (may_reference) "treated as unstained/control"
                                 else "FAILED (declared a sample, not a control)"),
                qc_status = if (may_reference) "control" else "failed",
                include = FALSE, is_control = may_reference,
                is_reference = may_reference))
  }
  if (identical(gate$cd45_source, "skipped"))
    return(list(pct_cd45_of_live = pct,
                verdict = "CD45 absent from panel - staining QC not applicable",
                qc_status = "pass",
                include = TRUE, is_control = FALSE, is_reference = FALSE))

  if (!is.finite(pct) || pct < min_cd45_pct) {
    if (force)
      return(list(pct_cd45_of_live = pct,
                  verdict = sprintf(paste("unstained/failed (%.2f%% CD45+ < %g%%) -",
                                         "included anyway (--include-qc-failed)"),
                                    pct, min_cd45_pct),
                  qc_status = "pass",
                  include = TRUE, is_control = FALSE, is_reference = FALSE))
    return(list(pct_cd45_of_live = pct,
                verdict = sprintf("unstained/failed (%.2f%% CD45+ < %g%%) - %s",
                                  pct, min_cd45_pct,
                                  if (may_reference) "excluded"
                                  else "FAILED (declared a sample, not a control)"),
                qc_status = if (may_reference) "control" else "failed",
                include = FALSE, is_control = may_reference,
                is_reference = may_reference))
  }
  list(pct_cd45_of_live = pct, verdict = "stained sample - included",
       qc_status = "pass",
       include = TRUE, is_control = FALSE, is_reference = FALSE)
}

# =============================================================================
