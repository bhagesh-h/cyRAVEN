# SECTION 10 -- LABELS FROM ANOTHER TOOL, TURNED INTO A GATE
# =============================================================================
#
# THE PROBLEM THIS SOLVES. A clustering tool finds a population. It exists in a
# data frame of cluster assignments and nowhere else. Nobody can sort on it,
# nobody can draw it at the instrument, and the next cohort has to be clustered
# again to find out whether it is there. explain_cluster() already converts a
# label into a sequence of two-marker polygon gates; what was missing was a way
# to hand it a label that did not come from this package, a way to find out
# whether the gate works on a donor it was not fitted to, and a file format the
# instrument can read.
#
# WHY LEAVE-ONE-DONOR-OUT AND NOT THE EXISTING HELD-OUT CELLS. explain_cluster()
# reserves a fraction of CELLS for evaluation, stratified by label. Those cells
# come from the same donors as the training cells, acquired in the same tubes on
# the same day, so they share every source of between-donor variation the gate
# will actually meet. A gate can score an excellent held-out F1 and still fail on
# the next patient. Refitting with one donor withheld and scoring on that donor
# measures the thing a reader wants to know, and the SPREAD across donors matters
# more than the mean: a gate with F1 0.9 on every donor and a gate averaging 0.9
# because it scores 1.0 on nine and 0.1 on one are not the same gate.
#
# PRIOR ART. Hypergate (Becht et al. 2019, Bioinformatics 35:301) established
# cluster-to-gate conversion, fitting a hyperrectangle per cluster. What is added
# here is the polygon hierarchy, which is the shape a sorter is actually driven
# in, the per-donor transferability distribution, and the gate file.

#' Read cell labels produced by another tool
#'
#' Accepts a CSV carrying one row per cell with a sample identifier, the event's
#' index within that sample's file, and a label. Column names are matched
#' loosely, so a table exported from a clustering package usually needs no
#' editing.
#'
#' WHY THE KEY MUST BE EXPLICIT AND NOT POSITIONAL. Tools subsample. cyCONDOR
#' takes `max_cell` events per file and cyRAVEN takes its own cap for the
#' embedding, so row *i* of one table is not row *i* of the other and a
#' positional join silently mislabels every cell. Requiring the event index makes
#' the mismatch an error instead of a result.
#'
#' @param path CSV path
#' @param sample_col,event_col,label_col column names; NULL matches by convention
#' @return data.frame(sample_id, event_index, label)
#' @export
read_external_labels <- function(path, sample_col = NULL, event_col = NULL,
                                 label_col = NULL) {
  if (!file.exists(path))
    stop("external label file not found: ", path, call. = FALSE)
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
  pick <- function(given, candidates, what) {
    if (!is.null(given)) {
      if (!given %in% names(d))
        stop("column '", given, "' not in ", basename(path), ". Columns present: ",
             paste(names(d), collapse = ", "), call. = FALSE)
      return(given)
    }
    hit <- intersect(tolower(candidates), tolower(names(d)))
    if (!length(hit))
      stop("no ", what, " column in ", basename(path),
           ". Expected one of: ", paste(candidates, collapse = ", "),
           ". Columns present: ", paste(names(d), collapse = ", "), call. = FALSE)
    names(d)[match(hit[1], tolower(names(d)))]
  }
  sc <- pick(sample_col, c("sample_id", "sample", "file", "filename"), "sample")
  ec <- pick(event_col, c("event_index", "event", "index", "cell_id", "event_id"),
             "event index")
  lc <- pick(label_col, c("label", "cluster", "metacluster", "cell_type",
                          "celltype", "population"), "label")
  out <- data.frame(sample_id = as.character(d[[sc]]),
                    event_index = suppressWarnings(as.integer(d[[ec]])),
                    label = as.character(d[[lc]]),
                    row.names = NULL, stringsAsFactors = FALSE)
  out <- out[!is.na(out$event_index) & !is.na(out$label) & nzchar(out$label), ,
             drop = FALSE]
  if (!nrow(out))
    stop("no usable rows in ", basename(path), call. = FALSE)
  log_msg("  external labels: ", nrow(out), " cell(s), ",
          length(unique(out$label)), " label(s), ",
          length(unique(out$sample_id)), " sample(s), from columns ",
          sc, "/", ec, "/", lc)
  out
}

#' Attach external labels to the embedding cell table
#'
#' Reports the overlap rather than assuming it. A join that matched a handful of
#' cells is the normal consequence of two tools having subsampled differently,
#' and it has to be visible before anything is fitted to the result.
#'
#' @param cells embedding cell table with sample_id and event_index
#' @param labels the table from [read_external_labels()]
#' @return `cells` with an `external_label` column, or NULL when nothing matched
#' @export
join_external_labels <- function(cells, labels) {
  if (is.null(cells) || !nrow(cells)) return(NULL)
  if (!all(c("sample_id", "event_index") %in% names(cells)))
    stop("the embedding cell table has no sample_id/event_index to join on",
         call. = FALSE)
  k1 <- paste(cells$sample_id, cells$event_index, sep = "\r")
  k2 <- paste(labels$sample_id, labels$event_index, sep = "\r")
  cells$external_label <- labels$label[match(k1, k2)]
  n <- sum(!is.na(cells$external_label))
  log_msg("  external labels matched ", n, " of ", nrow(cells),
          " embedded cell(s) (", round(100 * n / nrow(cells), 1), "%)")
  if (n < 100L) {
    log_msg("  NOTE too few cells matched to fit a gate. The join is on ",
            "sample_id and event_index; check that the label file indexes ",
            "events within each FILE, from the same numbering this run used, ",
            "and that the sample identifiers agree.")
    return(NULL)
  }
  if (n < 0.2 * nrow(cells))
    log_msg("  NOTE fewer than a fifth of embedded cells carry a label. Both ",
            "tools subsample independently, so this is expected; the gate is ",
            "fitted on the intersection.")
  cells
}

#' Apply a learned gating strategy to new cells
#'
#' Runs the polygons in order, each on the cells the previous one kept. This is
#' what makes a proposal executable: the same function scores the cells the gate
#' was fitted on, a held-out donor, and next year's cohort.
#'
#' @param X marker matrix carrying at least the columns each level names
#' @param strategy the list returned by [explain_cluster()]
#' @param depth stop after this many levels; NULL uses all of them
#' @return logical vector, TRUE for cells inside every gate
#' @export
apply_gate_strategy <- function(X, strategy, depth = NULL) {
  X <- as.matrix(X)
  alive <- rep(TRUE, nrow(X))
  lv <- strategy$levels
  if (!is.null(depth)) lv <- utils::head(lv, as.integer(depth))
  for (g in lv) {
    if (!all(g$pair %in% colnames(X)))
      stop("marker(s) missing from the matrix: ",
           paste(setdiff(g$pair, colnames(X)), collapse = ", "), call. = FALSE)
    inside <- point_in_polygon(X[alive, g$pair, drop = FALSE], g$polygon)
    alive[alive] <- inside
  }
  alive
}

#' Does the gate work on a donor it was not fitted to
#'
#' Refits the whole strategy with one donor withheld and scores it on that donor,
#' once per donor. Returns the per-donor scores and their spread.
#'
#' READ THE MINIMUM, NOT THE MEAN. A gate that transfers is one whose worst donor
#' is acceptable. The mean hides exactly the failure this function exists to
#' expose.
#'
#' COST, AND THE TWO CAPS THAT BOUND IT. Every fold refits the whole strategy, so
#' the work is one fit per donor per label and grows with both. On a cohort of
#' twenty-odd donors and half a dozen labels an uncapped version runs for hours,
#' which is not a validation statistic anyone will wait for.
#'
#' `max_donors` bounds the number of folds. The subset is drawn at random from the
#' fixed seed rather than taken as the largest donors, because a worst-donor
#' statistic computed only on the best-represented donors is optimistic in exactly
#' the direction that matters. `max_cells` subsamples each fold's TRAINING rows,
#' stratified by label; the held-out donor is always scored in full. Both caps are
#' logged when they bind, so a summary over eight donors is never mistaken for one
#' over twenty.
#'
#' @param X marker matrix
#' @param y 0/1 target indicator
#' @param donor grouping vector, one entry per row of `X`
#' @param min_donors refuse below this many donors, where the statistic would be
#'   an anecdote
#' @param max_donors ceiling on folds
#' @param max_cells ceiling on training rows per fold
#' @param seed seed for the local RNG stream
#' @param ... passed to [explain_cluster()]
#' @return list(per_donor = data.frame, summary = data.frame), or NULL
#' @export
gate_transferability <- function(X, y, donor, min_donors = 3L, max_donors = 8L,
                                 max_cells = 20000L, seed = 42L, ...) {
  X <- as.matrix(X); y <- as.integer(as.logical(y))
  donor <- as.character(donor)
  ds <- sort(unique(donor[!is.na(donor)]))
  if (length(ds) < min_donors) {
    log_msg("  transferability: ", length(ds), " donor(s), need at least ",
            min_donors, " -- skipped")
    return(NULL)
  }

  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  n_all <- length(ds)
  if (n_all > max_donors) {
    ds <- sort(sample(ds, max_donors))
    log_msg("  transferability: ", n_all, " donors, scoring ", max_donors,
            " drawn at random (--transfer-max-donors); the summary describes ",
            "those ", max_donors, ", not the cohort")
  }
  .t0 <- Sys.time()

  rows <- list()
  for (d in ds) {
    tr <- donor != d; te <- !tr
    if (sum(tr) > max_cells) {
      # Stratified, so a rare target is not thinned out of existence by a cap
      # aimed at the background.
      itr <- which(tr)
      pos <- itr[y[itr] == 1]; neg <- itr[y[itr] == 0]
      k_pos <- min(length(pos), max(20L, round(max_cells * length(pos) / length(itr))))
      k_neg <- min(length(neg), max_cells - k_pos)
      keep <- c(sample(pos, k_pos), sample(neg, k_neg))
      tr <- rep(FALSE, length(y)); tr[keep] <- TRUE
    }
    if (sum(y[tr] == 1) < 20L || sum(y[te] == 1) < 10L) {
      rows[[length(rows) + 1L]] <- data.frame(
        donor = d, n_cells = sum(te), n_targets = sum(y[te] == 1),
        precision = NA_real_, recall = NA_real_, f1 = NA_real_,
        note = "too few target cells to fit or to score",
        stringsAsFactors = FALSE)
      next
    }
    st <- tryCatch(explain_cluster(X[tr, , drop = FALSE], y[tr], ...),
                   error = function(e) NULL)
    if (is.null(st)) {
      rows[[length(rows) + 1L]] <- data.frame(
        donor = d, n_cells = sum(te), n_targets = sum(y[te] == 1),
        precision = NA_real_, recall = NA_real_, f1 = NA_real_,
        note = "no gate could be fitted without this donor",
        stringsAsFactors = FALSE)
      next
    }
    inside <- tryCatch(apply_gate_strategy(X[te, , drop = FALSE], st),
                       error = function(e) NULL)
    if (is.null(inside)) next
    m <- gate_metrics(inside, y[te], sum(y[te] == 1))
    rows[[length(rows) + 1L]] <- data.frame(
      donor = d, n_cells = sum(te), n_targets = sum(y[te] == 1),
      precision = round(unname(m[["precision"]]), 4),
      recall = round(unname(m[["recall"]]), 4),
      f1 = round(unname(m[["f1"]]), 4),
      note = paste(vapply(st$levels, function(l) paste(l$pair, collapse = "/"),
                          character(1)), collapse = " then "),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  per <- do.call(rbind, rows)
  f <- per$f1[is.finite(per$f1)]
  summ <- data.frame(
    n_donors_in_cohort = n_all, n_donors = nrow(per), n_scored = length(f),
    f1_min = if (length(f)) round(min(f), 4) else NA_real_,
    f1_median = if (length(f)) round(median(f), 4) else NA_real_,
    f1_max = if (length(f)) round(max(f), 4) else NA_real_,
    f1_iqr = if (length(f)) round(diff(quantile(f, c(0.25, 0.75))), 4) else NA_real_,
    seconds = round(as.numeric(difftime(Sys.time(), .t0, units = "secs"))),
    stringsAsFactors = FALSE)
  list(per_donor = per, summary = summ)
}

#' Learn and validate a gate for each externally supplied label
#'
#' @param cells embedding cell table carrying `external_label`
#' @param features marker columns to gate on
#' @param min_cells smallest label worth gating
#' @param max_labels ceiling on how many labels to process
#' @param donor_col column holding the donor identity for the transferability
#'   split
#' @param max_donors,transfer_max_cells caps passed to [gate_transferability()]
#' @param seed seed for the local RNG stream
#' @param ... passed to [explain_cluster()]
#' @return list(summary, polygons, transfer, transfer_summary, strategies), or NULL
#' @export
explain_external_labels <- function(cells, features, min_cells = 200L,
                                    max_labels = 6L, donor_col = "sample_id",
                                    max_donors = 8L, transfer_max_cells = 20000L,
                                    seed = 42L, ...) {
  feats <- intersect(features, names(cells))
  if (length(feats) < 2L) return(NULL)
  lab <- cells$external_label
  tab <- sort(table(lab[!is.na(lab)]), decreasing = TRUE)
  tab <- tab[tab >= min_cells]
  if (!length(tab)) {
    log_msg("  external labels: none reaches ", min_cells, " cells -- nothing to gate")
    return(NULL)
  }
  n_drop <- max(0L, length(tab) - max_labels)
  if (n_drop > 0L)
    log_msg("  external labels: ", length(tab), " qualify; gating the ",
            max_labels, " largest and skipping ", n_drop)
  tab <- utils::head(tab, max_labels)

  X <- as.matrix(cells[, feats, drop = FALSE])
  donor <- if (donor_col %in% names(cells)) as.character(cells[[donor_col]]) else NULL

  srows <- list(); prows <- list(); trows <- list(); tsum <- list()
  strategies <- list()
  for (L in names(tab)) {
    y <- as.integer(!is.na(lab) & lab == L)
    st <- tryCatch(explain_cluster(X, y, ...), error = function(e) NULL)
    if (is.null(st)) {
      log_msg("  external labels: '", L, "' -- no gate could be fitted")
      next
    }
    strategies[[L]] <- st
    s <- st$summary; s$label <- L; s$n_cells <- as.integer(tab[[L]])
    srows[[length(srows) + 1L]] <- s
    for (lv in st$levels) {
      p <- lv$polygon
      if (is.null(p)) next
      prows[[length(prows) + 1L]] <- data.frame(
        label = L, depth = lv$depth, marker_x = lv$pair[1L],
        marker_y = lv$pair[2L], vertex = seq_len(nrow(p)),
        x_transformed = round(p[, 1L], 6), y_transformed = round(p[, 2L], 6),
        stringsAsFactors = FALSE)
    }
    b <- s[s$depth == st$best_depth, , drop = FALSE]
    log_msg("  external labels: '", L, "' -> ",
            paste(vapply(st$levels, function(l) paste(l$pair, collapse = "/"),
                         character(1)), collapse = " then "),
            "  (held-out F1 ", format(round(b$cumulative_f1, 3), nsmall = 3), ")")

    if (!is.null(donor)) {
      tf <- tryCatch(gate_transferability(X, y, donor, max_donors = max_donors,
                                          max_cells = transfer_max_cells,
                                          seed = seed, ...),
                     error = function(e) NULL)
      if (!is.null(tf)) {
        pd <- tf$per_donor; pd$label <- L
        trows[[length(trows) + 1L]] <- pd
        ts <- tf$summary; ts$label <- L
        tsum[[length(tsum) + 1L]] <- ts
        log_msg("    across donors: F1 median ", ts$f1_median, ", worst ",
                ts$f1_min, " over ", ts$n_scored, " donor(s), ", ts$seconds, "s")
      }
    }
  }
  if (!length(srows)) return(NULL)
  list(summary = do.call(rbind, srows),
       polygons = if (length(prows)) do.call(rbind, prows) else NULL,
       transfer = if (length(trows)) do.call(rbind, trows) else NULL,
       transfer_summary = if (length(tsum)) do.call(rbind, tsum) else NULL,
       strategies = strategies)
}
