# Explore mode, part 2: the orchestrator and its figures.
#
# The helpers are in R/explore.R, which also carries the design note for the
# whole feature. This file assembles them into one stage and writes everything
# under <outdir>/explore/.

#' Run explore mode
#'
#' Unsupervised discovery over every eligible channel, writing into
#' `<outdir>/explore/`. Every argument beyond `reads`, `fpr`, `opt` and `outdir`
#' is optional: supplied, it makes the output better; absent, the run degrades to
#' what a standalone clusterer would produce and records that it did so in
#' `explore_provenance.csv`.
#'
#' @param reads Named list of reads.
#' @param fpr Panel resolution.
#' @param opt Option list.
#' @param outdir Run output directory.
#' @param transforms Per-panel transforms. Derived here if absent.
#' @param gates Per-sample gate objects, for thresholds and the viability marker.
#' @param verdicts Per-sample staining QC verdicts. Recorded, NOT used to
#'   exclude: explore does not depend on the CD45 gate, so a sample the declared
#'   path had to drop still contributes here. That is a capability rather than an
#'   oversight, which is why it becomes a column instead of a filter.
#' @param pops Per-sample scored populations, for the bridge tables.
#' @param group_of Named character vector, sample id -> group.
#' @param confounding One-row data frame from the batch/group confounding check.
#' @param file_paths Full paths of the FCS files the run was given. Needed only
#'   when the declared pipeline has already freed the event matrices, because
#'   `reads[[s]]` records a basename and cannot be re-opened on its own.
#' @return Invisibly, `list(dir, files, gaps)`. `gaps` is what explore has to say
#'   about the declared specification -- populations spanning several clusters,
#'   and clusters no population covers. The caller writes it to `spec_gaps.csv`
#'   only under `--maybe-learn`; this function never writes outside `explore/`.
#' @keywords internal
run_explore <- function(reads, fpr, opt, outdir, transforms = NULL, gates = NULL,
                        verdicts = NULL, pops = NULL, group_of = NULL,
                        confounding = NULL, file_paths = character(0),
                        total_counts = NULL) {
  ex_dir <- file.path(outdir, "explore")
  dir.create(ex_dir, showWarnings = FALSE, recursive = TRUE)
  seed <- opt$seed %||% 42L
  prov <- list()
  note <- function(item, value, why) {
    prov[[length(prov) + 1L]] <<- data.frame(
      item = item, value = as.character(value), basis = why,
      stringsAsFactors = FALSE)
  }
  gaps <- list()
  written <- character(0)
  wcsv <- function(df, name) {
    utils::write.csv(df, file.path(ex_dir, name), row.names = FALSE)
    written <<- c(written, name)
  }

  for (p in fpr$panels) {
    sids <- p$samples[vapply(p$samples, function(s) !is.null(reads[[s]]), TRUE)]
    if (!length(sids)) next
    tag <- if (length(fpr$panels) > 1L) paste0("_", p$name) else ""
    log_msg("explore: panel ", p$name, ", ", length(sids), " sample(s)")

    tr <- transforms[[p$name]]
    if (is.null(tr)) {
      cf <- derive_cofactor_pooled(reads, sids)
      tr <- make_transform("arcsinh", cofactor = cf)
      note("transform", sprintf("arcsinh cofactor %.4g", cf),
           "derived here; the declared pipeline did not run")
    } else {
      note("transform", tr$method,
           "reused from the declared pipeline, cofactor estimated on this panel")
    }

    scat <- names(reads[[sids[1]]]$scatter_cols %||% list())
    prefer <- .split_opt(opt[["explore_markers", exact = TRUE]])
    excl <- .split_opt(opt[["explore_exclude", exact = TRUE]])
    feats <- explore_features(p$markers, scat, prefer = prefer, exclude = excl)
    if (length(feats) < 2L) {
      log_msg("explore: fewer than 2 eligible channels, skipping ", p$name)
      next
    }
    note("features", paste(feats, collapse = "|"),
         if (is.null(prefer))
           "every eligible channel, including scatter and viability"
         else "--explore-markers")

    # The declared pipeline frees reads[[s]]$exprs after scoring, unless
    # --keep-exprs: holding every event matrix for a whole cohort is what makes
    # a run run out of memory. Explore therefore re-reads whatever it needs, one
    # file at a time, subsampling immediately so only the draw is ever held.
    # That costs I/O and saves memory, and it is also what lets explore run at
    # all in a pipeline whose primary path is done with the raw events.
    # rd$file is a BASENAME, not a path, so it cannot be re-opened on its own.
    # `file_paths` carries the full paths the run was given; match on basename.
    .path_of <- function(s) {
      bn <- reads[[s]]$file %||% NA_character_
      if (is.na(bn)) return(NA_character_)
      if (file.exists(bn)) return(bn)
      hit <- file_paths[basename(file_paths) == basename(bn)]
      if (length(hit)) hit[1] else NA_character_
    }
    .events <- function(s) {
      rd <- reads[[s]]
      if (!is.null(rd$exprs)) return(rd)
      f <- .path_of(s)
      if (is.na(f)) {
        log_msg("  explore: cannot locate the file for ", s,
                " to re-read (", reads[[s]]$file %||% "no filename recorded",
                "); pass --keep-exprs to keep the events in memory instead")
        return(NULL)
      }
      re <- try(read_fcs_resolved(f, sample_id = s,
                                  max_events = opt$max_events_per_file %||% 0L),
                silent = TRUE)
      if (inherits(re, "try-error")) {
        log_msg("  explore: re-reading ", basename(f), " failed: ",
                conditionMessage(attr(re, "condition")))
        return(NULL)
      }
      re$exprs <- maybe_compensate(re$exprs, re$keywords)
      re
    }
    n_of <- function(s) {
      rd <- reads[[s]]
      if (!is.null(rd$exprs)) return(nrow(rd$exprs))
      as.integer(rd$n_events %||% 0L)
    }
    navail <- vapply(sids, n_of, integer(1))
    if (all(navail <= 0L)) {
      log_msg("explore: no event counts available for ", p$name, ", skipping")
      next
    }
    planned <- plan_subsample(setNames(navail, sids),
                              cap = opt$explore_cells_per_sample %||% 6000L,
                              total_cap = opt$explore_max_cells %||% 150000L)
    n_reread <- 0L
    rows <- list()
    for (s in sids) {
      if (is.na(planned[[s]]) || planned[[s]] < 1L) next
      rd <- .events(s)
      if (is.null(rd)) {
        log_msg("  explore: no events available for ", s, ", skipped")
        next
      }
      if (is.null(reads[[s]]$exprs)) n_reread <- n_reread + 1L
      n <- nrow(rd$exprs)
      take <- withr::with_seed(seed, sort(sample.int(n, min(n, planned[[s]]))))
      # `take` is passed in rather than applied to the result: see the note in
      # explore_matrix(). The draw is identical either way; only the peak memory
      # differs, and on a full-depth run the difference is what decides whether
      # the container survives.
      X <- explore_matrix(rd, feats, tr, rows = take)
      if (is.null(X) || !ncol(X)) next
      # check.names = FALSE throughout: "HLA-DR" and "FSC-A" are legal marker
      # names and R would silently rewrite them to "HLA.DR" and "FSC.A", after
      # which every lookup against `feats` misses and the features vanish.
      d <- as.data.frame(X, check.names = FALSE)
      d$sample_id <- s
      d$event_index <- take
      rows[[s]] <- d
      rm(rd, X)
    }
    if (!length(rows)) {
      log_msg("explore: no sample yielded events for ", p$name, ", skipping")
      next
    }
    if (n_reread)
      note("events_reread", n_reread,
           "the declared run had already freed the event matrices; re-read one file at a time")
    cells <- as.data.frame(data.table::rbindlist(rows, use.names = TRUE, fill = TRUE),
                           check.names = FALSE)
    rm(rows)
    fcols <- intersect(feats, colnames(cells))
    note("cells_assembled", nrow(cells),
         sprintf("%d sample(s), UNGATED, equalised at %d events per sample",
                 length(sids), min(unlist(planned))))

    # Per-sample thresholds come from thresholds_used.csv, not from `gates`: the
    # gate object carries masks and geometry, while the resolved per-marker cuts
    # are produced later in scoring and written to that table. Reading the table
    # also keeps this decoupled from the gate object's internals.
    thr_by_sample <- list()
    if (!is.null(gates)) {
      tp <- file.path(outdir, "thresholds_used.csv")
      if (file.exists(tp)) {
        tt <- try(utils::read.csv(tp, stringsAsFactors = FALSE,
                                  check.names = FALSE), silent = TRUE)
        if (!inherits(tt, "try-error") &&
            all(c("sample_id", "marker", "threshold") %in% names(tt))) {
          tt <- tt[is.finite(tt$threshold), , drop = FALSE]
          if ("panel" %in% names(tt)) tt <- tt[tt$panel == p$name, , drop = FALSE]
          for (s in intersect(sids, unique(tt$sample_id))) {
            sub <- tt[tt$sample_id == s, , drop = FALSE]
            v <- setNames(sub$threshold, sub$marker)
            if (length(v)) thr_by_sample[[s]] <- v
          }
        }
      }
    }
    have_thr <- length(thr_by_sample) > 0L
    note("positivity_basis",
         if (have_thr) "per-sample thresholds" else "pooled per-feature medians",
         if (have_thr)
           sprintf("thresholds for %d of %d sample(s); a call is against that sample's own cut",
                   length(thr_by_sample), length(sids))
         else "no declared pipeline ran; this is what a standalone clusterer has")

    leuko <- if ("CD45" %in% fcols) "CD45" else NULL
    viab <- NULL
    if (!is.null(gates)) {
      vm <- unique(stats::na.omit(vapply(sids, function(s)
        gates[[s]]$viability_marker %||% NA_character_, character(1))))
      vm <- intersect(vm, fcols)
      if (length(vm)) viab <- vm[1]
    }

    qc_tab <- NULL
    keep_idx <- seq_len(nrow(cells))
    if (!isTRUE(opt$no_explore_qc)) {
      cc <- run_unsupervised_clusters(cells, fcols,
                                      n_clusters = opt$explore_qc_k %||% 20L,
                                      grid = opt$explore_grid %||% 10L,
                                      seed = seed, max_cells = nrow(cells))
      if (is.null(cc)) {
        note("qc_gate", "skipped", "coarse clustering returned nothing")
      } else {
        cl0 <- cc$cluster
        Xm <- as.matrix(cells[, fcols, drop = FALSE])
        uk <- sort(unique(cl0))
        pos0 <- explore_positivity(Xm, cl0, cells$sample_id, thr_by_sample)
        med0 <- t(vapply(uk, function(k)
          apply(Xm[cl0 == k, , drop = FALSE], 2, stats::median, na.rm = TRUE),
          numeric(length(fcols))))
        dimnames(med0) <- list(paste0("k", uk), fcols)
        sizes0 <- as.integer(table(cl0)[as.character(uk)])
        qc_tab <- explore_qc_gate(pos0, med0, sizes0, leukocyte = leuko,
                                  viability = viab, have_thresholds = have_thr,
                                  drop_dead = !isTRUE(opt$explore_keep_dead),
                                  seed = seed)
        qc_tab <- cbind(qc_tab, round(med0[qc_tab$cluster, , drop = FALSE], 4))
        wcsv(qc_tab, sprintf("explore_qc_clusters%s.csv", tag))

        kept <- as.integer(sub("^k", "", qc_tab$cluster[qc_tab$call == "keep"]))
        keep_idx <- which(cl0 %in% kept)
        frac <- length(keep_idx) / nrow(cells)
        log_msg(sprintf("explore: QC gate keeps %d of %d cells (%.1f%%)",
                        length(keep_idx), nrow(cells), 100 * frac))
        for (cn in setdiff(unique(qc_tab$call), "keep"))
          log_msg("  dropped as ", cn, ": ",
                  sum(qc_tab$cells[qc_tab$call == cn]), " cells")
        note("qc_gate", sprintf("%d of %d kept (%.1f%%)", length(keep_idx),
                                nrow(cells), 100 * frac),
             if (have_thr)
               "whole clusters judged against each sample's own thresholds"
             else "whole clusters judged by 2-means on per-cluster medians")

        min_frac <- opt$explore_min_retained %||% 0.05
        if (frac < min_frac) {
          log_msg("explore: the gate would keep under ", round(100 * min_frac),
                  "%. Markers are likelier mis-named than the data bad; ",
                  "keeping every cell and flagging it.")
          keep_idx <- seq_len(nrow(cells))
          note("qc_gate_override", "gate ignored",
               sprintf("retention below --explore-min-retained (%.2f)", min_frac))
        }
      }
    } else {
      note("qc_gate", "disabled", "--no-explore-qc")
    }

    cells <- cells[keep_idx, , drop = FALSE]
    if (nrow(cells) < 50L) {
      log_msg("explore: too few cells after the gate, skipping ", p$name)
      next
    }

    # Re-equalise so no donor dominates the embedding. cyCONDOR does not, and
    # its faceted panels then differ in density for a reason that is group size
    # rather than biology.
    if (!isTRUE(opt$no_explore_equalise)) {
      per <- table(cells$sample_id)
      m <- min(per)
      idx <- unlist(lapply(names(per), function(s) {
        w <- which(cells$sample_id == s)
        withr::with_seed(seed, sample(w, m))
      }))
      note("equalised", sprintf("%d cells per sample", m),
           sprintf("%d dropped so no donor dominates the embedding",
                   nrow(cells) - length(idx)))
      cells <- cells[sort(idx), , drop = FALSE]
    }

    M <- as.matrix(cells[, fcols, drop = FALSE])
    emb <- run_umap(M, seed = seed,
                    n_neighbors = opt$explore_neighbors %||% 30L,
                    min_dist = opt$explore_min_dist %||% 0.1,
                    n_threads = max(1L, opt$threads %||% 1L))
    # run_umap() returns $coords, not $layout. Assigning NULL to a data.frame
    # column is a silent no-op, so getting this wrong loses the columns without
    # an error until something later selects them.
    cells$umap_1 <- emb$coords[, 1]
    cells$umap_2 <- emb$coords[, 2]

    cl <- run_unsupervised_clusters(cells, fcols,
                                    n_clusters = opt$explore_k %||% 20L,
                                    grid = opt$explore_grid %||% 10L,
                                    seed = seed, max_cells = nrow(cells))
    if (is.null(cl)) {
      log_msg("explore: clustering returned nothing for ", p$name)
      next
    }
    cells$cluster <- paste0("k", cl$cluster)
    note("clusters", length(unique(cells$cluster)),
         sprintf("SOM %dx%d + consensus metaclustering, k = %d",
                 opt$explore_grid %||% 10L, opt$explore_grid %||% 10L,
                 opt$explore_k %||% 20L))

    if (!is.null(group_of)) cells$group <- unname(group_of[cells$sample_id])
    if (!is.null(verdicts)) cells$staining_qc_verdict <-
      vapply(cells$sample_id, function(s)
        if (isTRUE(verdicts[[s]]$include)) "pass" else "failed", character(1))

    keepcols <- c("sample_id", "event_index", "cluster", "umap_1", "umap_2",
                  intersect(c("group", "staining_qc_verdict"), colnames(cells)))
    wcsv(cells[, keepcols, drop = FALSE], sprintf("explore_cells%s.csv", tag))

    Xk <- as.matrix(cells[, fcols, drop = FALSE])
    kv <- cells$cluster
    ks <- sort(unique(kv))
    pos <- explore_positivity(Xk, kv, cells$sample_id, thr_by_sample)
    med <- t(vapply(ks, function(k)
      apply(Xk[kv == k, , drop = FALSE], 2, stats::median, na.rm = TRUE),
      numeric(length(fcols))))
    dimnames(med) <- list(ks, fcols)
    pheno <- explore_phenotype(pos[ks, , drop = FALSE])
    nk <- as.integer(table(kv)[ks])

    prof <- data.frame(cluster = ks, cells = nk,
                       pct_of_gated = round(100 * nk / length(kv), 3),
                       phenotype = unname(pheno[ks]),
                       phenotype_basis = if (have_thr) "per-sample thresholds"
                                         else "pooled medians",
                       stringsAsFactors = FALSE)
    fp <- as.data.frame(round(pos[ks, , drop = FALSE], 4), check.names = FALSE)
    names(fp) <- paste0("frac_pos.", fcols)
    mp <- as.data.frame(round(med[ks, , drop = FALSE], 4), check.names = FALSE)
    names(mp) <- paste0("median.", fcols)
    prof <- cbind(prof, fp, mp)
    wcsv(prof, sprintf("explore_cluster_profile%s.csv", tag))

    # ---- the subset each cluster's marker profile matches --------------------
    # Written here, from the profile alone, so it exists whether or not a
    # declared specification was scored: an explore-only run has no populations
    # to cross-tabulate against, and "which subset is this cluster" is the
    # question it most needs answered. The declared analysis is not touched --
    # this labels the CLUSTER, never the cell. See R/subsets.R.
    subs <- tryCatch(annotate_clusters_with_subsets(prof),
                     error = function(e) {
                       log_msg("  WARNING cluster subset annotation failed: ",
                               conditionMessage(e), ", every other output is ",
                               "unaffected")
                       NULL
                     })
    if (!is.null(subs) && nrow(subs)) {
      wcsv(subs, sprintf("explore_cluster_subsets%s.csv", tag))
      .named <- sum(!is.na(subs$subset_label))
      log_msg("explore: ", .named, " of ", nrow(subs), " cluster(s) match an ",
              "immune subset on their marker profile; the rest match none of ",
              "the definitions this panel can express")
      for (.u in unreachable_subsets(sub("^frac_pos[.]", "",
                                         grep("^frac_pos[.]", names(prof),
                                              value = TRUE))))
        log_msg("  NOT definable on this panel -- ", .u)
    }

    tab <- table(cells$sample_id, kv)
    ab <- do.call(rbind, lapply(rownames(tab), function(s) {
      n <- sum(tab[s, ])
      data.frame(sample_id = s, population = colnames(tab),
                 count = as.integer(tab[s, ]), n_parent_events = n,
                 pct_of_gated = 100 * as.integer(tab[s, ]) / n,
                 stringsAsFactors = FALSE)
    }))
    cu <- counting_uncertainty(ab$count, ab$n_parent_events)
    ab$u_counting_pct_points <- round(cu$u, 4)
    ab$detection <- cu$detection
    ab$qc_status <- "pass"
    ab$phenotype <- prof$phenotype[match(ab$population, prof$cluster)]
    if (!is.null(verdicts)) ab$staining_qc_verdict <-
      vapply(ab$sample_id, function(s)
        if (isTRUE(verdicts[[s]]$include)) "pass" else "failed", character(1))

    # --total-counts: re-express each cluster as a cell number. pct_of_gated is
    # a share and shares sum to 100, so a frequency table cannot separate "this
    # cluster expanded" from "everything else contracted". Multiplying by the
    # acquisition's own yield lifts that constraint. Added as extra columns
    # only -- pct_of_gated above is untouched, because the product carries the
    # external instrument's error as well as this pipeline's.
    if (!is.null(total_counts) && nrow(total_counts)) {
      ab <- attach_absolute_cells(ab, total_counts, share_col = "pct_of_gated")
      n_cov <- length(intersect(unique(ab$sample_id), total_counts$sample_id))
      log_msg("explore: --total-counts covers ", n_cov, " of ",
              length(unique(ab$sample_id)), " embedded sample(s); ",
              "cells_absolute written alongside pct_of_gated")
      if (n_cov < length(unique(ab$sample_id)))
        note("absolute_counts_coverage", "partial",
             paste0(length(unique(ab$sample_id)) - n_cov, " embedded sample(s) ",
                    "have no external total; their cells_absolute is NA and they ",
                    "drop out of the absolute test only"))
    }
    wcsv(ab, sprintf("explore_cluster_abundance%s.csv", tag))

    # --no-group-tests is honoured here as well as in the declared path. It is
    # documented as skipping every between-group test, and a user who set it
    # because the design cannot carry one does not mean "except on the
    # clusters" -- explore's tests are the same test on a different partition
    # of the same cells, and are just as unsupportable when the donors run out.
    if (isTRUE(opt$no_group_tests)) {
      note("group_statistics", "not run", "--no-group-tests was set")
    } else if (!is.null(group_of) &&
        length(unique(stats::na.omit(unname(group_of[sids])))) > 1L) {
      st <- try(stats_group_comparison(ab, group_of,
                                       reference = opt$reference_group,
                                       value_col = "pct_of_gated"), silent = TRUE)
      if (!inherits(st, "try-error") && !is.null(st) && nrow(st)) {
        st$phenotype <- prof$phenotype[match(st$population, prof$cluster)]
        # The confounding verdict travels WITH the statistics. Without it a
        # reader sees q < 0.05 and believes it, which is the failure the
        # declared pipeline refuses to permit.
        if (!is.null(confounding) && nrow(confounding)) {
          st$batch_group_cramers_v <- confounding$cramers_v[1]
          st$batch_group_verdict <- confounding$verdict[1]
        }
        wcsv(st, sprintf("explore_cluster_stats%s.csv", tag))
        nsig <- sum(st$significant_BH %in% c(TRUE, "TRUE"), na.rm = TRUE)
        log_msg("explore: ", nrow(st), " cluster test(s), ", nsig,
                " surviving BH")
        if (!is.null(confounding) && nrow(confounding) &&
            grepl("SEVERE", confounding$verdict[1] %||% "")) {
          log_msg("  NOTE: ", confounding$verdict[1])
          log_msg("  Every one of these is inseparable from acquisition date.")
        }
      }

      # THE SAME COMPARISON ON CELL NUMBERS RATHER THAN SHARES.
      # Written as its own file and never merged into the one above. A cluster
      # that moves on cells_absolute but not on pct_of_gated is the whole
      # reason the external total was supplied: it is the signature of every
      # population changing together, which a share cannot represent. The
      # converse -- significant on the share, not on the number -- says the
      # composition shifted while the compartment did not, so the result is
      # about redistribution. Merging the two tables would erase that
      # distinction, which is the only thing this route adds.
      if ("cells_absolute" %in% names(ab) && any(is.finite(ab$cells_absolute))) {
        sa <- try(stats_group_comparison(ab, group_of,
                                         reference = opt$reference_group,
                                         value_col = "cells_absolute"), silent = TRUE)
        if (!inherits(sa, "try-error") && !is.null(sa) && nrow(sa)) {
          sa$phenotype <- prof$phenotype[match(sa$population, prof$cluster)]
          sa$count_basis <- "dual-platform: share x external total"
          if (!is.null(confounding) && nrow(confounding)) {
            sa$batch_group_cramers_v <- confounding$cramers_v[1]
            sa$batch_group_verdict <- confounding$verdict[1]
          }
          wcsv(sa, sprintf("explore_cluster_stats_absolute%s.csv", tag))
          nsa <- sum(sa$significant_BH %in% c(TRUE, "TRUE"), na.rm = TRUE)
          log_msg("explore: ", nrow(sa), " absolute-count test(s), ", nsa,
                  " surviving BH")

          # The concordance is the readable summary of the paragraph above, so
          # it is computed rather than left to the reader to join by hand.
          if (exists("st", inherits = FALSE) && !inherits(st, "try-error") &&
              !is.null(st) && nrow(st)) {
            # EVERY cluster, not only the significant ones. On a cohort this
            # size nothing survives BH, and a table that is written only when
            # something does would be absent exactly when the comparison is
            # most worth seeing. Both adjusted p-values sit side by side so the
            # divergence is readable whether or not either crosses a threshold.
            best_q <- function(d) {
              q <- suppressWarnings(as.numeric(d$p_adj_BH))
              tapply(q, as.character(d$population),
                     function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
            }
            qa <- best_q(sa); qr <- best_q(st)
            cl <- union(names(qa), names(qr))
            q_abs <- unname(qa[cl]); q_shr <- unname(qr[cl])
            sig_a <- !is.na(q_abs) & q_abs < 0.05
            sig_r <- !is.na(q_shr) & q_shr < 0.05
            cc <- data.frame(
              cluster    = cl,
              q_share    = round(q_shr, 5),
              q_absolute = round(q_abs, 5),
              verdict    = ifelse(sig_a & sig_r, "both",
                           ifelse(sig_a, "absolute only",
                           ifelse(sig_r, "share only", "neither"))),
              stringsAsFactors = FALSE)
            cc$reading <- c(
              "absolute only" =
                "compartment changed size; composition did not. Only an external total can show this",
              "share only" =
                "composition shifted without the compartment changing size; a redistribution",
              "both" = "moved on both measures",
              "neither" = "no difference resolved on either measure")[cc$verdict]
            cc$phenotype <- prof$phenotype[match(cc$cluster, prof$cluster)]
            cc <- cc[order(pmin(cc$q_absolute, cc$q_share, na.rm = TRUE)), , drop = FALSE]
            wcsv(cc, sprintf("explore_cluster_count_concordance%s.csv", tag))
            only_abs <- cc$cluster[cc$verdict == "absolute only"]
            if (length(only_abs))
              log_msg("  ", length(only_abs), " cluster(s) significant on cell ",
                      "number but NOT on share -- the case a frequency table ",
                      "cannot express: ", paste(utils::head(only_abs, 8), collapse = ", "))
          }
        }
      } else if (!is.null(total_counts)) {
        note("absolute_group_statistics", "not run",
             "no embedded sample carried an external total")
      }
    } else {
      note("group_statistics", "not run",
           "no --group-column, or fewer than two group levels")
    }

    if (!is.null(pops)) {
      lab <- rep(NA_character_, nrow(cells))
      for (s in unique(cells$sample_id)) {
        sc <- pops[[s]]$scored$labels
        if (is.null(sc)) next
        w <- which(cells$sample_id == s)
        ei <- cells$event_index[w]
        ok <- ei <= length(sc)
        lab[w[ok]] <- sc[ei[ok]]
      }
      if (any(!is.na(lab))) {
        ct <- as.data.frame(table(cluster = kv, population = lab),
                            stringsAsFactors = FALSE)
        ct <- ct[ct$Freq > 0, , drop = FALSE]
        tot <- tapply(ct$Freq, ct$cluster, sum)
        ct$pct_of_cluster <- round(100 * ct$Freq / tot[ct$cluster], 2)
        ct$phenotype <- prof$phenotype[match(ct$cluster, prof$cluster)]
        names(ct)[names(ct) == "Freq"] <- "cells"
        wcsv(ct[order(ct$cluster, -ct$cells), ],
             sprintf("explore_vs_populations%s.csv", tag))

        # ---- what each cluster corresponds to --------------------------------
        # The cross-tab above is the data; on this cohort it is 189 rows of
        # long-format counts, which answers "cluster 12 is what?" only after the
        # reader has done the pivot themselves. The identity table scores the
        # match in both directions and the figure draws the correspondence
        # matrix. See R/explore-identity.R for why F1 rather than the plurality
        # label, and why the catch-all is ranked apart.
        ident <- cluster_identity_table(ct)
        if (!is.null(ident)) {
          # The subset label answers a different question from the declared
          # correspondence beside it: correspondence is measured by which cells
          # overlap, the subset by what the cluster EXPRESSES. They can
          # disagree, and a disagreement is informative -- a cluster that
          # corresponds to no declared population can still have a clean subset
          # phenotype, which is precisely the gap explore mode exists to find.
          if (!is.null(subs) && nrow(subs)) {
            .m <- match(ident$cluster, subs$cluster)
            ident$subset_label   <- subs$subset_label[.m]
            ident$subset_markers <- subs$subset_markers[.m]
            ident$subset_margin  <- subs$margin[.m]
          }
          wcsv(ident, sprintf("explore_cluster_identity%s.csv", tag))
          fi <- file.path(ex_dir, sprintf("explore_cluster_identity%s.png", tag))
          ok <- tryCatch({
            fig_cluster_identity(ct, ident, fi, panel_label = p$name); TRUE
          }, error = function(e) {
            log_msg("  WARNING explore cluster identity figure failed: ",
                    conditionMessage(e), ", every other output is unaffected")
            FALSE
          })
          if (ok && file.exists(fi)) written <- c(written, basename(fi))
          .nd <- sum(ident$call == "undescribed by the specification")
          .nc <- sum(ident$call == "confident")
          log_msg("explore: cluster identity -- ", .nc, " of ", nrow(ident),
                  " cluster(s) match a declared population confidently (F1 >= ",
                  "0.5), ", .nd, " are >=70% catch-all and describe cells the ",
                  "specification does not name")
        }

        # THE CATCH-ALL COUNTS AS UNLABELLED. This list used to be spelled out
        # here and did not contain "Other CD45+" -- the name this package
        # actually assigns to cells inside the parent gate that match no
        # declared population. Every cluster therefore scored 0% unlabelled and
        # every verdict read "covered by the declared specification", including
        # clusters that were 100% remainder, on runs whose specification
        # described under a tenth of the cells. This table exists to report what
        # the specification is MISSING, so that failure inverted its only
        # output. is_catch_all_label() is now the single definition, shared with
        # the scoring step that assigns the label.
        undecl <- vapply(ks, function(k) {
          sub <- ct[ct$cluster == k, , drop = FALSE]
          if (!nrow(sub)) return(100)
          100 * sum(sub$cells[is_catch_all_label(sub$population)]) / sum(sub$cells)
        }, numeric(1))
        fnd <- data.frame(
          cluster = ks, phenotype = unname(pheno[ks]),
          pct_unlabelled = round(undecl, 2),
          verdict = ifelse(undecl >= 70,
            "UNDECLARED - no declared population covers most of this cluster",
            "covered by the declared specification"),
          stringsAsFactors = FALSE)
        wcsv(fnd[order(-fnd$pct_unlabelled), ], sprintf("explore_findings%s.csv", tag))

        # Blank labels only. The catch-all is deliberately KEPT here: the note
        # below explains that its spread across clusters is worth seeing in
        # explore_population_split.csv, and it is filtered out one step later
        # where it would masquerade as a finding. Dropping it at this point
        # would remove the row that table is partly written for.
        real <- ct[!is.na(ct$population) & nzchar(trimws(ct$population)), ,
                   drop = FALSE]
        if (nrow(real)) {
          # Counting DISTINCT clusters is nearly useless as a flag: a coarse
          # label like "Lymphocytes" legitimately touches most clusters, so a
          # ">= 2 clusters" rule marks almost everything and the table becomes
          # noise. What matters is CONCENTRATION -- how much of the label sits
          # in its single largest cluster. A label with 90% in one cluster is
          # coherent however many stragglers it has; one with 30% in its largest
          # is genuinely several phenotypes wearing one name.
          agg <- stats::aggregate(cells ~ population + cluster, data = real, FUN = sum)
          tot <- tapply(agg$cells, agg$population, sum)
          dom <- tapply(agg$cells, agg$population, max)
          eff <- tapply(seq_len(nrow(agg)), agg$population, function(i) {
            pr <- agg$cells[i] / sum(agg$cells[i])
            round(1 / sum(pr^2), 2)          # inverse Simpson: effective clusters
          })
          nsp <- tapply(agg$cluster, agg$population, function(x) length(unique(x)))
          pops_nm <- names(tot)
          pct_dom <- round(100 * as.numeric(dom[pops_nm]) / as.numeric(tot[pops_nm]), 1)
          sp <- data.frame(
            population = pops_nm,
            n_clusters_touched = as.integer(nsp[pops_nm]),
            effective_clusters = as.numeric(eff[pops_nm]),
            pct_in_dominant_cluster = pct_dom,
            verdict = ifelse(pct_dom < 50,
              "SPLIT - no single cluster holds half of this label; it may lump distinct phenotypes",
              "coherent - one cluster holds most of it"),
            stringsAsFactors = FALSE)
          wcsv(sp[order(sp$pct_in_dominant_cluster), ],
               sprintf("explore_population_split%s.csv", tag))

          # What explore has to say about the declared specification, in the
          # shape the declared side can use. Written back into the run
          # directory as spec_gaps.csv, but ONLY under --maybe-learn; the
          # caller decides, this only assembles it.
          # "Other CD45+" is excluded from the SPLIT verdict. It is not a
          # population, it is the complement of the specification -- the cells
          # no definition matched -- so "it spans several clusters" is true by
          # construction and says nothing a reader can act on. Left in, it
          # appears beside the real findings as though it were one of them, and
          # on a specification with poor coverage it is the single largest
          # label, so it is the row most likely to be read first.
          #
          # It stays in explore_population_split.csv, where the per-population
          # numbers are the point and the catch-all's spread is worth seeing.
          g1 <- sp[sp$pct_in_dominant_cluster < 50 &
                     !is_catch_all_label(sp$population), , drop = FALSE]
          if (nrow(g1)) gaps[[length(gaps) + 1L]] <- data.frame(
            panel = p$name,
            issue = "population spans several clusters",
            subject = g1$population,
            detail = sprintf("%.0f%% in its largest cluster; %.1f effective clusters",
                             g1$pct_in_dominant_cluster, g1$effective_clusters),
            action = paste("the declared label may lump distinct phenotypes;",
                           "its total can be flat while a subset inside it moves"),
            stringsAsFactors = FALSE)
          g2 <- fnd[fnd$pct_unlabelled >= 70, , drop = FALSE]
          if (nrow(g2)) gaps[[length(gaps) + 1L]] <- data.frame(
            panel = p$name,
            issue = "cluster no declared population covers",
            subject = g2$cluster,
            detail = sprintf("%.0f%% unlabelled; %s", g2$pct_unlabelled,
                             g2$phenotype),
            action = paste("nothing in the specification describes these cells;",
                           "see explore_suggested_spec.yaml for a starting",
                           "definition"),
            stringsAsFactors = FALSE)
        }
      }
    }

    if (!isTRUE(opt$no_explore_spec))
      written <- c(written, explore_write_spec(pos, prof, ks, ex_dir, tag))

    written <- c(written,
                 explore_figures(cells, fcols, prof, ex_dir, tag = tag,
                                 group_col = if ("group" %in% colnames(cells))
                                   "group" else NULL))

    # --total-counts figures. The QC one is drawn whenever a total was supplied,
    # including when no group exists to compare across: checking that the
    # external numbers are sane is not conditional on wanting a test.
    if (!is.null(total_counts) && nrow(total_counts)) {
      # Subset to THIS panel. The figure is written per panel, so drawing every
      # acquisition in the study on it puts samples that are not in this
      # embedding beside the ones that are, and the reader cannot tell which
      # bar belongs to the table underneath.
      tc_panel <- total_counts[total_counts$sample_id %in% sids, , drop = FALSE]
      # Both figures are drawn with the RNG stream restored afterwards. Without
      # that, drawing panel 1's figures re-rolls panel 2's embedding, and a run
      # with --total-counts would stop matching one without it.
      without_spending_draws({
        # tag is "" on a single-panel run and "_panel_1" otherwise; the label is
        # that, minus the leading underscore, so the panel reaches the image and
        # not only the filename.
        plab <- if (nzchar(tag)) sub("^_", "", tag) else ""
        fq <- file.path(ex_dir, sprintf("explore_total_counts_qc%s.png", tag))
        if (nrow(tc_panel) &&
            !is.null(try(fig_total_counts_qc(tc_panel, fq, group_of = group_of,
                                             all_samples = sids,
                                             panel_label = plab), silent = TRUE)))
          written <<- c(written, basename(fq))
        if ("cells_absolute" %in% names(ab) && any(is.finite(ab$cells_absolute))) {
          fa <- file.path(ex_dir, sprintf("explore_absolute_vs_share%s.png", tag))
          if (!is.null(try(fig_absolute_vs_share(ab, fa, group_of = group_of,
                                                 panel_label = plab),
                           silent = TRUE)))
            written <<- c(written, basename(fa))
        }
      })
    }
  }

  if (length(prov)) wcsv(do.call(rbind, prov), "explore_provenance.csv")

  # Its own report, not a section in report.html. See R/explore-report.R for
  # why: report.html is a declared deliverable and --explore must not change it.
  rp <- try(write_explore_report(ex_dir, opt), silent = TRUE)
  if (!inherits(rp, "try-error") && !is.null(rp))
    written <- c(written, basename(rp))

  log_msg("explore: ", length(unique(written)), " file(s) in ", ex_dir)
  invisible(list(dir = ex_dir, files = unique(written),
                 gaps = if (length(gaps)) do.call(rbind, gaps) else NULL))
}

#' @keywords internal
.split_opt <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  trimws(strsplit(x, ",")[[1]])
}

#' Write a draft population specification from the clusters
#'
#' Closes the loop: discovery becomes a declaration you can curate and re-run
#' through the supervised path, where it acquires per-sample thresholds,
#' uncertainty and the concordance checks. Never adopted automatically.
#'
#' @keywords internal
explore_write_spec <- function(pos, prof, ks, ex_dir, tag = "") {
  lines <- c(
    "# Draft population specification, written by cyRAVEN explore mode.",
    "#",
    "# These are SUGGESTIONS from unsupervised clusters, not a validated",
    "# specification. Each entry is one cluster, named by the markers most of",
    "# its cells were positive or negative for, judged against each sample's",
    "# own threshold where one was available.",
    "#",
    "# Curate before use: merge clusters that are one population, drop debris",
    "# and doublets, give them real names. Then run the supervised path with",
    "# --config, which is where they gain per-sample thresholds, propagated",
    "# uncertainty and the six specification checks.",
    "populations:")
  for (k in ks) {
    r <- pos[k, ]
    up <- names(r)[is.finite(r) & r >= 0.65]
    dn <- names(r)[is.finite(r) & r <= 0.20]
    up <- utils::head(up[order(-r[up])], 4L)
    dn <- utils::head(dn[order(r[dn])], 3L)
    if (!length(up)) next
    lines <- c(lines,
               sprintf("  # %s, %.2f%% of gated cells", k,
                       prof$pct_of_gated[match(k, prof$cluster)]),
               sprintf("  %s:", k),
               # "above"/"below", the vocabulary score_populations() actually
               # reads. This emitted "pos"/"neg" until now, which no parser in
               # the package accepts -- so the file the header invites you to
               # curate and pass to --config could never have been run, and
               # would have reported every population UNAVAILABLE without an
               # error. The whole point of the file is that it is a starting
               # specification, so it has to parse.
               sprintf("    \"%s\": above", up),
               if (length(dn)) sprintf("    \"%s\": below", dn))
  }
  nm <- sprintf("explore_suggested_spec%s.yaml", tag)
  writeLines(lines, file.path(ex_dir, nm))
  nm
}
