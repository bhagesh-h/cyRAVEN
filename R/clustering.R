# =============================================================================
# cyRAVEN -- unsupervised clustering and embedding-model persistence
#   * unsupervised clustering (FlowSOM-style SOM + consensus metaclustering)
#   * gate-vs-cluster agreement: the check that supervised gating lacks
#   * data-driven k for the reference subclustering (silhouette)
#   * persisting the UMAP model so later batches land in the same space
#
# None of this feeds back into gating: the clusters are a second, independent
# reading of the same cells, kept separate precisely so they can contradict the
# population labels.
# =============================================================================


# =============================================================================
# UNSUPERVISED CLUSTERING
# =============================================================================
#
# THE GAP THIS CLOSES, AND WHY IT IS THE BIGGEST ONE. Every population in the
# baseline is a Boolean conjunction of per-marker thresholds written down in
# advance (config_cohorts.yaml). That has real advantages -- the populations are
# named, reproducible, and mean exactly what the gating strategy document says --
# and one decisive weakness: it CANNOT FIND ANYTHING THE SPEC DOES NOT DESCRIBE.
# Cells matching no definition become "Other CD45+", which is excluded from the
# UMAPs by default, so unexpected biology is not merely unlabelled -- it is not
# drawn.
#
# It also has no way to be wrong out loud. If a threshold lands badly, the
# affected population comes out the wrong size and there is nothing to contradict
# it: a mis-set CD4 cut and a genuine absence of CD4 T cells produce identical
# output. Unsupervised clustering is what breaks that tie, because a CD4 T-cell
# cluster appears at its true size whether or not the CD4 threshold was set well,
# and the DISAGREEMENT between cluster and gate is the diagnosis (see
# cluster_gate_agreement() below).
#
# THE ALGORITHM, AND WHY THIS ONE. FlowSOM (Van Gassen et al. 2015) is the
# field standard and won the Weber & Robinson (2016) benchmark on both accuracy
# and runtime. It has two stages:
#   1. a self-organising map -- a grid of prototype vectors, trained so
#      neighbouring nodes hold similar cells. This is a fast, order-preserving
#      quantisation of the marker space into ~100 micro-clusters.
#   2. metaclustering -- hierarchical clustering of those 100 prototypes down to
#      the requested number of populations.
# The two-stage design is what makes it scale: the expensive step touches every
# cell once per epoch, and the clustering step operates on 100 vectors, not 10^6.
#
# WHY BOTH A PACKAGE PATH AND A BUILT-IN PATH. When the FlowSOM package is
# installed it is used, because a reference implementation is worth more than a
# re-derivation and its results are directly comparable to published work. When
# it is not, the built-in below runs the same two-stage algorithm rather than
# skipping the analysis, so this pipeline keeps its property of adding no
# mandatory dependencies. The two are algorithmically equivalent and will NOT be
# numerically identical -- different initialisation, different neighbourhood decay
# -- and `method` in the output records which one ran, because a cluster number
# from one is not a cluster number from the other.

#' Train a self-organising map on the marker matrix
#'
#' WHAT: batch-SOM training. Each epoch assigns every cell to its best-matching
#' node, then moves every node toward the weighted mean of the cells assigned
#' near it, with the neighbourhood radius shrinking over epochs.
#'
#' WHY BATCH AND NOT ONLINE: batch SOM is deterministic given the initialisation,
#' order-independent (an online SOM's result depends on the order cells arrive,
#' which here is the order files were read), and vectorises. Reproducibility of
#' cluster identity between runs is a hard requirement in this pipeline -- the
#' whole subcluster-lettering machinery exists for that reason.
#'
#' @param X numeric matrix, cells x markers, already scaled
#' @param xdim,ydim SOM grid dimensions
#' @param epochs training epochs
#' @param seed Random seed. The RNG stream is restored afterwards, so this cannot perturb later sampling. Default `42L`.
#' @return list(codes = nodes x markers prototypes, mapping = node per cell)
#' @keywords internal
som_train <- function(X, xdim = 10L, ydim = 10L, epochs = 10L, seed = 42L) {
  n <- nrow(X); p <- ncol(X)
  nodes <- xdim * ydim
  stopifnot(n >= nodes)
  grid <- cbind(rep(seq_len(xdim), times = ydim), rep(seq_len(ydim), each = xdim))
  gd2 <- as.matrix(stats::dist(grid))^2

  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  # Initialise on a random sample of real cells rather than on noise: the map
  # starts inside the data manifold, which is what lets a modest epoch count
  # converge.
  codes <- X[sample.int(n, nodes), , drop = FALSE]

  # Radius decays from covering roughly a third of the grid to a single node.
  r0 <- max(xdim, ydim) / 3
  radii <- seq(r0, 1, length.out = max(2L, epochs))
  mapping <- integer(n)
  block <- max(1L, min(4096L, n))
  for (ep in seq_len(epochs)) {
    sq_c <- rowSums(codes^2)
    # Assign every cell to its best-matching node, blocked so the cell x node
    # distance matrix never exists in full.
    for (start in seq(1L, n, by = block)) {
      idx <- start:min(n, start + block - 1L)
      d2 <- outer(rowSums(X[idx, , drop = FALSE]^2), sq_c, "+") -
            2 * (X[idx, , drop = FALSE] %*% t(codes))
      mapping[idx] <- max.col(-d2, ties.method = "first")
    }
    # Node update: neighbourhood-weighted mean of assigned cells. Gaussian kernel
    # over grid distance, the standard SOM neighbourhood.
    hood <- exp(-gd2 / (2 * radii[ep]^2))
    cnt  <- tabulate(mapping, nbins = nodes)
    # Per-node column sums in ONE pass over the data.
    #
    # The loop this replaces ran `mapping == j` and subset X for each of the 100
    # nodes, every epoch: 100 full-length logical scans plus 100 matrix
    # allocations per epoch, i.e. O(nodes * n) work and O(n) garbage, to compute
    # something rowsum() does in a single O(n) C pass. At 100k cells and 10
    # epochs that is 100 million redundant comparisons.
    #
    # rowsum() returns only the groups that OCCUR, labelled by group value, so
    # the result is placed back by name -- an empty node must stay a zero row,
    # and silently shifting rows up would corrupt the codebook.
    sums <- matrix(0, nodes, p)
    rs <- rowsum(X, mapping, reorder = FALSE)
    sums[as.integer(rownames(rs)), ] <- rs
    num <- hood %*% sums
    den <- as.vector(hood %*% cnt)
    keep <- den > 0
    codes[keep, ] <- num[keep, , drop = FALSE] / den[keep]
  }
  # Final assignment against the trained codebook.
  sq_c <- rowSums(codes^2)
  for (start in seq(1L, n, by = block)) {
    idx <- start:min(n, start + block - 1L)
    d2 <- outer(rowSums(X[idx, , drop = FALSE]^2), sq_c, "+") -
          2 * (X[idx, , drop = FALSE] %*% t(codes))
    mapping[idx] <- max.col(-d2, ties.method = "first")
  }
  list(codes = codes, mapping = mapping, xdim = xdim, ydim = ydim)
}

#' Unsupervised clustering of the embedded cells
#'
#' WHY IT CLUSTERS THE MARKER MATRIX AND NOT THE UMAP COORDINATES: UMAP is a
#' visualisation. Its distances are not metric, it does not preserve density, and
#' the gaps between its islands are partly an artefact of min_dist. Clustering it
#' would give clusters of the PICTURE. Every established tool -- FlowSOM and
#' Phenograph among them -- clusters the high-dimensional space and uses the
#' embedding only to display the result, which is what this does.
#'
#' @param cells embedding cell table
#' @param markers marker columns to cluster on (normally the UMAP feature set)
#' @param n_clusters number of metaclusters
#' @param grid SOM grid side length; grid^2 nodes
#' @param scale_method matched to run_umap()'s default so clusters and the
#'   embedding see the same geometry
#' @param epochs Number of SOM training epochs. Default `10L`.
#' @param seed Random seed. The RNG stream is restored afterwards, so this cannot perturb later sampling. Default `42L`.
#' @param max_cells Ceiling on the number of cells used. Default `200000L`.
#' @return list(cluster = integer per cell, method, codes, node)
#' @export
run_unsupervised_clusters <- function(cells, markers, n_clusters = 12L,
                                      grid = 10L, epochs = 10L,
                                      scale_method = "robust", seed = 42L,
                                      max_cells = 200000L) {
  markers <- intersect(markers, names(cells))
  markers <- markers[vapply(cells[markers], is.numeric, logical(1))]
  if (length(markers) < 2L) {
    log_msg("  unsupervised clustering skipped: fewer than 2 numeric markers")
    return(NULL)
  }
  M <- as.matrix(cells[, markers, drop = FALSE])
  ok <- stats::complete.cases(M)
  if (sum(ok) < grid * grid * 5L) {
    log_msg("  unsupervised clustering skipped: too few complete cells (",
            sum(ok), ")")
    return(NULL)
  }
  # Same scaling the embedding uses, for the reason given above: robust
  # (median/MAD) because marker distributions are bimodal and an SD is dominated
  # by the separation between the two modes rather than by within-mode spread.
  Xs <- M[ok, , drop = FALSE]
  if (scale_method == "robust") {
    med <- apply(Xs, 2, stats::median)
    mad_ <- apply(Xs, 2, function(x) { m <- stats::mad(x); if (m == 0 || !is.finite(m)) stats::sd(x) else m })
    Xs <- sweep(sweep(Xs, 2, med, "-"), 2, mad_, "/")
  } else if (scale_method == "zscore") {
    Xs <- scale(Xs)
  }
  Xs[!is.finite(Xs)] <- 0

  cl <- rep(NA_integer_, nrow(cells))
  node <- rep(NA_integer_, nrow(cells))
  method <- NA_character_
  codes <- NULL

  # Seed and stream guard for BOTH branches.
  #
  # WHY IT IS HERE AND NOT ONLY IN som_train(). The guard used to sit inside
  # som_train(), which covers the built-in path and nothing else. FlowSOM's
  # metaClustering_consensus() takes a seed, but BuildSOM() does not: it draws
  # its node initialisation from the global stream. So with the FlowSOM package
  # installed, this function was neither reproducible (two calls with the same
  # `seed` argument saw different stream states and built different maps) nor
  # stream-safe (it left the stream advanced, which silently changes which cells
  # every later stage subsamples). Neither showed up while FlowSOM was absent
  # from the test environment, and both are exactly what a Suggests-dependent
  # code path can hide.
  #
  # Setting the seed here makes the two branches behave the same way, which is
  # the contract the `seed` argument already advertised.
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  if (requireNamespace("FlowSOM", quietly = TRUE)) {
    res <- try({
      fs <- FlowSOM::ReadInput(Xs, transform = FALSE, scale = FALSE)
      fs <- FlowSOM::BuildSOM(fs, colsToUse = seq_len(ncol(Xs)),
                              xdim = grid, ydim = grid, silent = TRUE)
      mc <- FlowSOM::metaClustering_consensus(fs$map$codes, k = n_clusters, seed = seed)
      list(node = fs$map$mapping[, 1], meta = as.integer(mc), codes = fs$map$codes)
    }, silent = TRUE)
    if (!inherits(res, "try-error")) {
      node[ok] <- res$node
      cl[ok] <- res$meta[res$node]
      codes <- res$codes
      method <- paste0("FlowSOM ", utils::packageVersion("FlowSOM"),
                       " (SOM ", grid, "x", grid, " + consensus metaclustering)")
    }
  }
  if (is.na(method)) {
    som <- som_train(Xs, xdim = grid, ydim = grid, epochs = epochs, seed = seed)
    # Metacluster the codebook by hierarchical clustering. Ward's criterion on
    # Euclidean distance is what FlowSOM's consensus step converges to in
    # practice and, unlike k-means on the codes, is deterministic -- which matters
    # here because a cluster number has to mean the same thing between runs.
    hc <- stats::hclust(stats::dist(som$codes), method = "ward.D2")
    meta <- stats::cutree(hc, k = min(n_clusters, nrow(som$codes)))
    node[ok] <- som$mapping
    cl[ok] <- meta[som$mapping]
    codes <- som$codes
    method <- paste0("built-in batch SOM ", grid, "x", grid,
                     " + Ward.D2 metaclustering (FlowSOM package not installed)")
  }

  # Renumber clusters by descending size so cluster 1 is always the largest.
  # Arbitrary label order would otherwise permute between runs and make
  # "cluster 7" meaningless across two results folders.
  tb <- sort(table(cl[!is.na(cl)]), decreasing = TRUE)
  remap <- setNames(seq_along(tb), names(tb))
  cl[!is.na(cl)] <- unname(remap[as.character(cl[!is.na(cl)])])

  log_msg("  unsupervised clustering: ", method, " -> ",
          length(unique(cl[!is.na(cl)])), " clusters over ", sum(ok), " cells")
  list(cluster = cl, node = node, method = method, codes = codes,
       markers = markers, n_clusters = n_clusters)
}


# =============================================================================
# GATE vs CLUSTER AGREEMENT: the cross-check supervised gating lacks
# =============================================================================
#
# WHAT IT ANSWERS. For every unsupervised cluster: which gate label do its cells
# carry, and how pure is that? And for every gate label: is it recovered as a
# cluster, or is it scattered across several?
#
# HOW TO READ IT, which is the part that matters:
#
#   A cluster that is >80% one gate label      the gate and the data agree.
#   A cluster that is mostly "Other CD45+"     a real population the spec does
#                                              not describe. This is the finding
#                                              supervised gating structurally
#                                              cannot produce.
#   A gate label split across many clusters    the label is a union of distinct
#                                              phenotypes; the gate is coarser
#                                              than the biology.
#   A gate label with far fewer cells than the
#   cluster it dominates                       THE THRESHOLD IS WRONG. The cells
#                                              are there and cluster together;
#                                              the Boolean rule is rejecting
#                                              them. This is the diagnosis that
#                                              a frequency table alone cannot
#                                              deliver, and the specific reason
#                                              this function exists.
#
# The last case is worth stating concretely, because it is what motivated this:
# if CD4 T cells score 0.29% of cells while a cluster of ~20% of cells is
# CD4-bright and labelled "Other CD45+", the population is not absent -- the cut
# is misplaced, and the two numbers side by side say so immediately.

#' Cross-tabulate unsupervised clusters against gate labels
#'
#' @param cells embedding cell table with population_label
#' @param cluster integer vector from run_unsupervised_clusters()
#' @param other_pattern regex identifying the catch-all label
#' @return list(per_cluster, per_population)
#' @export
cluster_gate_agreement <- function(cells, cluster,
                                   other_pattern = "^Other|unclassified") {
  if (is.null(cluster) || !"population_label" %in% names(cells)) return(NULL)
  lab <- as.character(cells$population_label)
  ok <- !is.na(cluster) & !is.na(lab) & nzchar(trimws(lab))
  if (!any(ok)) return(NULL)
  cl <- cluster[ok]; lb <- lab[ok]
  tb <- table(cl, lb)
  if (!nrow(tb)) return(NULL)

  n_by_lab <- table(lb)
  per_cluster <- do.call(rbind, lapply(rownames(tb), function(k) {
    row <- tb[k, ]
    tot <- sum(row)
    o <- sort(row, decreasing = TRUE)
    dom <- names(o)[1]
    purity <- as.numeric(o[1]) / tot
    second <- if (length(o) > 1L) names(o)[2] else NA_character_
    data.frame(
      cluster = as.integer(k), n_cells = as.integer(tot),
      pct_of_all_cells = round(100 * tot / sum(tb), 2),
      dominant_gate_label = dom, purity_pct = round(100 * purity, 1),
      second_label = second,
      second_pct = if (length(o) > 1L) round(100 * as.numeric(o[2]) / tot, 1) else NA_real_,
      n_labels_present = sum(row > 0),
      # Share of the DOMINANT label's total cells that this cluster holds. Low
      # values mean the label is fragmented across clusters.
      pct_of_that_label_captured =
        round(100 * as.numeric(o[1]) / as.numeric(n_by_lab[dom]), 1),
      interpretation =
        if (grepl(other_pattern, dom))
          "UNDESCRIBED POPULATION - the gate spec has no definition for these cells"
        else if (purity >= 0.8) "clean - cluster and gate agree"
        else if (purity >= 0.5) "mixed - gate label covers more than one phenotype here"
        else "fragmented - no gate label dominates this cluster",
      stringsAsFactors = FALSE)
  }))

  per_population <- do.call(rbind, lapply(colnames(tb), function(p) {
    col <- tb[, p]
    tot <- sum(col)
    o <- sort(col, decreasing = TRUE)
    # Effective number of clusters the label occupies (inverse Simpson): 1 means
    # it lands in one cluster, 4 means it is genuinely spread over about four.
    prop <- col / tot
    eff <- 1 / sum(prop^2)
    dom_cl <- as.integer(names(o)[1])
    dom_cl_size <- sum(tb[as.character(dom_cl), ])
    data.frame(
      population = p, n_cells = as.integer(tot),
      dominant_cluster = dom_cl,
      pct_in_dominant_cluster = round(100 * as.numeric(o[1]) / tot, 1),
      effective_n_clusters = round(eff, 2),
      dominant_cluster_total_cells = as.integer(dom_cl_size),
      # The threshold-failure signal: the label holds a small share of a cluster
      # it nonetheless dominates, i.e. most of that cluster's cells failed the
      # Boolean rule while sitting in the same phenotypic neighbourhood.
      pct_of_dominant_cluster_captured =
        round(100 * as.numeric(o[1]) / dom_cl_size, 1),
      interpretation =
        if (eff >= 3) "fragmented - this label spans several distinct phenotypes"
        else if (100 * as.numeric(o[1]) / dom_cl_size < 40 && dom_cl_size > tot)
          "SUSPECT THRESHOLD - most cells of this cluster are phenotypically alike but fail the gate"
        else "coherent",
      stringsAsFactors = FALSE)
  }))

  list(per_cluster = per_cluster[order(-per_cluster$n_cells), , drop = FALSE],
       per_population = per_population[order(-per_population$n_cells), , drop = FALSE])
}

#' Unsupervised-clustering figure: the embedding by cluster, next to the same
#' embedding by gate label
#'
#' WHY SIDE BY SIDE AND NOT TWO FILES: the entire value is the comparison. A
#' cluster panel alone is a pretty picture; against the gate panel it is an
#' audit, and a reader can see in one glance which islands the spec named and
#' which it missed.
#' @param cells Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`, `sample_id`, `population_label` and one column per marker.
#' @param cluster Integer vector of cluster assignments, one per cell.
#' @param outfile Path to write the figure to.
#' @param agreement The agreement.
#' @param panel_label Marker-panel name added to the figure title; empty for none. Default `""`.
#' @param dpi Resolution in dots per inch. May be reduced automatically to respect the raster ceiling; see [safe_ggsave()]. Default `200`.
#' @param colors Named list of colours; defaults to the package palette. See [fcs_colors()]. Default `fcs_colors()`.
#' @export
fig_unsupervised_clusters <- function(cells, cluster, outfile, agreement = NULL,
                                      panel_label = "", dpi = 200,
                                      colors = fcs_colors()) {
  if (is.null(cluster) || !all(c("umap_1", "umap_2") %in% names(cells)))
    return(invisible(NULL))
  d <- cells
  d$.cluster <- factor(cluster)
  d <- d[!is.na(d$.cluster) & is.finite(d$umap_1) & is.finite(d$umap_2), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  aes_pt <- auto_point_aes(nrow(d))

  lv <- levels(droplevels(d$.cluster))
  ccols <- setNames(pop_palette(length(lv), colors = colors), lv)

  # Cluster number printed at each cluster's median position, matching the way
  # the multigraph overlay labels its clusters -- plain text with a halo, no box.
  cc <- do.call(rbind, lapply(lv, function(k) {
    s <- d$.cluster == k
    data.frame(.cluster = k, x = stats::median(d$umap_1[s]),
               y = stats::median(d$umap_2[s]), stringsAsFactors = FALSE)
  }))

  p1 <- ggplot(d, aes(umap_1, umap_2, colour = .cluster)) +
    geom_point(size = aes_pt$size, alpha = aes_pt$alpha, stroke = 0) +
    scale_colour_manual(values = ccols, name = "cluster") +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 2))
  for (dx in c(-0.08, 0.08)) for (dy in c(-0.08, 0.08))
    p1 <- p1 + annotate("text", x = cc$x + dx, y = cc$y + dy, label = cc$.cluster,
                        size = 3.1, colour = "white", fontface = "bold")
  p1 <- p1 +
    annotate("text", x = cc$x, y = cc$y, label = cc$.cluster, size = 3.1,
             colour = colors$bracket, fontface = "bold") +
    labs(x = "UMAP 1", y = "UMAP 2", title = "Unsupervised clusters") +
    theme_cyto(colors = colors) +
    theme(legend.position = "right", legend.key.size = unit(9, "pt"),
          legend.text = element_text(size = 7))

  p2 <- if ("population_label" %in% names(d)) {
    lv2 <- sort(unique(as.character(d$population_label)))
    ggplot(d, aes(umap_1, umap_2, colour = population_label)) +
      geom_point(size = aes_pt$size, alpha = aes_pt$alpha, stroke = 0) +
      scale_colour_manual(values = population_colours(lv2, colors = colors),
                          name = "gate label") +
      guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 1)) +
      labs(x = "UMAP 1", y = "UMAP 2", title = "Gate labels (config spec)") +
      theme_cyto(colors = colors) +
      theme(legend.position = "right", legend.key.size = unit(9, "pt"),
            legend.text = element_text(size = 7))
  } else NULL

  cap <- paste0(
    "Clusters are computed on the marker matrix, not on these coordinates \u2014 UMAP is ",
    "the display, not the input.\nA cluster with no matching gate label is a population ",
    "the config does not describe; a gate label scattered across clusters is a label ",
    "covering several phenotypes.")
  if (!is.null(agreement)) {
    sus <- agreement$per_cluster$cluster[
      grepl("^UNDESCRIBED", agreement$per_cluster$interpretation)]
    bad <- agreement$per_population$population[
      grepl("^SUSPECT", agreement$per_population$interpretation)]
    if (length(sus))
      cap <- paste0(cap, "\nUndescribed cluster(s): ", paste(sus, collapse = ", "), ".")
    if (length(bad))
      cap <- paste0(cap, "\nPossible threshold problem for: ", paste(bad, collapse = ", "),
                    " \u2014 see cluster_gate_agreement_populations.csv.")
  }

  plots <- if (is.null(p2)) list(p1) else list(p1, p2)
  fig <- patchwork::wrap_plots(plots, ncol = length(plots)) +
    patchwork::plot_annotation(
      title = paste0("Unsupervised clustering vs. the gate spec",
                     if (nzchar(panel_label)) paste0(" \u2014 ", panel_label) else ""),
      caption = cap,
      theme = ggplot2::theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.caption = element_text(size = 7, hjust = 0, colour = colors$caption_text)))
  safe_ggsave(outfile, plot = fig, width = 6.6 * length(plots), height = 5.4,
              dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile, " (", length(lv), " clusters)")
  invisible(fig)
}


# =============================================================================
# DATA-DRIVEN k FOR THE REFERENCE SUBCLUSTERING
# =============================================================================
#
# THE GAP THIS CLOSES. subcluster_by_reference() fixes k = 3 for every
# population. Nothing checks whether three compartments is what the data
# contains, and every subcluster label in umap_multigraph_overlay.png and every
# row of subcluster_marker_shifts.csv rests on that number. A population with two
# real compartments gets one of them split arbitrarily; a population with five
# gets four of them merged.
#
# THE METHOD. Mean silhouette width over a range of k, computed on the reference
# group's cells only -- the same cells the subclustering is fitted to. Silhouette
# compares each cell's mean distance to its own cluster against its mean distance
# to the nearest other cluster, so it rewards partitions that are both tight and
# separated, and needs no ground truth.
#
# HONEST LIMITS, because a "chosen" k invites more trust than a fixed one:
# silhouette systematically favours small k and compact, roughly spherical
# clusters, which is also what k-means favours -- so the criterion and the
# algorithm share a bias. It is a better default than a hardcoded 3, not an
# oracle. The chosen k and the full score curve are both written out so the
# choice is inspectable rather than implicit, and --subcluster-k still overrides
# it outright.

#' Mean silhouette width for a k-means partition
#'
#' Computed on a subsample: silhouette is O(n^2) in distances, and the mean over
#' a few thousand cells is stable to well within the differences between adjacent
#' k that this is used to resolve.
#' @param X Numeric matrix, cells x features, already scaled.
#' @param k Neighbourhood size, or number of clusters, depending on the function.
#' @param sample_n Number of cells to subsample. Default `1500L`.
#' @param seed Random seed. The RNG stream is restored afterwards, so this cannot perturb later sampling. Default `42L`.
#' @param nstart Number of k-means restarts. Default `5L`.
#' @keywords internal
mean_silhouette <- function(X, k, sample_n = 1500L, seed = 42L, nstart = 5L) {
  n <- nrow(X)
  if (n < k * 10L || k < 2L) return(NA_real_)
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  km <- try(stats::kmeans(X, centers = k, nstart = nstart, iter.max = 50L), silent = TRUE)
  if (inherits(km, "try-error")) return(NA_real_)
  idx <- if (n > sample_n) sort(sample.int(n, sample_n)) else seq_len(n)
  Xs <- X[idx, , drop = FALSE]; cl <- km$cluster[idx]
  if (length(unique(cl)) < 2L) return(NA_real_)
  D <- as.matrix(stats::dist(Xs))
  sil <- vapply(seq_along(cl), function(i) {
    own <- cl == cl[i]
    if (sum(own) <= 1L) return(0)
    a <- mean(D[i, own & seq_along(cl) != i])
    b <- min(vapply(setdiff(unique(cl), cl[i]), function(g)
      mean(D[i, cl == g]), numeric(1)))
    if (!is.finite(a) || !is.finite(b) || max(a, b) == 0) return(0)
    (b - a) / max(a, b)
  }, numeric(1))
  mean(sil, na.rm = TRUE)
}

#' Choose k per population from the reference group's cells
#'
#' @param cells Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`, `sample_id`, `population_label` and one column per marker.
#' @param markers Character vector of marker names to use.
#' @param group_col Name of the column holding the biological grouping (cohort). Default `"cohort"`.
#' @param reference The group every other group is compared against.
#' @param k_range Candidate values of k to score. Default `2:5`.
#' @param min_ref Minimum reference-group cells required to fit. Default `150L`.
#' @param seed Random seed. The RNG stream is restored afterwards, so this cannot perturb later sampling. Default `42L`.
#' @return list(k = named integer population -> k, curve = data.frame of scores)
#' @export
choose_subcluster_k <- function(cells, markers, group_col = "cohort",
                                reference = NULL, k_range = 2:5,
                                min_ref = 150L, seed = 42L) {
  markers <- intersect(markers, names(cells))
  markers <- markers[vapply(cells[markers], is.numeric, logical(1))]
  if (length(markers) < 2L || is.null(reference)) return(NULL)
  if (!all(c(group_col, "population_label") %in% names(cells))) return(NULL)
  grp <- as.character(cells[[group_col]])
  if (!any(grp == reference, na.rm = TRUE)) return(NULL)
  M <- as.matrix(cells[, markers, drop = FALSE])

  ks <- integer(0); curve <- list()
  for (pp in sort(unique(as.character(cells$population_label)))) {
    r <- which(!is.na(cells$population_label) & cells$population_label == pp &
               !is.na(grp) & grp == reference)
    Xr <- M[r, , drop = FALSE]
    Xr <- Xr[stats::complete.cases(Xr), , drop = FALSE]
    if (nrow(Xr) < min_ref) next
    kk <- k_range[k_range <= floor(nrow(Xr) / 75)]
    if (!length(kk)) next
    sc <- vapply(kk, function(k) mean_silhouette(Xr, k, seed = seed), numeric(1))
    if (all(!is.finite(sc))) next
    best <- kk[which.max(sc)]
    ks[pp] <- as.integer(best)
    curve[[length(curve) + 1L]] <- data.frame(
      population = pp, k = kk, mean_silhouette = round(sc, 4),
      chosen = kk == best, n_reference_cells = nrow(Xr), stringsAsFactors = FALSE)
  }
  if (!length(ks)) return(NULL)
  list(k = ks, curve = do.call(rbind, curve))
}


# =============================================================================
# PERSIST THE UMAP MODEL
# =============================================================================
#
# THE GAP THIS CLOSES. run_umap() discards uwot's model, so adding one sample
# re-embeds everything and every coordinate moves. Cluster 4 in this run is not
# cluster 4 in the next, and two results folders from the same study cannot be
# laid side by side. The README documents adding batches over time, which makes
# this a live problem rather than a hypothetical one.
#
# The established remedy is to train once with uwot's ret_model = TRUE, then
# project later data through umap_transform() into the SAME space.
#
# WHAT PROJECTION IS AND IS NOT. Projected cells are placed by the existing
# model; they do not move the manifold. That is the point -- coordinates stay
# comparable -- and it is also the limitation: a population present only in the
# new samples has no region of its own to land in and will be placed among
# whatever it is nearest. Projection is right for adding more of the same kind of
# sample, and wrong for adding a new panel or a new cell type. Retrain when the
# biology changes; project when the cohort grows.
#
# WHY save_uwot() AND NOT saveRDS(): a uwot model holds an external pointer to
# its nearest-neighbour index, which saveRDS serialises as a null pointer that
# fails on load -- usually at transform time, well after the file looked fine.
# save_uwot() bundles the index properly. saveRDS is used only as a labelled
# fallback on uwot versions that lack it.

#' Save a trained UMAP model plus the scaling needed to reproduce its input
#'
#' The scaling parameters travel WITH the model. Re-deriving median/MAD from a
#' new batch would scale it against itself, so the projection would land in a
#' space subtly different from the one the model was trained in -- a silent error
#' that looks like a batch effect.
#' @param model A trained uwot model.
#' @param path File path.
#' @param scale_params Scaling constants that must accompany the model for a projection to land in the same space.
#' @param features Character vector of feature names the model was trained on.
#' @param meta Named list of metadata stored alongside the model. Default `list()`.
#' @export
save_umap_model <- function(model, path, scale_params, features, meta = list()) {
  if (is.null(model)) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  side <- paste0(path, ".meta.rds")
  saveRDS(list(scale_params = scale_params, features = features, meta = meta,
               uwot_version = as.character(utils::packageVersion("uwot")),
               saved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), side)
  ok <- FALSE
  if ("save_uwot" %in% getNamespaceExports("uwot")) {
    r <- try(uwot::save_uwot(model, file = path, verbose = FALSE), silent = TRUE)
    ok <- !inherits(r, "try-error")
  }
  if (!ok) {
    saveRDS(model, paste0(path, ".rds"))
    log_msg("  NOTE uwot::save_uwot() unavailable \u2014 model written with saveRDS to ",
            basename(paste0(path, ".rds")), ". It may not reload on another ",
            "machine or uwot version; retrain rather than trusting a failed load.")
  }
  log_msg("  saved UMAP model to ", basename(path), " (+ .meta.rds)")
  invisible(path)
}

#' Load a saved UMAP model
#' @param path File path.
#' @export
load_umap_model <- function(path) {
  side <- paste0(path, ".meta.rds")
  if (!file.exists(side)) {
    log_msg("  NOTE no ", basename(side), " beside the model; cannot reproduce its ",
            "scaling, so projection is refused rather than done wrongly")
    return(NULL)
  }
  meta <- readRDS(side)
  m <- NULL
  if (file.exists(path) && "load_uwot" %in% getNamespaceExports("uwot"))
    m <- try(uwot::load_uwot(file = path, verbose = FALSE), silent = TRUE)
  if ((is.null(m) || inherits(m, "try-error")) && file.exists(paste0(path, ".rds")))
    m <- try(readRDS(paste0(path, ".rds")), silent = TRUE)
  if (is.null(m) || inherits(m, "try-error")) {
    log_msg("  NOTE UMAP model at ", path, " could not be loaded \u2014 embedding fresh")
    return(NULL)
  }
  list(model = m, scale_params = meta$scale_params, features = meta$features,
       meta = meta$meta)
}

#' Project new cells into a saved embedding
#'
#' Refuses rather than improvises when the feature sets differ: a model trained
#' on 11 markers cannot place cells described by 9, and quietly filling the gap
#' with zeros would produce coordinates that look plausible and mean nothing.
#' @param saved A model loaded by [load_umap_model()].
#' @param mat Numeric matrix, cells x features.
#' @param n_threads Number of threads. Default `1L`.
#' @export
project_umap <- function(saved, mat, n_threads = 1L) {
  if (is.null(saved)) return(NULL)
  miss <- setdiff(saved$features, colnames(mat))
  if (length(miss)) {
    log_msg("  NOTE projection refused: the saved model needs marker(s) absent ",
            "from this run (", paste(miss, collapse = ", "), "). Embedding fresh.")
    return(NULL)
  }
  X <- mat[, saved$features, drop = FALSE]
  sp <- saved$scale_params
  if (!is.null(sp) && !is.null(sp$center))
    X <- sweep(sweep(X, 2, sp$center, "-"), 2, sp$scale, "/")
  X[!is.finite(X)] <- 0
  emb <- try(uwot::umap_transform(X, saved$model, n_threads = n_threads), silent = TRUE)
  if (inherits(emb, "try-error")) {
    log_msg("  NOTE umap_transform() failed (", conditionMessage(attr(emb, "condition")),
            ") \u2014 embedding fresh")
    return(NULL)
  }
  colnames(emb) <- c("umap_1", "umap_2")
  emb
}
