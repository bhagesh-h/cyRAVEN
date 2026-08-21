# =============================================================================
# Clustering, subclustering, embedding models, and RNG hygiene
# =============================================================================

test_that("unsupervised clustering is deterministic and size-ordered", {
  cells <- synth_cells()
  uc1 <- suppressMessages(
    run_unsupervised_clusters(cells, MARKERS, n_clusters = 6L, grid = 4L, seed = 1L))
  uc2 <- suppressMessages(
    run_unsupervised_clusters(cells, MARKERS, n_clusters = 6L, grid = 4L, seed = 1L))

  expect_length(uc1$cluster, nrow(cells))
  # Same seed, same answer. Cluster identity has to mean the same thing between
  # runs or none of the downstream labelling survives a re-run.
  expect_identical(uc1$cluster, uc2$cluster)
  expect_true(nzchar(uc1$method))

  # Renumbered largest-first, so "cluster 1" is comparable across result sets.
  tb <- table(uc1$cluster)
  tb <- tb[order(as.integer(names(tb)))]
  expect_true(all(diff(as.integer(tb)) <= 0))
})

test_that("the SOM's rowsum node update matches an explicit per-node loop", {
  # This is the optimisation guard. The node update was a loop of
  # colSums(X[mapping == j, ]); it is now one rowsum() pass. Empty nodes must
  # stay zero rows, rowsum() only returns groups that occur, and shifting rows
  # up would silently corrupt the codebook.
  withr::local_seed(5)
  X <- matrix(stats::rnorm(600 * 4), ncol = 4)
  mapping <- sample(c(1L, 2L, 5L, 6L), 600, TRUE)   # nodes 3 and 4 deliberately empty
  nodes <- 6L

  slow <- matrix(0, nodes, ncol(X))
  for (j in seq_len(nodes)) {
    sel <- mapping == j
    if (any(sel)) slow[j, ] <- colSums(X[sel, , drop = FALSE])
  }
  fast <- matrix(0, nodes, ncol(X))
  rs <- rowsum(X, mapping, reorder = FALSE)
  fast[as.integer(rownames(rs)), ] <- rs

  expect_equal(fast, slow, tolerance = 1e-12)
  expect_true(all(fast[3:4, ] == 0))
})

test_that("clustering leaves the RNG stream exactly as it found it", {
  # These functions run mid-pipeline. Consuming draws would silently change which
  # cells later stages subsample, which is the kind of bug that only shows up as
  # an irreproducible figure months later.
  cells <- synth_cells()
  set.seed(99); before <- stats::runif(3)
  set.seed(99)
  invisible(suppressMessages(
    run_unsupervised_clusters(cells, MARKERS, n_clusters = 4L, grid = 3L)))
  after <- stats::runif(3)
  expect_identical(before, after)
})

test_that("gate/cluster agreement identifies undescribed and suspect populations", {
  cells <- synth_cells()
  uc <- suppressMessages(
    run_unsupervised_clusters(cells, MARKERS, n_clusters = 6L, grid = 4L, seed = 1L))
  ag <- cluster_gate_agreement(cells, uc$cluster)

  expect_false(is.null(ag))
  expect_true(all(ag$per_cluster$purity_pct >= 0 & ag$per_cluster$purity_pct <= 100))
  # CD4 T cells were planted as separable, so they must land coherently.
  r <- ag$per_population[ag$per_population$population == "CD4 T cells", ]
  expect_gt(r$pct_in_dominant_cluster, 60)
  # Any cluster whose dominant label is the catch-all is a population the spec
  # does not describe, and must be called that.
  oth <- ag$per_cluster[grepl("Other", ag$per_cluster$dominant_gate_label), ]
  if (nrow(oth)) expect_true(all(grepl("^UNDESCRIBED", oth$interpretation)))
})

test_that("silhouette k-selection returns one choice per population", {
  withr::local_seed(2)
  n <- 400L
  big <- do.call(rbind, lapply(1:6, function(i) {
    pl <- rep(c("A", "B"), each = n / 2)
    d <- data.frame(sample_id = sprintf("H%d", i), cohort = "HC",
                    population_label = pl, umap_1 = stats::rnorm(n),
                    umap_2 = stats::rnorm(n), stringsAsFactors = FALSE)
    for (m in MARKERS) d[[m]] <- stats::rnorm(n)
    d
  }))
  ck <- choose_subcluster_k(big, MARKERS, group_col = "cohort",
                            reference = "HC", k_range = 2:4)
  expect_false(is.null(ck))
  expect_true(all(tapply(ck$curve$chosen, ck$curve$population, sum) == 1))
  expect_true(all(ck$k >= 2))
})

test_that("subcluster_by_reference accepts a scalar k and a named vector", {
  cells <- synth_cells()
  s1 <- subcluster_by_reference(cells, MARKERS, group_col = "cohort",
                                reference = "HC", k = 3L, min_ref = 50L)
  s2 <- subcluster_by_reference(cells, MARKERS, group_col = "cohort",
                                reference = "HC", min_ref = 50L,
                                k = c("CD4 T cells" = 2L, "B cells" = 4L))
  expect_length(s1, nrow(cells))
  expect_length(s2, nrow(cells))
  expect_true(all(s1 >= 1) && all(s2 >= 1))
})

test_that("run_umap keeps its default return shape and honours ret_model", {
  withr::local_seed(4)
  M <- matrix(stats::rnorm(600 * 5), ncol = 5,
              dimnames = list(NULL, paste0("m", 1:5)))

  a <- run_umap(M, n_neighbors = 10L, n_epochs = 20L, n_threads = 1L, seed = 3L)
  expect_true(is.matrix(a$coords))
  expect_identical(colnames(a$coords), c("umap_1", "umap_2"))
  expect_null(a$model)

  # Each setting must be reproducible with ITSELF...
  b1 <- run_umap(M, n_neighbors = 10L, n_epochs = 20L, n_threads = 1L,
                 seed = 3L, ret_model = TRUE)
  b2 <- run_umap(M, n_neighbors = 10L, n_epochs = 20L, n_threads = 1L,
                 seed = 3L, ret_model = TRUE)
  expect_equal(b1$coords, b2$coords)
  expect_false(is.null(b1$model))
  # ...and the scaling constants must travel with the model, or a projection
  # would rescale the new batch against itself and land in a different space.
  expect_true(all(c("center", "scale") %in% names(b1$scale_params)))
})

test_that("projection is refused rather than improvised when markers are missing", {
  withr::local_seed(6)
  M <- matrix(stats::rnorm(600 * 5), ncol = 5,
              dimnames = list(NULL, paste0("m", 1:5)))
  e <- run_umap(M, n_neighbors = 10L, n_epochs = 20L, n_threads = 1L,
                ret_model = TRUE)
  p <- file.path(withr::local_tempdir(), "m.uwot")
  suppressMessages(save_umap_model(e$model, p, e$scale_params, colnames(M)))
  saved <- suppressMessages(load_umap_model(p))

  expect_false(is.null(saved))
  expect_false(is.null(suppressMessages(project_umap(saved, M[1:50, , drop = FALSE],
                                                     n_threads = 1L))))
  # Filling the gap with zeros would give coordinates that look plausible and
  # mean nothing, so the refusal is the feature.
  expect_null(suppressMessages(
    project_umap(saved, M[, 1:3, drop = FALSE], n_threads = 1L)))
})
