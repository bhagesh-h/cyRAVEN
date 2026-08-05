# 3. EMBEDDING
# =============================================================================

#' Build the UMAP input matrix and run the embedding
#'
#' WHAT: assembles an asinh-transformed, per-feature-scaled matrix from the
#'       per-sample subsamples and returns UMAP coordinates.
#' WHY:  scaling is computed on the GATED, subsampled matrix only. The template
#'       z-scored a debris-dominated matrix, so its feature means and SDs
#'       described noise rather than cell populations. Robust scaling (median /
#'       MAD) is used by default because marker distributions are bimodal, and
#'       an outlier-sensitive mean/SD compresses the informative separation.
#'
#' @param mat numeric matrix, cells x features, already asinh-transformed.
#' @param scale_method "robust" (median/MAD), "zscore", or "none".
#' @param n_neighbors,min_dist,metric,n_epochs uwot parameters.
#' @param seed RNG seed for reproducible coordinates.
#' @param n_threads threads for uwot.
#' @return list(coords = 2-column matrix, params = list of settings used)
#' @param ret_model if TRUE, ask uwot to keep the trained model so later batches
#'   can be PROJECTED into this same embedding instead of re-embedded from
#'   scratch (see save_umap_model()/project_umap() below).
#'
#'   DEFAULT FALSE, AND THE DEFAULT PATH IS UNCHANGED from the baseline -- same
#'   call, same arguments, same coordinates.
#'
#'   TURNING IT ON CHANGES THE COORDINATES. Measured on 3,000 x 6 synthetic
#'   cells at a fixed seed: each setting is perfectly reproducible with itself,
#'   but ret_model = TRUE and ret_model = FALSE give different embeddings of the
#'   same data (max coordinate difference 3.7, per-axis correlation 0.96). uwot
#'   takes a different internal path when it has to build and retain the
#'   nearest-neighbour index, and that path consumes the RNG stream differently.
#'   The embedding is equally valid -- it is not a degraded one -- but it is not
#'   the SAME picture, so clusters may be numbered differently and islands may
#'   sit elsewhere on the page.
#'
#'   CONSEQUENCE FOR --save-umap-model: adding that flag to a finished analysis
#'   re-draws the UMAP. Decide to persist the model at the START of a study, not
#'   after the figures have been circulated.
#' @export
run_umap <- function(mat, scale_method = "robust", n_neighbors = 30L,
                     min_dist = 0.3, metric = "euclidean", n_epochs = 200L,
                     seed = 42L, n_threads = max(1L, parallel::detectCores() - 1L),
                     ret_model = FALSE) {
  stopifnot(is.matrix(mat), ncol(mat) >= 2L)
  keep_col <- apply(mat, 2, function(x) is.finite(sd(x)) && sd(x) > 0)
  if (any(!keep_col)) {
    message("[umap] dropping zero-variance features: ",
            paste(colnames(mat)[!keep_col], collapse = ", "))
    mat <- mat[, keep_col, drop = FALSE]
  }
  # The scaling constants are CAPTURED, not just applied. A model is only usable
  # for projection together with the transform that produced its input: deriving
  # median/MAD afresh from a later batch would scale that batch against itself
  # and land it in a subtly different space, which shows up as a batch effect
  # that is really an arithmetic error.
  scale_params <- NULL
  if (scale_method == "robust") {
    med <- apply(mat, 2, median)
    mad_ <- apply(mat, 2, function(x) { m <- mad(x); if (m == 0 || !is.finite(m)) sd(x) else m })
    mat <- sweep(sweep(mat, 2, med, "-"), 2, mad_, "/")
    scale_params <- list(center = med, scale = mad_, method = "robust")
  } else if (scale_method == "zscore") {
    mat <- scale(mat)
    scale_params <- list(center = attr(mat, "scaled:center"),
                         scale = attr(mat, "scaled:scale"), method = "zscore")
  }
  mat[!is.finite(mat)] <- 0
  n_neighbors <- min(n_neighbors, max(2L, nrow(mat) - 1L))
  set.seed(seed)
  res <- uwot::umap(mat, n_neighbors = n_neighbors, min_dist = min_dist,
                    metric = metric, n_epochs = n_epochs, n_components = 2L,
                    verbose = FALSE, n_threads = n_threads,
                    ret_model = isTRUE(ret_model))
  # With ret_model the return is a list whose $embedding holds the coordinates;
  # without it, it is the coordinate matrix itself. Unwrap so every caller sees
  # the same `coords` it always did.
  emb <- if (isTRUE(ret_model)) res$embedding else res
  colnames(emb) <- c("umap_1", "umap_2")
  list(coords = emb,
       model = if (isTRUE(ret_model)) res else NULL,
       scale_params = scale_params,
       params = list(n_cells = nrow(mat), n_features = ncol(mat),
                     features = colnames(mat), scale_method = scale_method,
                     n_neighbors = n_neighbors, min_dist = min_dist,
                     metric = metric, n_epochs = n_epochs, seed = seed))
}

# =============================================================================
