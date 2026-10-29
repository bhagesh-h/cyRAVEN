# 3. EMBEDDING
# =============================================================================

#' The marker scaling the embedding and the clustering use
#'
#' WHAT: `--umap-scale` when it names a method. Under `auto`, the default: no
#' scaling where the transform already puts every marker on one scale (logicle,
#' or arcsinh with per-marker cofactors), and the percentile range otherwise.
#'
#' WHY NO SCALING ON A SHARED SCALE. Logicle puts every channel on the same
#' 4.5-decade display scale by construction, and per-marker cofactors give every
#' marker's negative population the same width, so the markers are already
#' comparable. Dividing each by its spread (MAD) then only reweights them by how
#' narrow their negative population is: a mostly negative marker such as Vd1 has
#' a small MAD and gets its noise multiplied. Measured on 40,000 cells of one
#' cohort, embedding unscaled logicle values kept 97.0% of each cell's 15
#' nearest neighbours in its own cluster; median/MAD scaling kept 94.8%, and
#' merged islands that unscaled data keeps apart. FlowSOM itself defaults to no
#' scaling.
#'
#' WHY THE PERCENTILE RANGE, NOT MEDIAN/MAD, UNDER ONE PANEL COFACTOR. One
#' cofactor does not put markers on a common scale, so some scaling is needed;
#' but the MAD of a mostly negative marker is the width of its negative
#' population, so median/MAD scaling inflates exactly the background and
#' spreading noise that distort embeddings (Mazza et al. 2018), and does so
#' more the rarer the positive cells are. Dividing by the 0.5th-99.5th
#' percentile range instead gives every marker the same span from its negatives
#' to its brightest cells, whatever share is positive. `robust` remains
#' available by name.
#'
#' @param opt The run's option list.
#' @return "none", "range", "robust" or "zscore".
#' @keywords internal
umap_scale_method <- function(opt) {
  s <- opt$umap_scale %||% "auto"
  if (!identical(s, "auto")) return(match.arg(s, c("none", "range", "robust", "zscore")))
  if (identical(opt$transform %||% "arcsinh", "logicle") ||
      isTRUE(opt$per_marker_cofactor)) "none" else "range"
}

#' Scale each feature before embedding or clustering
#'
#' One implementation for [run_umap()] and [run_unsupervised_clusters()], so
#' clusters and coordinates are computed in the same space.
#'
#' @param mat numeric matrix, cells x features.
#' @param method "none", "range" (0.5th-99.5th percentile range), "robust"
#'   (median/MAD) or "zscore".
#' @return list(mat, params); `params` carries `center` and `scale` so a saved
#'   model can be applied to later data in the same space, or NULL for "none".
#' @keywords internal
scale_features <- function(mat, method = "range") {
  spread <- function(x, f) { m <- f(x); if (!is.finite(m) || m == 0) stats::sd(x) else m }
  if (identical(method, "range")) {
    lo <- apply(mat, 2, stats::quantile, 0.005, names = FALSE)
    rg <- apply(mat, 2, function(x) spread(x, function(v)
      diff(stats::quantile(v, c(0.005, 0.995), names = FALSE))))
    return(list(mat = sweep(sweep(mat, 2, lo, "-"), 2, rg, "/"),
                params = list(center = lo, scale = rg, method = "range")))
  }
  if (identical(method, "robust")) {
    med <- apply(mat, 2, stats::median)
    mad_ <- apply(mat, 2, spread, stats::mad)
    return(list(mat = sweep(sweep(mat, 2, med, "-"), 2, mad_, "/"),
                params = list(center = med, scale = mad_, method = "robust")))
  }
  if (identical(method, "zscore")) {
    z <- scale(mat)
    return(list(mat = z, params = list(center = attr(z, "scaled:center"),
                                       scale = attr(z, "scaled:scale"),
                                       method = "zscore")))
  }
  list(mat = mat, params = NULL)
}

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
#' @param scale_method "range" (percentile range), "robust" (median/MAD),
#'   "zscore", or "none". See [scale_features()].
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
run_umap <- function(mat, scale_method = "range", n_neighbors = 30L,
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
  sf <- scale_features(mat, scale_method)
  mat <- sf$mat
  scale_params <- sf$params
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
