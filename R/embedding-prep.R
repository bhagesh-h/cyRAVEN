# 1. FEATURE SELECTION
# =============================================================================

#' Choose the markers that define the embedding space
#'
#' WHAT: returns the subset of available markers to use as UMAP input features.
#' WHY:  the embedding must be driven by lineage/phenotype identity, not by
#'       (a) scatter channels, which encode size/granularity and are already
#'       used for gating, (b) the viability dye, whose variation is technical,
#'       (c) height/width channels, which are redundant duplicates of area and
#'       would double-weight every marker, or (d) time. The original template
#'       selected features by channel-name regex and so retained 26 redundant
#'       height channels -- the single largest cause of its uninformative UMAP.
#'
#' @param marker_cols named integer vector: marker symbol -> column index
#'   (AREA channels only, as produced by the reading module).
#' @param prefer character vector of preferred lineage markers; those present
#'   are used. If none of them are present, falls back to all eligible markers.
#' @param exclude character vector of marker symbols to drop unconditionally
#'   (e.g. the detected viability marker).
#' @param lineage_only logical; if FALSE, every eligible marker is used.
#' @return character vector of marker symbols, in stable sorted order.
#' @export
select_umap_features <- function(marker_cols, prefer = NULL, exclude = NULL,
                                 lineage_only = TRUE) {
  avail <- names(marker_cols)
  # structural / technical channels are never features
  drop_pat <- "^FSC|^SSC|^Time|Time Stamp|^Event|Width$|Height$|-H$|-W$"
  eligible <- avail[!grepl(drop_pat, avail, ignore.case = TRUE)]
  if (length(exclude)) eligible <- setdiff(eligible, exclude)

  if (lineage_only && length(prefer)) {
    feats <- intersect(prefer, eligible)
    missing <- setdiff(prefer, eligible)
    if (length(missing))
      message("[features] requested lineage markers absent from this panel: ",
              paste(missing, collapse = ", "))
    if (length(feats) < 2L) {
      warning("[features] fewer than 2 preferred lineage markers present; ",
              "falling back to all ", length(eligible), " eligible markers")
      feats <- eligible
    }
  } else {
    feats <- eligible
  }
  sort(unique(feats))
}

# =============================================================================
# 2. BALANCED SUBSAMPLING
# =============================================================================

#' Draw a size-balanced subsample of gated cells across samples
#'
#' WHAT: returns row indices, per sample, capped so no sample dominates.
#' WHY:  event counts differ by orders of magnitude between files. Embedding all
#'       cells lets the largest file dictate the manifold and makes "sample"
#'       structure indistinguishable from acquisition depth. Equal-N per sample
#'       (or the smallest sample's N, whichever is smaller) makes cross-sample
#'       comparison of the shared embedding meaningful.
#'
#' @param n_per_sample named integer vector: sample_id -> number of gated cells.
#' @param cap maximum cells to take from any one sample.
#' @param equalise if TRUE, take min(cap, smallest sample's N) from every sample
#'   so all samples contribute equally; if FALSE, take min(cap, N) per sample.
#' @param total_cap optional overall ceiling on embedded cells; the per-sample
#'   allowance is reduced proportionally if the total would exceed it.
#' @return named integer vector: sample_id -> number of cells to draw.
#' @export
plan_subsample <- function(n_per_sample, cap = 20000L, equalise = TRUE,
                           total_cap = NULL) {
  n_per_sample <- n_per_sample[n_per_sample > 0]
  if (!length(n_per_sample)) stop("no sample has any gated cells to embed")
  take <- if (equalise) {
    rep(min(cap, min(n_per_sample)), length(n_per_sample))
  } else {
    pmin(cap, n_per_sample)
  }
  names(take) <- names(n_per_sample)
  if (!is.null(total_cap) && sum(take) > total_cap) {
    take <- pmax(1L, floor(take * total_cap / sum(take)))
    message("[subsample] total cap ", total_cap, " applied -> ",
            unique(take)[1], " cells/sample")
  }
  as.integer(take) -> v; names(v) <- names(take); v
}

#' Sample row indices reproducibly
#' @param gated_idx Integer indices of gated cells.
#' @param n_take Number of cells to draw.
#' @param seed Random seed. The RNG stream is restored afterwards, so this cannot perturb later sampling. Default `1L`.
#' @return integer vector of selected indices into the per-sample gated index set
#' @keywords internal
draw_subsample <- function(gated_idx, n_take, seed = 1L) {
  set.seed(seed)
  if (length(gated_idx) <= n_take) return(gated_idx)
  sort(sample(gated_idx, n_take))
}

# =============================================================================
# 2b. DENSITY-AWARE SUBSAMPLING
# =============================================================================
#
# WHY AN ALTERNATIVE TO UNIFORM DRAWING EXISTS. The uniform draw above is correct
# for the question "what does this sample look like": it reproduces the
# composition, so the embedding shows the cohort in the proportions it actually
# has. It is the wrong draw for the question `--cluster` asks.
#
# A population at 0.3% contributes about 60 cells out of 20,000. That is enough
# to be a smudge and not enough to be a cluster, so the unsupervised cross-check
# cannot recover it, and `cluster_gate_agreement()` cannot report a population
# the specification missed. README section 4.6 names finding exactly that as the
# purpose of the check, which makes the sampling the binding constraint on the
# package's own falsification claim for rare populations.
#
# WHAT THIS DOES INSTEAD. Each cell is weighted by the inverse of the local
# density of the cells around it in marker space, so cells in sparse regions are
# more likely to be drawn. This is the principle behind geometric sketching (Hie
# et al. 2019, Cell Syst 8:483) and density-dependent downsampling as SPADE uses
# it, implemented here by a distance-to-k-th-neighbour estimate rather than by
# adding a dependency.
#
# WHAT IT COSTS, AND WHY IT IS OPT-IN. The embedded set is no longer a random
# sample of the sample, so any quantity computed from `cells_umap.csv` without
# reweighting is biased toward rare populations. The frequencies, the MFI tables
# and every test are unaffected, because none of them read the embedded subset:
# they are computed per sample from all gated events. What changes is which cells
# appear in the UMAP, and therefore every embedding figure. `sampling_weight` is
# written alongside so anything derived from the embedded cells can be weighted
# back to the true composition.

#' Inverse-density sampling weights in marker space
#'
#' The weight is the distance to the k-th nearest neighbour, which is a standard
#' non-parametric density estimate: large in sparse regions, small in dense ones.
#' Distances are computed on a bounded random reference subset rather than
#' between all pairs, which makes the cost linear in cells rather than quadratic.
#'
#' @param X numeric matrix of features, one row per cell
#' @param k neighbour rank used for the density estimate
#' @param n_ref reference cells the distances are measured against
#' @param seed seed for the local RNG stream, which is restored on exit
#' @return numeric vector of non-negative weights, one per row of `X`
#' @export
density_weights <- function(X, k = 20L, n_ref = 2000L, seed = 1L) {
  if (is.null(X) || !nrow(X)) return(numeric(0))
  n <- nrow(X)
  if (n <= k + 1L) return(rep(1, n))

  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  # Scale each feature to unit MAD so one wide channel does not define the
  # geometry on its own, which is the same reason the embedding standardises.
  X <- as.matrix(X)
  s <- apply(X, 2, function(v) max(stats::mad(v, na.rm = TRUE), 1e-6))
  X <- sweep(X, 2, s, "/")
  X[!is.finite(X)] <- 0

  ref <- X[sample.int(n, min(n, as.integer(n_ref))), , drop = FALSE]
  kk <- min(as.integer(k), nrow(ref))

  # Blocked, because an n x n_ref distance matrix at 200,000 cells is 3 GB.
  w <- numeric(n)
  step <- max(1L, floor(2e6 / max(1L, nrow(ref))))
  for (from in seq(1L, n, by = step)) {
    to <- min(n, from + step - 1L)
    d <- as.matrix(stats::dist(rbind(X[from:to, , drop = FALSE], ref)))
    m <- nrow(X[from:to, , drop = FALSE])
    d <- d[seq_len(m), (m + 1L):(m + nrow(ref)), drop = FALSE]
    w[from:to] <- apply(d, 1, function(r) sort(r, partial = kk)[kk])
  }
  w[!is.finite(w) | w < 0] <- 0
  w
}

#' Draw a subsample that preserves sparse regions of marker space
#'
#' @param gated_idx Integer indices of gated cells.
#' @param X feature matrix for those cells, in the same row order
#' @param n_take Number of cells to draw.
#' @param seed Random seed. The stream is restored afterwards.
#' @param k neighbour rank for [density_weights()]
#' @return list(idx = drawn indices, weight = sampling weight per drawn cell)
#' @export
draw_subsample_rare <- function(gated_idx, X, n_take, seed = 1L, k = 20L) {
  if (length(gated_idx) <= n_take)
    return(list(idx = gated_idx, weight = rep(1, length(gated_idx))))

  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  w <- density_weights(X, k = k, seed = seed)
  if (!length(w) || all(w <= 0)) {
    idx <- sort(sample(gated_idx, n_take))
    return(list(idx = idx, weight = rep(1, length(idx))))
  }
  # A floor on the probability, so no cell is unreachable and the draw remains a
  # sample of the whole gate rather than of its outskirts only.
  p <- w + stats::quantile(w[w > 0], 0.05, names = FALSE)
  pick <- sample(seq_along(gated_idx), n_take, replace = FALSE, prob = p)
  ord <- order(gated_idx[pick])
  pick <- pick[ord]
  # The weight to reweight BY is the reciprocal of the inclusion probability,
  # normalised so an unbiased draw would give 1 throughout.
  pr <- p[pick] / sum(p)
  wt <- (1 / pr); wt <- wt / stats::median(wt)
  list(idx = gated_idx[pick], weight = wt)
}

# =============================================================================
