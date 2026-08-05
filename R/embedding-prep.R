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
