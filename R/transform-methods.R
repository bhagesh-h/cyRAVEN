# =============================================================================
# cyRAVEN -- choosing the intensity transform
# =============================================================================
#
# WHY THERE IS A CHOICE AT ALL. Every downstream number in this package -- the
# density valley a threshold sits in, the median a differential-state test
# compares, the distance UMAP embeds -- is computed on TRANSFORMED intensity.
# The transform is therefore not a display convenience; it decides where the
# modes of a distribution fall and how far apart they are, and a threshold is
# placed between modes.
#
# arcsinh with a derived cofactor is the mass-cytometry convention and is what
# this package did exclusively. It is defensible on flow data and it has one
# concrete virtue here: the cofactor is derived from the data rather than
# assumed, so it adapts to instrument gain. But the convention in FLOW cytometry
# is the logicle (biexponential) family, because it was designed for exactly the
# problem flow data has and mass cytometry does not: compensated fluorescence
# spreads into NEGATIVE values, and a log axis cannot show them while a linear
# axis compresses the positive decades into nothing. Logicle is linear near zero
# and logarithmic away from it, so the negative population stays visible and
# on-scale instead of piling against an axis limit.
#
# Both are offered because neither is right for every panel, and a package that
# silently picks one is making a decision that belongs to whoever knows what
# instrument produced the file.
#
# WHY THE PARAMETERS ARE POOLED PER PANEL AND NOT FITTED PER FILE. Automatic
# logicle estimation is normally run per file, which gives each sample its own
# axis. Every cross-sample quantity this package computes -- a median compared
# between cohorts, a threshold checked for drift, a shared embedding -- then
# compares numbers measured on different rulers, and a difference in staining
# becomes indistinguishable from a difference in the transform fitted to it.
# Deriving one set of parameters per PANEL, from a pooled subsample, keeps every
# sample on the same ruler. It is the same reason derive_cofactor_pooled()
# exists for the arcsinh path.


#' Estimate logicle parameters for one marker
#'
#' WHAT: the standard automatic-logicle rule (Parks, Roederer & Moore 2006;
#' Moore & Parks 2012). `t` is the top of scale, `m` the number of decades, and
#' `w` the width of the linear region, chosen so that the linear region just
#' covers the spread of the NEGATIVE population:
#'
#'     w = (m - log10(t / |r|)) / 2
#'
#' where `r` is a low quantile of the values below zero.
#'
#' WHY A QUANTILE OF THE NEGATIVES AND NOT THE MINIMUM: the minimum is one event.
#' On a few hundred thousand events it is whatever the noisiest cell in the tube
#' did, and `w` derived from it stretches the linear region until the positive
#' decades are squashed. The 5th percentile of the negative values describes the
#' bulk of that population instead of its worst member.
#'
#' Returns `w = 0` when a channel has no negative values at all -- with nothing
#' below zero there is nothing for a linear region to rescue, and logicle
#' degenerates to a log transform, which is the right answer rather than a
#' failure.
#'
#' WHY T IS THE CHANNEL'S RANGE AND NOT ITS MAXIMUM. T is "the top of scale"
#' (Moore & Parks 2012): a property of the instrument, the same for every file
#' acquired on it. The largest event in a sample is a property of that sample,
#' so a scale topped by it moves whenever one bright cell is or is not sampled,
#' and with it the width `w` and every cut on the channel. `t_range` is the
#' channel's `$PnR`, linearised where the channel was log-amplified (see
#' [linearise_channels()]); the data maximum is used only where no trustworthy
#' range was declared.
#'
#' @param x numeric vector of raw intensities for one marker
#' @param m decades on the display scale
#' @param a additional negative decades to display
#' @param neg_q quantile of the negative values used as `r`
#' @param t_range top of the channel's declared range, or NA. Default NA.
#' @return list(w, t, m, a)
#' @keywords internal
logicle_params_one <- function(x, m = 4.5, a = 0, neg_q = 0.05, t_range = NA_real_) {
  x <- x[is.finite(x)]
  if (!length(x)) return(list(w = 0.5, t = 262144, m = m, a = a))
  t <- if (length(t_range) == 1L && is.finite(t_range) && t_range > 0) t_range
       else max(x, na.rm = TRUE)
  if (!is.finite(t) || t <= 0) t <- 262144
  neg <- x[x < 0]
  if (!length(neg)) return(list(w = 0, t = t, m = m, a = a))
  r <- stats::quantile(neg, neg_q, na.rm = TRUE, names = FALSE)
  if (!is.finite(r) || r >= 0) return(list(w = 0, t = t, m = m, a = a))
  w <- (m - log10(t / abs(r))) / 2
  # Clamp: w < 0 means the negatives are narrower than the display can resolve
  # (nothing to widen), w > m/2 would consume the whole scale with linear region
  # and leave no decades. Both are the documented bounds of the method.
  w <- max(0, min(w, m / 2))
  list(w = w, t = t, m = m, a = a)
}


#' Derive one set of transform parameters for a panel
#'
#' Pooled across up to `max_samples` files so that every sample in a panel is
#' transformed on the same scale. See the note at the top of this file for why
#' that matters more than fitting each file well.
#'
#' @param reads list of read_fcs_resolved() results
#' @param sids sample ids belonging to one panel
#' @param markers marker names to parameterise
#' @param m,a logicle display parameters
#' @param max_samples files pooled
#' @param max_events events drawn per file
#' @param seed RNG seed; the stream is restored on exit
#' @return named list of per-marker logicle parameters
#' @export
derive_logicle_pooled <- function(reads, sids, markers = NULL, m = 4.5, a = 0,
                                  max_samples = 8L, max_events = 50000L,
                                  seed = 42L) {
  sids <- intersect(sids, names(reads))
  if (!length(sids)) return(NULL)
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  pick <- function(x) if (length(x) > max_samples)
    x[round(seq(1, length(x), length.out = max_samples))] else x
  use <- pick(sids)
  # EVERY MARKER OF THE PANEL, POOLED FROM THE FILES THAT CARRY IT. The list used
  # to come from the first file alone, and marker_cols is a named integer vector,
  # where [[ on an absent name is an error rather than NULL. Under
  # --panel-optional-markers some files lack a marker, so a run with logicle
  # stopped with "subscript out of bounds", and a marker absent from the first
  # file got no parameters at all. For a marker every file carries the files
  # pooled are exactly those pooled before.
  has_mk <- function(s, mk) mk %in% names(reads[[s]]$marker_cols)
  if (is.null(markers))
    markers <- unique(unlist(lapply(sids, function(s) names(reads[[s]]$marker_cols))))

  out <- list()
  for (mk in markers) {
    carriers <- sids[vapply(sids, has_mk, logical(1), mk = mk)]
    from <- if (length(carriers) == length(sids)) use else pick(carriers)
    vals <- unlist(lapply(from, function(s) {
      rd <- reads[[s]]
      if (is.null(rd$exprs) || !has_mk(s, mk)) return(NULL)
      v <- rd$exprs[, rd$marker_cols[[mk]]]
      if (length(v) > max_events) v <- v[sort(sample(length(v), max_events))]
      v
    }), use.names = FALSE)
    if (!length(vals)) next
    # The widest declared range among the files pooled, so no file's data sit
    # above the top of the shared scale.
    tops <- vapply(from, function(s) {
      rd <- reads[[s]]
      if (!has_mk(s, mk) || is.null(rd$channel_range)) return(NA_real_)
      as.numeric(rd$channel_range[rd$marker_cols[[mk]]])
    }, numeric(1))
    t_rng <- if (all(is.finite(tops))) max(tops) else NA_real_
    out[[mk]] <- logicle_params_one(vals, m = m, a = a, t_range = t_rng)
  }
  if (!length(out)) return(NULL)
  attr(out, "n_samples_pooled") <- length(use)
  out
}


#' Build the intensity transform the whole run uses
#'
#' WHAT IT RETURNS: an object with `fn(x, marker)`, so every call site applies
#' the transform the same way and none of them needs to know which method is in
#' force. Adding a method here reaches the gating, scoring, embedding and
#' figure code without touching any of them.
#'
#' @param method "arcsinh", "logicle" or "none"
#' @param cofactor arcsinh cofactor (required for method = "arcsinh")
#' @param logicle named list of per-marker parameters from
#'   [derive_logicle_pooled()] (required for method = "logicle")
#' @param per_marker named numeric vector of arcsinh cofactors by marker, or
#'   NULL. A marker named here is transformed with its own cofactor; any other,
#'   and scatter, with `cofactor`. See [derive_cofactor()].
#' @return list(method, fn, params, label)
#' @export
make_transform <- function(method = c("arcsinh", "logicle", "none"),
                           cofactor = NULL, logicle = NULL, per_marker = NULL) {
  method <- match.arg(method)

  if (method == "arcsinh") {
    if (!is.numeric(cofactor) || !is.finite(cofactor) || cofactor <= 0)
      stop("arcsinh needs a positive cofactor, got: ",
           paste(utils::capture.output(str(cofactor)), collapse = " "),
           call. = FALSE)
    if (!is.null(per_marker)) {
      per_marker <- per_marker[is.finite(per_marker) & per_marker > 0]
      if (!length(per_marker)) per_marker <- NULL
    }
    cof <- function(marker)
      if (!is.null(marker) && !is.null(per_marker) && marker %in% names(per_marker))
        per_marker[[marker]] else cofactor
    return(list(method = "arcsinh",
                params = list(cofactor = cofactor, per_marker = per_marker),
                label = "asinh intensity",
                cofactor_of = cof,
                fn = function(x, marker = NULL) asinh(x / cof(marker)),
                inv = function(y, marker = NULL) sinh(y) * cof(marker)))
  }

  if (method == "none")
    return(list(method = "none", params = list(), label = "raw intensity",
                fn = function(x, marker = NULL) x,
                inv = function(y, marker = NULL) y))

  if (is.null(logicle) || !length(logicle))
    stop("method = \"logicle\" needs parameters from derive_logicle_pooled()",
         call. = FALSE)
  # One flowCore transform closure per marker, built once. Rebuilding them per
  # call would dominate the runtime of every gating step.
  fns <- lapply(logicle, function(p)
    flowCore::logicleTransform(w = p$w, t = p$t, m = p$m, a = p$a))
  # Inverses, built by monotone interpolation rather than by calling back into
  # flowCore. The logicle is strictly increasing, so a dense grid inverts it to
  # far better precision than the only consumer needs (polygon vertices on their
  # way into a gate file). Doing it this way keeps the inverse working whatever
  # flowCore names its inverse constructor in a given release, and costs one
  # 4001-point evaluation per marker at construction rather than per call.
  invs <- lapply(names(logicle), function(mk) {
    p <- logicle[[mk]]; f <- fns[[mk]]
    hi <- max(p$t, 1)
    pos <- exp(seq(log(1e-2), log(hi), length.out = 2000L))
    xg  <- unique(sort(c(-rev(pos), 0, pos)))
    yg  <- f(xg)
    ok  <- is.finite(xg) & is.finite(yg)
    xg  <- xg[ok]; yg <- yg[ok]
    function(y) stats::approx(yg, xg, xout = y, rule = 2L)$y
  })
  names(invs) <- names(logicle)
  list(method = "logicle", params = logicle, label = "logicle intensity",
       fn = function(x, marker = NULL) {
         f <- if (!is.null(marker)) fns[[marker]] else NULL
         # A marker with no fitted parameters is a caller error, not something
         # to paper over with a default scale: a silently different transform on
         # one channel is exactly the failure this file exists to prevent.
         if (is.null(f))
           stop("no logicle parameters for marker '",
                marker %||% "<unnamed>", "'", call. = FALSE)
         f(x)
       },
       inv = function(y, marker = NULL) {
         g <- if (!is.null(marker)) invs[[marker]] else NULL
         if (is.null(g))
           stop("no logicle parameters for marker '",
                marker %||% "<unnamed>", "'", call. = FALSE)
         g(y)
       })
}
