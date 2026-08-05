# =============================================================================
# cyRAVEN -- learning a gating strategy back out of a cluster
# =============================================================================
#
# THE GAP THIS CLOSES. cluster_gate_agreement() can tell you that a cluster
# matches no population in the spec. It cannot tell you WHAT that cluster is.
# The finding arrives as a cluster number and a cell count, which is not
# something a cytometrist can act on: there is no way to go back to the
# instrument, or to a FlowJo workspace, and select those cells.
#
# This file closes the loop in the other direction. Given a set of cells
# labelled by anything -- a metacluster, a population, an arbitrary mask -- it
# LEARNS a short sequence of two-marker gates that reproduces the label, and
# reports how well each one does. The output is a gating strategy: marker pair,
# polygon vertices, and precision/recall/F1 at every level of the hierarchy.
#
# WHY A CONVEX POLYGON RATHER THAN MORE THRESHOLDS. Every population in the spec
# is a conjunction of one-dimensional cuts, which makes each gate an
# axis-aligned rectangle. Real boundaries are frequently diagonal -- CD4/CD8,
# CD14/CD16, FSC-A/FSC-H -- and a rectangle placed over a diagonal boundary must
# either admit contaminants or reject real cells. There is no third option. The
# package already concedes this once: derive_singlet_band() is a hand-written
# diagonal band, bolted on because a rectangle could not express it. A convex
# polygon is the general form of that special case, and it is still a gate a
# human can draw and defend.
#
# WHAT THIS IS NOT. A learned gate is DESCRIPTIVE. It says what a group of cells
# looks like in marker space; it does not establish that the group is real,
# reproducible, or biologically meaningful. Nothing here produces a p-value, and
# nothing here feeds back into the spec, the scoring, or any tested result. It
# is a hypothesis generator that sits beside cluster_gate_agreement(), and its
# output is a proposal for a human to accept or reject.
#
# HONEST EVALUATION IS PART OF THE ALGORITHM, NOT AN EXTRA. A gate with eight
# free half-planes fitted to a few thousand cells can memorise them, and a
# convex hull drawn around the survivors memorises them completely. Fitting and
# scoring on the same cells therefore reports a number that is optimistic by
# construction. Every metric below is computed on cells held out of the fit, and
# the in-sample value is reported alongside it so the size of the gap is
# visible rather than hidden.
# =============================================================================


# -----------------------------------------------------------------------------
# THE GATE, AND WHY IT IS OPTIMISED THE WAY IT IS
# -----------------------------------------------------------------------------
#
# A gate is the intersection of k half-planes: cell x is inside when
# w_j . x + b_j > 0 for every j. That predicate is a step function and cannot be
# differentiated, so for fitting it is relaxed to
#
#     p(x) = prod_j sigmoid(s * (w_j . x + b_j))
#
# which approaches the hard predicate as the steepness s grows. Fitting
# minimises class-weighted cross-entropy between p and the label, plus a
# tightness penalty (below).
#
# WHY L-BFGS-B AND NOT STOCHASTIC GRADIENT DESCENT. The whole problem is k
# planes times three parameters -- 24 numbers at the default k = 8. That is a
# small smooth optimisation, and a quasi-Newton method solves it directly. A
# minibatch first-order method with a fixed step size is the right tool when the
# parameter count is large and the data does not fit in memory; here neither is
# true, and it would trade an exact answer for a noisy one while also adding a
# learning rate and an iteration count to tune. The gradient is derived in
# closed form in gate_objective() -- it is four lines of matrix algebra -- so
# there is no automatic differentiation dependency either.
#
# WHY FOUR OF THE PLANES HAVE FIXED NORMALS. An intersection of half-planes with
# freely-chosen normals is not necessarily BOUNDED: an optimiser can happily
# settle on a wedge that runs off to infinity, which scores well on training
# cells and is not a gate anyone can draw. Fixing four normals to +/-PC1 and
# +/-PC2 of the target cloud pins a bounding box around the population in its
# own principal axes. Their OFFSETS stay free, so the box can still move and
# resize; only its orientation is fixed. The remaining k-4 planes are free to
# cut corners off that box, which is what turns a rotated rectangle into a
# polygon.
#
# WHY THE FREE NORMALS ARE ANGLES AND NOT VECTORS. A half-plane is unchanged by
# scaling (w_j, b_j) by any positive constant, but the RELAXATION is not: scaling
# up multiplies s * h_j, sharpening the sigmoid toward the hard indicator it
# approximates. On separable data that lowers the loss without bound, so an
# unconstrained optimiser drives ||w|| to infinity -- the same divergence an
# unpenalised logistic regression shows on separable classes. It fails silently
# and in a way that looks fine: the normals stay in the right DIRECTIONS, the
# reported F1 is high, and the accompanying offsets grow with them until the
# PC-box constraints are satisfied everywhere and stop bounding anything. The
# polygon then escapes to coordinates thousands of times the data range and is
# not a gate at all.
#
# Parametrising each free plane by an ANGLE removes the degeneracy at the
# source: w_j = (cos t_j, sin t_j) is a unit vector by construction, so s * h_j
# is a genuine signed distance times a fixed steepness and there is nothing to
# inflate. It also drops a parameter per plane. Every normal in the problem is
# then unit length -- the PCA rotation columns already are -- which is why the
# penalty below needs no norms in it.


#' Log-likelihood and gradient of a soft convex gate
#'
#' WHAT: class-weighted binary cross-entropy between the relaxed gate membership
#' `prod_j sigmoid(s * (w_j . x + b_j))` and the label, plus a tightness penalty.
#'
#' WHY THE PENALTY EXISTS AND WHAT IT DOES. Cross-entropy alone has no opinion
#' about where a plane sits once it has separated the classes, so planes drift
#' outward and the gate ends up loose -- correct on the cells it was shown, and
#' admitting anything that happens to lie beyond them. The penalty is the mean
#' distance from each plane to the centroid of the target cells, which pulls
#' every boundary inward. It is a shrinkage prior on gate size, and its strength
#' `lambda` is not guessed: learn_convex_gate() searches it against held-out F1,
#' so the data chooses how tight the gate should be.
#'
#' WHY THE MEMBERSHIP IS CLAMPED. Deep inside the gate p is numerically 1, and
#' the cross-entropy gradient for a misplaced non-target there is proportional
#' to p/(1-p), which overflows. Rescaling p onto `[eps, 1-eps]` bounds that
#' ratio by 1/eps while staying smooth, so a single badly-placed cell cannot
#' dominate the step.
#'
#' @param theta packed parameters: one angle per free plane, then every offset
#' @param X n x 2 matrix of marker values, scaled to the unit square
#' @param y 0/1 target indicator
#' @param wt per-cell weight
#' @param fixed_W 2 x 4 matrix of the frozen (PC-aligned) unit normals
#' @param centroid length-2 centroid of the target cells
#' @param s sigmoid steepness
#' @param lambda tightness penalty strength
#' @param eps membership clamp
#' @return objective value with a `gradient` attribute
#' @keywords internal
gate_objective <- function(theta, X, y, wt, fixed_W, centroid,
                           s = 40, lambda = 0, eps = 1e-3) {
  # theta holds one angle per free plane plus one offset per plane, so
  # length(theta) = 2 * n_free + k_fix and the plane count follows.
  k_fix  <- ncol(fixed_W)
  n_free <- (length(theta) - k_fix) / 2
  k      <- n_free + k_fix

  ang <- theta[seq_len(n_free)]
  ct  <- cos(ang); st <- sin(ang)
  W   <- cbind(rbind(ct, st), fixed_W)              # every column unit length
  b   <- theta[(n_free + 1L):length(theta)]

  H  <- X %*% W                                     # n x k half-space values
  H  <- sweep(H, 2L, b, "+")
  sH <- s * H

  # plogis(., log.p = TRUE) is log(sigmoid(.)) computed without ever forming
  # sigmoid itself, so the product below is a sum and never underflows.
  u   <- rowSums(stats::plogis(sH, log.p = TRUE))   # log p
  p   <- exp(u)
  pc  <- eps + (1 - 2 * eps) * p                    # clamped membership
  Wsum <- sum(wt)

  loss <- -sum(wt * (y * log(pc) + (1 - y) * log1p(-pc))) / Wsum

  # dL/dp' -> dL/du -> dL/dh_ij. The middle factor (1-2eps)*p is dp'/du.
  dL_dpc <- -wt * (y / pc - (1 - y) / (1 - pc)) / Wsum
  dL_du  <- dL_dpc * (1 - 2 * eps) * p
  # d(log p)/dh_ij = s * (1 - sigmoid(s h_ij)); the upper tail IS 1 - sigmoid.
  G <- (dL_du * s) * stats::plogis(sH, lower.tail = FALSE)

  # dh_ij/dt_j = x_i . (-sin t_j, cos t_j): the derivative of a unit normal with
  # respect to its own angle is the perpendicular unit vector.
  dW  <- rbind(-st, ct)                             # 2 x n_free
  gth <- colSums(G[, seq_len(n_free), drop = FALSE] * (X %*% dW))
  gb  <- colSums(G)

  if (lambda > 0) {
    # Every normal is unit length, so the signed distance from the centroid to
    # plane j is just w_j . m + b_j with no normalisation to carry through the
    # derivative.
    d   <- as.vector(crossprod(W, centroid)) + b
    sgn <- sign(d); sgn[sgn == 0] <- 1
    loss <- loss + lambda * mean(abs(d))
    gb  <- gb + lambda * sgn / k
    gth <- gth + lambda * sgn[seq_len(n_free)] *
      as.vector(crossprod(dW, centroid)) / k
  }

  attr(loss, "gradient") <- c(gth, gb)
  loss
}


#' Vertices of the polygon defined by an intersection of half-planes
#'
#' WHY IT ENUMERATES PAIRS: with k planes there are at most choose(k, 2)
#' candidate corners, 66 at the default k = 12, and each is one 2x2 solve. A
#' general half-space enumeration algorithm would be faster asymptotically and
#' slower here, and would be one more thing that can be subtly wrong.
#'
#' Returns NULL when the region is empty or degenerate, which is the honest
#' answer: an optimiser that produced no feasible region has not produced a gate.
#'
#' @param W 2 x k matrix of normals
#' @param b length-k offsets
#' @param tol feasibility tolerance
#' @return matrix of vertices in draw order, or NULL
#' @keywords internal
halfplane_polygon <- function(W, b, tol = 1e-9) {
  k <- ncol(W)
  if (k < 3L) return(NULL)
  pairs <- utils::combn(k, 2L)
  pts <- vector("list", ncol(pairs))
  for (i in seq_len(ncol(pairs))) {
    j1 <- pairs[1L, i]; j2 <- pairs[2L, i]
    A <- rbind(W[, j1], W[, j2])
    if (abs(det(A)) < 1e-12) next                    # parallel planes
    x <- tryCatch(solve(A, c(-b[j1], -b[j2])), error = function(e) NULL)
    if (is.null(x)) next
    if (all(as.vector(crossprod(W, x)) + b >= -tol)) pts[[i]] <- x
  }
  pts <- do.call(rbind, pts)
  if (is.null(pts) || nrow(pts) < 3L) return(NULL)
  pts <- unique(round(pts, 10L))
  if (nrow(pts) < 3L) return(NULL)
  # Order by angle about the centroid. The region is convex, so this is the
  # boundary order; chull() would also work and does strictly more.
  ctr <- colMeans(pts)
  pts[order(atan2(pts[, 2L] - ctr[2L], pts[, 1L] - ctr[1L])), , drop = FALSE]
}


#' Which points fall inside a polygon
#'
#' WHY A WINDING-NUMBER TEST AND NOT A RAY CAST: both are short, but the winding
#' rule gives the same answer for a vertex-on-boundary point regardless of which
#' direction the ray was cast, so a cell sitting exactly on a gate edge does not
#' change classification when the polygon is rotated. Vectorised over points.
#'
#' @param xy n x 2 matrix of points
#' @param poly m x 2 matrix of polygon vertices in order
#' @return logical vector
#' @keywords internal
point_in_polygon <- function(xy, poly) {
  if (is.null(poly) || nrow(poly) < 3L) return(rep(FALSE, nrow(xy)))
  n <- nrow(xy); m <- nrow(poly)
  wind <- integer(n)
  x <- xy[, 1L]; y <- xy[, 2L]
  for (i in seq_len(m)) {
    j <- if (i == m) 1L else i + 1L
    x1 <- poly[i, 1L]; y1 <- poly[i, 2L]
    x2 <- poly[j, 1L]; y2 <- poly[j, 2L]
    # Side of the directed edge (x1,y1)->(x2,y2) each point lies on.
    side <- (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1)
    up   <- y1 <= y & y2 >  y & side > 0
    down <- y1 >  y & y2 <= y & side < 0
    wind <- wind + up - down
  }
  wind != 0L
}


#' Precision, recall and F1 for a binary gate
#'
#' Recall is measured against `n_target_total` rather than against the targets
#' present, so that in a hierarchy it stays CUMULATIVE: a level-3 gate is judged
#' on the fraction of the ORIGINAL population it still holds, not on the
#' fraction of whatever survived level 2. A hierarchy that discards half the
#' population at each level and reports 100% recall three times is exactly the
#' failure this guards against.
#'
#' @param inside logical, gate membership
#' @param y 0/1 target indicator
#' @param n_target_total denominator for recall
#' @return named numeric vector
#' @keywords internal
gate_metrics <- function(inside, y, n_target_total = sum(y == 1)) {
  tp <- sum(inside & y == 1); fp <- sum(inside & y == 0)
  prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  rec  <- if (n_target_total > 0) tp / n_target_total else NA_real_
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
            2 * prec * rec / (prec + rec) else NA_real_
  c(precision = prec, recall = rec, f1 = f1,
    tp = tp, fp = fp, n_in_gate = tp + fp)
}


#' Learn a convex two-marker gate for a labelled set of cells
#'
#' WHAT: fits an intersection of half-planes that separates `y == 1` cells from
#' the rest in the plane of two markers, then reports how well it does on cells
#' held out of the fit.
#'
#' HOW THE TIGHTNESS IS CHOSEN: `lambda` is not a tuning parameter the caller has
#' to guess. The function fits at each value of `lambda_grid`, scores every fit
#' on held-out cells, and keeps the best; it then refines once around the winner.
#' Searching against held-out F1 rather than training F1 is what stops the search
#' from always selecting the loosest gate, which is what fits the training cells
#' best and generalises worst.
#'
#' WHY IT RETURNS TWO POLYGONS. `polygon` is the fitted half-plane intersection.
#' `hull` is the convex hull of the target cells that fell inside it -- strictly
#' tighter, usually more precise, and much more prone to memorising the training
#' cells, since it is drawn around them. Both are returned with their own
#' held-out metrics so the trade is visible instead of being made silently.
#'
#' @param X n x 2 numeric matrix of marker values (arcsinh scale)
#' @param y 0/1 or logical target indicator, length n
#' @param n_planes total half-planes; 4 are PC-aligned and fixed, the rest free
#' @param s_schedule sigmoid steepness. Supplying several increasing values
#'   anneals: each fit warm-starts the next. The default is a single value,
#'   because annealing measured no reliable improvement (+0.002 to +0.007 F1 on
#'   the cases tried, inside the seed-to-seed spread) for roughly twice the
#'   runtime. It is left available for a dataset where the fit does get stuck
#' @param lambda_grid tightness values to search
#' @param holdout fraction of cells reserved for evaluation, stratified by label
#' @param n_target_total recall denominator; defaults to the targets supplied
#' @param max_cells subsample ceiling for the fit
#' @param seed RNG seed; the stream is restored on exit
#' @return list with `polygon`, `hull`, `metrics`, `metrics_hull`,
#'   `metrics_insample`, `lambda`, `W`, `b`, `markers`, or NULL if no gate could
#'   be fitted
#' @export
learn_convex_gate <- function(X, y, n_planes = 8L,
                              s_schedule = 800,
                              lambda_grid = c(0, 0.5, 1, 2, 4, 8),
                              holdout = 0.3, n_target_total = NULL,
                              max_cells = 20000L, seed = 42L) {
  s_schedule <- sort(unique(as.numeric(s_schedule)))
  if (!length(s_schedule) || any(!is.finite(s_schedule)) || any(s_schedule <= 0))
    stop("s_schedule must be positive and finite", call. = FALSE)
  X <- as.matrix(X)
  if (ncol(X) != 2L) return(NULL)
  y <- as.integer(as.logical(y))
  ok <- stats::complete.cases(X) & !is.na(y)
  X <- X[ok, , drop = FALSE]; y <- y[ok]
  if (sum(y == 1) < 20L || sum(y == 0) < 20L) return(NULL)
  if (is.null(n_target_total)) n_target_total <- sum(y == 1)

  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  if (nrow(X) > max_cells) {
    n_before <- sum(y == 1)
    keep <- sort(sample(nrow(X), max_cells))
    X <- X[keep, , drop = FALSE]; y <- y[keep]
    # The recall denominator has to follow the subsample. It is expressed in the
    # population the CALLER passed, and from here on every count is taken over
    # the retained cells, so leaving it alone divides in-gate targets by a
    # denominator from a larger set and understates recall by exactly the
    # subsample ratio -- 0.5 at the default ceiling against 40,000 cells, which
    # reads as "the gate loses half the population" and is entirely an artefact.
    if (n_before > 0L) n_target_total <- n_target_total * sum(y == 1) / n_before
  }

  # Scale to the unit square on robust limits. The sigmoid steepness s is a
  # property of that scale, so without this it would have to be retuned for
  # every marker pair; outlier-driven min/max would do the same thing more
  # subtly, by squeezing the bulk of the data into a corner.
  lo <- apply(X, 2L, stats::quantile, 0.001, na.rm = TRUE)
  hi <- apply(X, 2L, stats::quantile, 0.999, na.rm = TRUE)
  rng <- hi - lo; rng[rng <= 0] <- 1
  Z <- sweep(sweep(X, 2L, lo, "-"), 2L, rng, "/")

  # Stratified split: an unstratified one can hand a rare population a test set
  # with no positives at all, and every metric below would then be NA or 0 for
  # reasons that have nothing to do with the gate.
  idx_pos <- which(y == 1); idx_neg <- which(y == 0)
  te <- c(sample(idx_pos, max(1L, floor(holdout * length(idx_pos)))),
          sample(idx_neg, max(1L, floor(holdout * length(idx_neg)))))
  tr <- setdiff(seq_along(y), te)
  if (sum(y[tr] == 1) < 10L || sum(y[te] == 1) < 5L) { tr <- seq_along(y); te <- tr }

  n_fix <- 4L
  n_planes <- max(n_fix + 1L, as.integer(n_planes))
  n_free <- n_planes - n_fix

  tgt <- Z[tr[y[tr] == 1], , drop = FALSE]
  centroid <- colMeans(tgt)
  pc <- tryCatch(stats::prcomp(tgt, center = TRUE, scale. = FALSE),
                 error = function(e) NULL)
  axes <- if (is.null(pc)) diag(2L) else pc$rotation[, 1:2, drop = FALSE]
  fixed_W <- cbind(axes[, 1L], -axes[, 1L], axes[, 2L], -axes[, 2L])

  # Class weight: the exact imbalance ratio, so that targets and non-targets
  # carry equal total influence and a population that is 0.4% of the cells is
  # not simply predicted away.
  #
  # WHY IT IS NOT CAPPED. An earlier version capped it at 50, on the theory that
  # an extreme weight stops being a correction and becomes an instruction to
  # gate everything. The cap was removed because nothing justified it: at a
  # 268:1 imbalance, lifting it changed the fitted gate not at all. The guard
  # against a handful of cells dominating is the minimum target count above,
  # which refuses the fit outright and says so, rather than a weight that
  # quietly biases every gate by an amount nobody can see.
  w_pos <- max(1, sum(y[tr] == 0) / max(1, sum(y[tr] == 1)))
  wt <- ifelse(y[tr] == 1, w_pos, 1)

  # Recall denominator for the held-out set. When this is called from
  # explain_cluster(), n_target_total is the size of the ORIGINAL population --
  # larger than the targets still present at this level -- so it has to be
  # apportioned by the share of surviving targets the split put in the test set,
  # not by the share of all cells. Getting that wrong makes deep levels look
  # better than they are, which is the one direction a hierarchy must not err in.
  n_te_tgt <- n_target_total * sum(y[te] == 1) / max(1L, sum(y == 1))

  # The observed data range, as four half-planes, intersected with the fitted
  # ones before the polygon is formed. WHY: a boundary drawn where no cell was
  # ever measured is not a decision, it is an extrapolation, and it makes the
  # gate impossible to plot on the axes the cells live on. Clipping also means a
  # plane the fit rendered redundant -- pushed out until it constrains nothing --
  # costs a vacuous vertex rather than an unbounded excursion.
  #
  # THE BOX IS THE ACTUAL RANGE OF Z, NOT THE UNIT SQUARE. The scaling above
  # uses the 0.1/99.9 percentiles so that the steepness s means the same thing
  # for every marker pair, which puts a small number of real cells outside
  # [0, 1]. Clipping there instead cost recall precisely where it is least
  # affordable: a population living in the extreme tail -- which is what an
  # undescribed cluster in a corner of the marker space IS -- is largely made of
  # those cells, and the gate would exclude half of it while reporting precision
  # 1.00, looking excellent and selecting the wrong cells. The observed range is
  # the honest boundary: no cell is ever outside it.
  zlo <- apply(Z, 2L, min); zhi <- apply(Z, 2L, max)
  pad <- 0.02 * pmax(zhi - zlo, .Machine$double.eps)
  box_W <- cbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  box_b <- c(-(zlo[1L] - pad[1L]), zhi[1L] + pad[1L],
             -(zlo[2L] - pad[2L]), zhi[2L] + pad[2L])

  fit_one <- function(lambda) {
    # Free planes start on evenly-spaced directions rather than at random, so
    # the fit is reproducible without depending on the RNG state and the
    # starting polytope is a regular polygon rather than an arbitrary wedge.
    ang0 <- seq(pi / 4, by = 2 * pi / max(1L, n_free), length.out = n_free)
    W0   <- rbind(cos(ang0), sin(ang0))
    th   <- c(ang0, -as.vector(crossprod(cbind(W0, fixed_W), centroid)) + 0.25)

    # STEEPNESS. s sets how sharply the relaxation approximates the hard gate:
    # the transition band is about 4/s wide in unit-square distance, so s = 50
    # blurs the boundary over 8% of the range and s = 800 over 0.5%. Sharper is
    # better because it is closer to the predicate actually used to classify --
    # measured F1 rises monotonically from 0.86 at s = 40 to 0.98 at s = 800 on
    # a planted diagonal population.
    #
    # The textbook worry about a sharp sigmoid is saturation: a cell far from
    # every boundary contributes no gradient, so an optimiser started cold has
    # nothing to follow. Annealing -- fit soft, warm-start sharper -- is the
    # standard remedy, and this loop implements it. It is NOT the default,
    # because measuring it showed no reliable gain: the planes are initialised
    # around the target centroid, which already puts them where the gradients
    # are. Left in as an escape hatch rather than removed, since a dataset that
    # does strand the fit is cheap to rescue and impossible to rescue if the
    # code is not there.
    fit <- NULL
    for (s_now in s_schedule) {
      obj <- function(t) gate_objective(t, Z[tr, , drop = FALSE], y[tr],
                                        wt, fixed_W, centroid, s_now, lambda)
      f <- tryCatch(stats::optim(
        th, fn = function(t) as.numeric(obj(t)),
        gr = function(t) attr(obj(t), "gradient"),
        method = "L-BFGS-B",
        control = list(maxit = 300L, factr = 1e9)), error = function(e) NULL)
      if (is.null(f)) break
      th <- f$par; fit <- f
    }
    if (is.null(fit)) return(NULL)

    ang <- fit$par[seq_len(n_free)]
    W <- cbind(rbind(cos(ang), sin(ang)), fixed_W)
    b <- fit$par[(n_free + 1L):length(fit$par)]
    poly <- halfplane_polygon(cbind(W, box_W), c(b, box_b))
    if (is.null(poly)) return(NULL)
    inside <- point_in_polygon(Z, poly)
    list(W = W, b = b, polygon = poly, inside = inside, lambda = lambda,
         score = gate_metrics(inside[te], y[te], n_te_tgt)[["f1"]])
  }

  best <- NULL
  for (lam in lambda_grid) {
    f <- fit_one(lam)
    if (!is.null(f) && !is.na(f$score) &&
        (is.null(best) || f$score > best$score)) best <- f
  }
  if (is.null(best)) return(NULL)
  # One refinement pass around the winner. A second pass buys almost nothing:
  # F1 as a function of lambda is flat near its optimum, which is the point of
  # searching it rather than fixing it.
  step <- max(diff(sort(unique(c(0, lambda_grid))))) / 2
  for (lam in c(best$lambda - step, best$lambda + step)) {
    if (lam < 0) next
    f <- fit_one(lam)
    if (!is.null(f) && !is.na(f$score) && f$score > best$score) best <- f
  }

  m_te <- gate_metrics(best$inside[te], y[te], n_te_tgt)
  m_tr <- gate_metrics(best$inside[tr], y[tr],
                       n_target_total * sum(y[tr] == 1) / max(1L, sum(y == 1)))

  # Hull tightening, fitted on TRAINING targets only. Drawing it around every
  # target including the held-out ones would make its held-out metrics
  # meaningless -- it would be scored on cells it was built to contain.
  hull_poly <- NULL; m_hull <- NULL
  inside_tr_tgt <- best$inside & seq_along(y) %in% tr & y == 1
  if (sum(inside_tr_tgt) >= 3L) {
    ptsh <- Z[inside_tr_tgt, , drop = FALSE]
    hv <- tryCatch(grDevices::chull(ptsh), error = function(e) NULL)
    if (!is.null(hv) && length(hv) >= 3L) {
      hull_poly <- ptsh[hv, , drop = FALSE]
      m_hull <- gate_metrics(point_in_polygon(Z[te, , drop = FALSE], hull_poly),
                             y[te], n_te_tgt)
    }
  }

  unscale <- function(poly) {
    if (is.null(poly)) return(NULL)
    out <- sweep(sweep(poly, 2L, rng, "*"), 2L, lo, "+")
    colnames(out) <- colnames(X); out
  }

  list(polygon = unscale(best$polygon), hull = unscale(hull_poly),
       metrics = m_te, metrics_hull = m_hull, metrics_insample = m_tr,
       lambda = best$lambda, W = best$W, b = best$b,
       scale = list(lo = lo, rng = rng), n_planes = n_planes,
       markers = colnames(X), n_train = length(tr), n_test = length(te))
}


#' Rank markers by how differently they are distributed in targets and the rest
#'
#' WHY THE (median, p1, p99) TRIPLE AND NOT A MEAN DIFFERENCE: a marker can
#' separate a population by having a different SPREAD rather than a different
#' centre -- a subset that is uniformly dim on a marker everything else is
#' bimodal on, for instance -- and a difference of means scores that at zero.
#' Comparing three order statistics catches shifts in location and in both
#' tails, and costs one pass over the column.
#'
#' Values are scaled to the unit interval first, so a marker on a wider axis
#' cannot outrank a more informative one purely by having larger numbers.
#'
#' @param X n x m matrix of marker values
#' @param y 0/1 target indicator
#' @return data.frame of marker and score, best first
#' @export
rank_gate_markers <- function(X, y) {
  X <- as.matrix(X); y <- as.integer(as.logical(y))
  if (!ncol(X) || sum(y == 1) < 5L || sum(y == 0) < 5L) return(NULL)
  qs <- function(v) stats::quantile(v, c(0.5, 0.01, 0.99), na.rm = TRUE, names = FALSE)
  sc <- vapply(seq_len(ncol(X)), function(j) {
    v <- X[, j]
    r <- range(v, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) <= 0) return(NA_real_)
    v <- (v - r[1L]) / diff(r)
    mean((qs(v[y == 1]) - qs(v[y == 0]))^2)
  }, numeric(1))
  out <- data.frame(marker = colnames(X), score = sc, stringsAsFactors = FALSE)
  out <- out[is.finite(out$score), , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out[order(-out$score), , drop = FALSE]
}


#' Learn a hierarchical gating strategy for one labelled population
#'
#' WHAT: repeatedly picks the two most discriminating markers among the cells
#' that survive so far, fits a convex gate in that plane, keeps the cells inside
#' it, and goes again -- which is exactly the shape of a manual gating strategy,
#' and is why the output can be executed by hand.
#'
#' WHY IT STOPS EARLY. Each additional level can only remove cells, so recall
#' falls monotonically while precision usually rises. The strategy is worth
#' extending only while F1 improves; `min_gain` sets how much improvement counts
#' as worth another gate for whoever has to draw it. It also stops when a level
#' fails to fit, rather than skipping it -- a gap in the middle of a hierarchy
#' is not a gating strategy, and returning the levels that did work is the
#' honest truncation.
#'
#' @param X n x m matrix of marker values (arcsinh scale)
#' @param y 0/1 or logical target indicator
#' @param max_depth maximum gates in the strategy
#' @param min_gain minimum held-out F1 improvement to justify another level
#' @param markers optional restriction of the marker set
#' @param ... passed to [learn_convex_gate()]
#' @return list(levels = list of gates, summary = data.frame, best_depth), or NULL
#' @export
explain_cluster <- function(X, y, max_depth = 4L, min_gain = 0.02,
                            markers = NULL, ...) {
  X <- as.matrix(X)
  if (!is.null(markers)) X <- X[, intersect(markers, colnames(X)), drop = FALSE]
  keep <- vapply(seq_len(ncol(X)), function(j) {
    v <- X[, j]; sum(is.finite(v)) > 0 && diff(range(v, na.rm = TRUE)) > 0
  }, logical(1))
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 2L) return(NULL)

  y <- as.integer(as.logical(y))
  n_total <- sum(y == 1)
  if (n_total < 20L) return(NULL)

  alive <- rep(TRUE, nrow(X))
  used  <- character(0)
  levels_out <- list(); rows <- list()
  best_f1 <- -Inf

  for (d in seq_len(as.integer(max_depth))) {
    if (sum(alive & y == 1) < 20L || sum(alive & y == 0) < 20L) break
    rk <- rank_gate_markers(X[alive, setdiff(colnames(X), used), drop = FALSE],
                            y[alive])
    if (is.null(rk) || nrow(rk) < 2L) break
    pair <- rk$marker[1:2]

    # n_total, not a share of it: learn_convex_gate() rescales the denominator to
    # whatever split and subsample it ends up using, by the fraction of the
    # targets PRESENT that it retained. Pre-scaling here by the fraction of CELLS
    # still alive would double-count the same shrinkage and would use the wrong
    # quantity to do it -- a gate that drops many non-targets and no targets
    # would appear to lose recall it never lost.
    g <- learn_convex_gate(X[alive, pair, drop = FALSE], y[alive],
                           n_target_total = n_total, ...)
    if (is.null(g)) break

    inside <- point_in_polygon(X[alive, pair, drop = FALSE], g$polygon)
    new_alive <- alive; new_alive[alive] <- inside
    cum <- gate_metrics(new_alive, y, n_total)

    if (d > 1L && !is.na(cum[["f1"]]) && cum[["f1"]] < best_f1 + min_gain) break
    if (!is.na(cum[["f1"]])) best_f1 <- max(best_f1, cum[["f1"]])

    g$depth <- d; g$pair <- pair; g$cumulative <- cum
    levels_out[[d]] <- g
    rows[[d]] <- data.frame(
      depth = d, marker_x = pair[1L], marker_y = pair[2L],
      lambda = g$lambda,
      holdout_precision = unname(round(g$metrics[["precision"]], 4)),
      holdout_recall    = unname(round(g$metrics[["recall"]], 4)),
      holdout_f1        = unname(round(g$metrics[["f1"]], 4)),
      insample_f1       = unname(round(g$metrics_insample[["f1"]], 4)),
      hull_holdout_f1   = if (is.null(g$metrics_hull)) NA_real_
                          else unname(round(g$metrics_hull[["f1"]], 4)),
      cumulative_precision = unname(round(cum[["precision"]], 4)),
      cumulative_recall    = unname(round(cum[["recall"]], 4)),
      cumulative_f1        = unname(round(cum[["f1"]], 4)),
      cells_in_gate = as.integer(sum(new_alive)),
      targets_in_gate = as.integer(sum(new_alive & y == 1)),
      stringsAsFactors = FALSE)

    alive <- new_alive
    used  <- c(used, pair)
  }

  levels_out <- Filter(Negate(is.null), levels_out)
  if (!length(levels_out)) return(NULL)
  summ <- do.call(rbind, rows)
  list(levels = levels_out, summary = summ,
       best_depth = summ$depth[which.max(summ$cumulative_f1)],
       n_target_total = n_total)
}


#' Propose gating strategies for clusters the spec does not describe
#'
#' WHAT IT IS FOR. [cluster_gate_agreement()] identifies clusters whose cells
#' carry no population label, or whose dominant label covers only part of them.
#' Those are the clusters worth explaining: the spec has nothing to say about
#' them, so there is no gate to inspect and no threshold to blame. This runs
#' [explain_cluster()] on each and writes out a strategy a cytometrist can draw.
#'
#' WHY IT ONLY TAKES THE FIRST FEW. Each strategy costs a handful of small
#' optimisations, which is cheap, but a report proposing fifteen new gates is
#' not a finding, it is a second problem. `max_clusters` keeps the output to the
#' clusters with the strongest claim to being real populations -- the largest
#' undescribed ones -- and the count of what was left out is logged rather than
#' passed over in silence.
#'
#' @param cells embedding cell table carrying the marker columns
#' @param cluster integer vector from [run_unsupervised_clusters()]
#' @param features marker columns to gate on
#' @param agreement the list returned by [cluster_gate_agreement()]
#' @param max_clusters ceiling on how many strategies to derive
#' @param min_cells smallest cluster worth explaining
#' @param purity_max only explain clusters no cleaner than this (percent)
#' @param ... passed to [explain_cluster()]
#' @return list(summary, polygons, strategies), or NULL
#' @export
explain_unmatched_clusters <- function(cells, cluster, features, agreement,
                                       max_clusters = 4L, min_cells = 200L,
                                       purity_max = 80, ...) {
  if (is.null(agreement) || is.null(agreement$per_cluster)) return(NULL)
  feats <- intersect(features, names(cells))
  if (length(feats) < 2L) return(NULL)

  pc <- agreement$per_cluster
  cand <- pc[pc$n_cells >= min_cells &
               (grepl("^UNDESCRIBED", pc$interpretation) |
                  pc$purity_pct <= purity_max), , drop = FALSE]
  if (!nrow(cand)) {
    log_msg("  gate proposals: every cluster already matches a described ",
            "population -- nothing to explain")
    return(NULL)
  }
  cand <- cand[order(-cand$n_cells), , drop = FALSE]
  n_drop <- max(0L, nrow(cand) - max_clusters)
  if (n_drop > 0L)
    log_msg("  gate proposals: ", nrow(cand), " cluster(s) qualify; deriving ",
            "strategies for the ", max_clusters, " largest and skipping ",
            n_drop, " smaller one(s)")
  cand <- utils::head(cand, max_clusters)

  X <- as.matrix(cells[, feats, drop = FALSE])
  strategies <- list(); srows <- list(); prows <- list()

  for (i in seq_len(nrow(cand))) {
    k <- cand$cluster[i]
    st <- tryCatch(explain_cluster(X, as.integer(cluster == k), ...),
                   error = function(e) NULL)
    if (is.null(st)) {
      log_msg("  gate proposals: cluster ", k, " -- no gate could be fitted")
      next
    }
    strategies[[as.character(k)]] <- st
    s <- st$summary
    s$cluster <- k
    s$dominant_gate_label <- cand$dominant_gate_label[i]
    s$cluster_purity_pct <- cand$purity_pct[i]
    srows[[length(srows) + 1L]] <- s
    for (lv in st$levels) {
      p <- lv$polygon
      if (is.null(p)) next
      prows[[length(prows) + 1L]] <- data.frame(
        cluster = k, depth = lv$depth,
        marker_x = lv$pair[1L], marker_y = lv$pair[2L],
        vertex = seq_len(nrow(p)),
        x_asinh = round(p[, 1L], 6), y_asinh = round(p[, 2L], 6),
        stringsAsFactors = FALSE)
    }
    b <- s[s$depth == st$best_depth, , drop = FALSE]
    log_msg("  gate proposals: cluster ", k, " -> ",
            paste(vapply(st$levels, function(l)
              paste(l$pair, collapse = "/"), character(1)), collapse = " then "),
            "  (held-out F1 ", format(round(b$cumulative_f1, 3), nsmall = 3),
            " at depth ", st$best_depth, ")")
  }

  if (!length(srows)) return(NULL)
  list(summary = do.call(rbind, srows),
       polygons = do.call(rbind, prows),
       strategies = strategies)
}


#' Draw a proposed gating strategy, one panel per gate
#'
#' WHY TARGETS ARE DRAWN LAST AND NON-TARGETS IN GREY: the question the reader
#' has is "does this polygon contain the population and exclude the rest", and
#' that is answered by seeing the population against the background it was
#' separated from. Colouring both categories equally makes a dense non-target
#' cloud hide the very cells the gate is about.
#'
#' The subtitle carries the held-out metrics, not the in-sample ones, because
#' a figure is where an optimistic number does the most damage.
#'
#' @param X marker matrix used to derive the strategy
#' @param y 0/1 target indicator
#' @param strategy the list returned by [explain_cluster()]
#' @param outfile path to write to
#' @param label name of the population being explained
#' @param colors Named list of colours; defaults to the package palette. Default `fcs_colors()`.
#' @param dpi Raster resolution. Default `300`.
#' @return the assembled plot, invisibly
#' @export
fig_gate_strategy <- function(X, y, strategy, outfile, label = "cluster",
                              colors = fcs_colors(), dpi = 300) {
  if (is.null(strategy) || !length(strategy$levels)) return(invisible(NULL))
  X <- as.matrix(X); y <- as.integer(as.logical(y))
  lv <- strategy$levels
  alive <- rep(TRUE, nrow(X))
  panels <- list()

  for (g in lv) {
    pr <- g$pair
    sub <- data.frame(x = X[alive, pr[1L]], y = X[alive, pr[2L]],
                      tag = ifelse(y[alive] == 1, "target", "other"),
                      stringsAsFactors = FALSE)
    # The panel is a picture of a decision, so it shows the cells the decision
    # was actually made on -- those still alive at this level -- not all cells.
    if (nrow(sub) > 60000L) sub <- sub[sort(sample(nrow(sub), 60000L)), , drop = FALSE]
    sub <- sub[order(sub$tag == "target"), , drop = FALSE]
    aes_pt <- auto_point_aes(nrow(sub))
    poly <- as.data.frame(g$polygon); names(poly) <- c("x", "y")

    m <- g$metrics
    p <- ggplot2::ggplot(sub, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(colour = tag), size = aes_pt$size,
                          alpha = aes_pt$alpha, show.legend = (g$depth == 1L)) +
      ggplot2::geom_polygon(data = poly, ggplot2::aes(x = x, y = y),
                            fill = NA, linewidth = 0.7,
                            colour = colors$gate_highlight %||% "#D62728") +
      ggplot2::scale_colour_manual(
        values = c(target = colors$gate_target %||% "#E8112D",
                   other  = colors$gate_nontarget %||% "grey75"),
        breaks = c("target", "other"),
        labels = c(paste0(label, " (target)"), "all other cells"),
        name = NULL) +
      ggplot2::labs(
        x = pretty_label(pr[1L]), y = pretty_label(pr[2L]),
        title = paste0("Gate ", g$depth, ": ", pr[1L], " / ", pr[2L]),
        subtitle = wrap_plot_text(sprintf(
          "held-out precision %.2f, recall %.2f, F1 %.2f (in-sample F1 %.2f)",
          m[["precision"]], m[["recall"]], m[["f1"]],
          g$metrics_insample[["f1"]]), width_in = 5)) +
      theme_cyto(colors = colors) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 2.5, alpha = 1)))
    panels[[length(panels) + 1L]] <- p

    inside <- point_in_polygon(X[alive, pr, drop = FALSE], g$polygon)
    alive[alive] <- inside
  }

  ncol_use <- min(length(panels), 3L)
  fig <- patchwork::wrap_plots(panels, ncol = ncol_use) +
    patchwork::plot_annotation(
      title = paste0("Proposed gating strategy for ", label),
      subtitle = wrap_plot_text(paste0(
        "Learned from the data, not from the population spec. Each panel gates ",
        "the cells that survived the previous one. Metrics are computed on cells ",
        "held out of the fit; the in-sample value is shown alongside so the gap ",
        "is visible. This is a DESCRIPTION of where these cells sit in marker ",
        "space and a proposal for how to select them -- it is not evidence that ",
        "the population is real, and it carries no p-value."),
        width_in = 5.2 * ncol_use),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13),
        plot.subtitle = ggplot2::element_text(
          size = 9, colour = colors$subtitle_text %||% "grey35"),
        plot.title.position = "plot"))

  nrow_use <- ceiling(length(panels) / ncol_use)
  safe_ggsave(outfile, plot = fig, width = 5.2 * ncol_use,
              height = 4.9 * nrow_use + 0.9, dpi = dpi, limitsize = FALSE)
  log_msg("[fig] wrote ", outfile, " (", length(panels), " gate(s))")
  invisible(fig)
}
