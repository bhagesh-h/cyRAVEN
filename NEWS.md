# cyRAVEN 0.1.0

First packaged release. Previously a set of command-line R scripts; the analysis
is unchanged, the delivery is not.

## Standalone coverage

* New `--transform logicle` alongside the existing arcsinh path. Logicle is the
  flow-cytometry convention: linear near zero, so compensated negative values
  stay on scale rather than piling against an axis limit. Parameters are pooled
  per PANEL, not fitted per file, so every sample stays on one ruler and
  cross-sample medians remain comparable.
* The transform choice is not cosmetic. On a seven-colour T-cell panel the
  top-level populations agree to 0.1% between methods while the CCR7/CD45RA
  memory subsets move by tens of percent, because under arcsinh the CCR7
  threshold resolves to a density valley in some samples and a quantile fallback
  in others. Logicle finds a valley in every sample. `thresholds_used.csv`
  records which happened.
* New `--correct-batch`, which aligns each marker across batches by monotone
  quantile mapping. Being monotone means it cannot reorder cells within a batch,
  so it corrects location and spread without inventing structure.
* Batch correction REFUSES when Cramer's V between batch and group exceeds
  `--batch-max-cramers-v`. Where the two overlap that far, removing the batch
  and removing the finding are the same operation.
  `--force-batch-correction` overrides and is recorded in the run manifest.

## Learned gating strategies

* New: `--explain-clusters` derives a gating strategy for any cluster the
  population specification does not describe. `cluster_gate_agreement()` could
  already report that a cluster matched nothing; it could not say what the
  cluster *was*, which left the finding as a cluster number and a cell count --
  nothing anyone could act on at the instrument. `explain_cluster()` learns a
  short sequence of two-marker gates that selects those cells and reports how
  well each performs.
* Gates are convex polygons -- intersections of half-planes -- rather than the
  axis-aligned rectangles a conjunction of one-dimensional thresholds can
  express. Diagonal boundaries (CD4/CD8, CD14/CD16, FSC-A/FSC-H) cannot be
  captured by a rectangle without either admitting contaminants or rejecting
  real cells; `derive_singlet_band()` was already a hand-written special case of
  exactly this, and a polygon is its general form.
* No new dependencies. The objective's gradient is derived in closed form and
  the fit is 24 parameters, so `stats::optim`'s L-BFGS-B solves it directly with
  no automatic-differentiation package; `stats::prcomp` initialises and
  `grDevices::chull` tightens, both already imported.
* Every reported metric is computed on cells held out of the fit, with the
  in-sample value printed beside it. Eight free half-planes fitted to a few
  thousand cells can memorise them, so a resubstitution F1 is optimistic by
  construction and the size of that gap is information the reader needs.
* Strictly descriptive. A proposed gate says where cells sit in marker space,
  not that they are a real population; it carries no p-value, and the scored
  populations, frequencies and tested results are identical whether this runs or
  not.

## Packaging

* Reorganised into a standard R package: documented, namespaced, tested, and
  installable with `R CMD INSTALL`. The command-line entry point survives as
  `system.file("scripts", "cyraven.R", package = "cyRAVEN")` and calls the same
  `run_cyraven()` that is now exported.
* Mutable globals (`COLORS`, `REFERENCE_GROUP`) replaced by a private package
  environment with accessors — `fcs_colors()`/`set_fcs_colors()` and
  `reference_group()`/`set_reference_group()`. A namespace is locked after
  loading, so the previous pattern could not work in a package; assigning into
  the global environment instead is something a package must never do.
* All progress reporting moved from `cat()` to `message()`, so it can be
  suppressed, captured, and redirected. Verbosity is set once per session with
  `options(cyRAVEN.verbose = )` rather than threaded through every call.
* The figure raster ceiling is now `options(cyRAVEN.max_raster_px = )` instead of
  a hard-coded constant.

## Performance

All three changes below are exact: they compute the same numbers by a cheaper
route, and are covered by tests that assert bit-identity rather than closeness.

* **Cofactor derivation is roughly two orders of magnitude faster.** The
  bisection evaluated `asinh()` over the entire background vector on each of its
  40 iterations — around 264 million transcendental calls per sample at 300,000
  events. Because `asinh(./m)` is strictly increasing it maps order statistics to
  order statistics, and R's type-7 quantile interpolates between positions that
  depend only on `n` and `p`, never on `m`. The four order statistics the
  interquartile range needs are therefore selected once, and each iteration costs
  four scalar `asinh()` calls.
* **SOM node updates are O(n) per epoch instead of O(nodes x n).** The per-node
  `colSums(X[mapping == j, ])` loop is replaced by a single `rowsum()` pass,
  removing about 100 million redundant comparisons per `--cluster` run.
* **LISI neighbour selection uses a partial sort.** Only the identity of the k
  nearest neighbours affects the Simpson index, so the full `order()` of every
  distance row was unnecessary.

## Memory

* Each sample's raw expression matrix is released as soon as its transformed
  matrix exists, rather than being held for the whole run. Both were previously
  alive simultaneously during population scoring — which is exactly where an
  uncapped run was killed by the OOM reaper. Unstained controls are held back
  until every sample that could reference them has been scored, and
  `--keep-exprs` still retains everything.

## Fixes

* Long plot subtitles and captions are wrapped to the figure width. `ggplot2`
  does not wrap them, so a long one was drawn on a single line and silently
  clipped at the device edge, mid-sentence, with no warning.
* Heatmap tile labels choose black or white per tile by WCAG contrast ratio
  against the fill they sit on. A single fixed grey was unreadable on the dark
  end of a sequential palette, and which tiles were affected depended on the
  data.
* `--cofactor-from-first-sample` no longer collides with `--cofactor` through R's
  partial matching of `$` on lists. `optparse` omits `NULL`-default options from
  the list it returns, so `opt$cofactor` matched the longer name and returned
  `FALSE`, which was then used as a cofactor of zero: every marker became
  `asinh(x/0) = Inf` and every sample failed staining QC with a message that
  pointed nowhere near the cause.
