# =============================================================================
# Labels from another tool: joining them, gating them, and writing gates out
# =============================================================================

#' Cells from several donors with one planted population, so a gate that works
#' has something specific to find and a gate that memorises donors will show it.
planted_donor_cells <- function(n_donor = 5L, per = 600L, seed = 31) {
  withr::local_seed(seed)
  do.call(rbind, lapply(seq_len(n_donor), function(d) {
    lab <- sample(c("target", "other"), per, TRUE, prob = c(0.3, 0.7))
    data.frame(
      sample_id = sprintf("D%02d", d), event_index = seq_len(per),
      external_label = lab,
      CD3 = ifelse(lab == "target", stats::rnorm(per, 4, 0.5),
                   stats::rnorm(per, 0.5, 0.5)),
      CD4 = ifelse(lab == "target", stats::rnorm(per, 3.5, 0.5),
                   stats::rnorm(per, 0.5, 0.5)),
      CD8 = stats::rnorm(per, 1, 1), CD19 = stats::rnorm(per, 1, 1),
      stringsAsFactors = FALSE)
  }))
}

test_that("label columns are matched by convention and missing ones named", {
  p <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = "S1", event = 1:5, cluster = "c1"),
                   p, row.names = FALSE)
  d <- suppressMessages(read_external_labels(p))
  expect_identical(names(d), c("sample_id", "event_index", "label"))
  expect_identical(nrow(d), 5L)

  q <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(data.frame(a = 1, b = 2), q, row.names = FALSE)
  expect_error(suppressMessages(read_external_labels(q)), "no sample column")
})

test_that("the join is on the event key, so a subsampled table cannot mislabel", {
  cells <- planted_donor_cells()
  cells$external_label <- NULL
  lab <- data.frame(sample_id = cells$sample_id, event_index = cells$event_index,
                    label = "target", stringsAsFactors = FALSE)
  # Shuffle the label table: a positional join would now be wrong everywhere,
  # a keyed one is unaffected.
  lab <- lab[sample.int(nrow(lab)), ]
  lab$label <- ifelse(lab$event_index %% 2L == 0L, "even", "odd")

  out <- suppressMessages(join_external_labels(cells, lab))
  expect_identical(out$external_label,
                   ifelse(cells$event_index %% 2L == 0L, "even", "odd"))
})

test_that("a join that matches almost nothing returns NULL rather than a gate", {
  cells <- planted_donor_cells()
  lab <- data.frame(sample_id = "NOT_A_SAMPLE", event_index = 1:500,
                    label = "x", stringsAsFactors = FALSE)
  expect_null(suppressMessages(join_external_labels(cells, lab)))
})

test_that("applying a strategy reproduces the selection it was fitted with", {
  # This is what makes a proposal executable: the same polygons, run again by a
  # different function, must select the same cells.
  cells <- planted_donor_cells()
  X <- as.matrix(cells[, c("CD3", "CD4", "CD8", "CD19")])
  y <- as.integer(cells$external_label == "target")
  st <- suppressMessages(explain_cluster(X, y, max_depth = 2L))
  skip_if(is.null(st), "no gate fitted on this fixture")

  inside <- apply_gate_strategy(X, st)
  expect_length(inside, nrow(X))
  # The gate finds the planted population rather than a slice of the background.
  expect_gt(mean(y[inside] == 1), 0.8)
  expect_gt(sum(inside & y == 1) / sum(y == 1), 0.5)

  # Truncating the strategy can only keep more cells, never fewer.
  expect_gte(sum(apply_gate_strategy(X, st, depth = 1L)), sum(inside))
})

test_that("transferability scores each donor and refuses below three", {
  cells <- planted_donor_cells()
  X <- as.matrix(cells[, c("CD3", "CD4", "CD8", "CD19")])
  y <- as.integer(cells$external_label == "target")

  tf <- suppressMessages(gate_transferability(X, y, cells$sample_id,
                                              max_depth = 2L, seed = 3))
  skip_if(is.null(tf), "no gate fitted on this fixture")
  expect_identical(nrow(tf$per_donor), length(unique(cells$sample_id)))
  expect_true(all(c("f1_min", "f1_median", "f1_max") %in% names(tf$summary)))
  # The population is planted identically in every donor, so a gate learned
  # without one donor must still find it in that donor.
  expect_gt(tf$summary$f1_min, 0.5)

  two <- cells[cells$sample_id %in% c("D01", "D02"), ]
  expect_null(suppressMessages(gate_transferability(
    as.matrix(two[, c("CD3", "CD4", "CD8", "CD19")]),
    as.integer(two$external_label == "target"), two$sample_id)))
})

test_that("transferability leaves the RNG stream where it found it", {
  cells <- planted_donor_cells()
  X <- as.matrix(cells[, c("CD3", "CD4", "CD8", "CD19")])
  y <- as.integer(cells$external_label == "target")

  set.seed(99)
  before <- .Random.seed
  expected <- stats::runif(4)
  set.seed(99)
  invisible(suppressMessages(gate_transferability(X, y, cells$sample_id,
                                                  max_depth = 2L, seed = 3)))
  expect_identical(.Random.seed, before)
  expect_identical(stats::runif(4), expected)
})

test_that("the arcsinh inverse round-trips and densification tracks the edge", {
  tr <- make_transform("arcsinh", cofactor = 150)
  y <- c(-1, 0, 0.5, 2, 4.5)
  expect_equal(tr$fn(tr$inv(y)), y, tolerance = 1e-9)

  # A straight edge on the analysis scale is a curve in linear units. The
  # subdivided boundary must therefore have more vertices than the corners, and
  # every one of them must sit back on the fitted edge when transformed again.
  poly <- cbind(c(1, 4, 4, 1), c(1, 1, 3, 3))
  lin <- polygon_to_linear(poly, tr, "CD3", "CD4", n_per_edge = 8L)
  expect_identical(nrow(lin), 32L)
  back <- cbind(tr$fn(lin[, 1]), tr$fn(lin[, 2]))
  # Vertex 1 and vertex 9 are the first two corners of the original polygon.
  expect_equal(back[1, ], poly[1, ], tolerance = 1e-6)
  expect_equal(back[9, ], poly[2, ], tolerance = 1e-6)
})

test_that("Gating-ML output is well formed and chains the levels", {
  tr <- make_transform("arcsinh", cofactor = 150)
  polys <- do.call(rbind, lapply(1:2, function(d) data.frame(
    label = "target", depth = d, marker_x = "CD3", marker_y = "CD4",
    vertex = 1:4, x_transformed = c(1, 4, 4, 1) + d,
    y_transformed = c(1, 1, 3, 3), stringsAsFactors = FALSE)))
  p <- withr::local_tempfile(fileext = ".xml")
  suppressMessages(write_gating_ml(polys, p, transform = tr, n_per_edge = 4L,
                                   space = "linear"))

  x <- paste(readLines(p), collapse = "\n")
  expect_match(x, "Gating-ML/v2.0/gating")
  expect_identical(lengths(regmatches(x, gregexpr("<gating:PolygonGate", x)))[[1]], 2L)
  # Level 2 declares level 1 as its parent, which is what makes the file a
  # strategy rather than two unrelated regions.
  expect_match(x, 'gating:id="target_L2" gating:parent_id="target_L1"')
  # 4 edges x 4 subdivisions x 2 gates.
  expect_identical(lengths(regmatches(x, gregexpr("<gating:vertex>", x)))[[1]], 32L)

  tab <- polygons_linear_table(polys, tr, n_per_edge = 4L)
  expect_identical(nrow(tab), 32L)
  expect_true(all(c("x_linear", "y_linear") %in% names(tab)))
})

test_that("XML special characters in a label cannot break the document", {
  tr <- make_transform("arcsinh", cofactor = 150)
  polys <- data.frame(label = "CD4<T> & \"cells\"", depth = 1L,
                      marker_x = "HLA-DR", marker_y = "CD4", vertex = 1:4,
                      x_transformed = c(1, 4, 4, 1),
                      y_transformed = c(1, 1, 3, 3), stringsAsFactors = FALSE)
  p <- withr::local_tempfile(fileext = ".xml")
  suppressMessages(write_gating_ml(polys, p, transform = tr, n_per_edge = 2L))
  x <- paste(readLines(p), collapse = "\n")
  expect_false(grepl("<T>", x, fixed = TRUE))
  expect_match(x, 'gating:id="CD4_T_____cells__L1"', fixed = TRUE)
})

test_that("Gating-ML on the analysis scale states the fitted polygon exactly", {
  # asinh(x / c) is Gating-ML's fasinh with T = c sinh(M ln 10), A = 0, times
  # M ln 10. Checked against flowCore's own Gating-ML 2.0 implementation.
  tr <- make_transform("arcsinh", cofactor = 150, per_marker = c(CD4 = 400))
  polys <- data.frame(label = "g", depth = 1L, marker_x = "CD3", marker_y = "CD4",
                      vertex = 1:4, x_transformed = c(1, 4, 4, 1),
                      y_transformed = c(1, 1, 3, 3), stringsAsFactors = FALSE)
  p <- withr::local_tempfile(fileext = ".xml")
  suppressMessages(write_gating_ml(polys, p, transform = tr))
  x <- paste(readLines(p), collapse = "\n")
  expect_identical(lengths(regmatches(x, gregexpr("<gating:vertex>", x)))[[1]], 4L)
  expect_match(x, 'gating:compensation-ref="FCS" gating:transformation-ref="T_CD4"')
  t4 <- as.numeric(sub('.*transforms:fasinh transforms:T="([^"]+)".*', "\\1",
                       regmatches(x, regexpr('T_CD4">\n\\s*<transforms:fasinh[^/]*', x))))
  expect_equal(t4, 400 * sinh(4.5 * log(10)), tolerance = 1e-8)
  co <- as.numeric(sub('.*value="([^"]+)".*', "\\1",
                       regmatches(x, gregexpr('data-type:value="[^"]+"', x))[[1]]))
  expect_equal(co[1:2], c(1, 1) / (4.5 * log(10)), tolerance = 1e-8)
  # the round trip through flowCore's implementation
  skip_if_not_installed("flowCore")
  g <- flowCore::asinhtGml2(parameters = "CD4", T = t4, M = 4.5, A = 0)
  xv <- tr$inv(3, "CD4")
  expect_equal(as.numeric(flowCore::eval(g)(cbind(CD4 = xv))) * 4.5 * log(10), 3,
               tolerance = 1e-8)
})

test_that("a logicle gate states its polygon in Gating-ML units for any A", {
  # flowCore's logicle reaches M at T and Gating-ML's reaches 1, whatever A is,
  # so the unit is M; M + A was low by A / (M + A).
  skip_if_not_installed("flowCore")
  pr <- list(w = 0.5, t = 262144, m = 4.5, a = 0.3)
  tr <- make_transform("logicle", logicle = list(CD3 = pr, CD4 = pr))
  polys <- data.frame(label = "g", depth = 1L, marker_x = "CD3", marker_y = "CD4",
                      vertex = 1:4, x_transformed = c(1, 4, 4, 1),
                      y_transformed = c(1, 1, 3, 3), stringsAsFactors = FALSE)
  p <- withr::local_tempfile(fileext = ".xml")
  suppressMessages(write_gating_ml(polys, p, transform = tr))
  x <- paste(readLines(p), collapse = "\n")
  co <- as.numeric(sub('.*value="([^"]+)".*', "\\1",
                       regmatches(x, gregexpr('data-type:value="[^"]+"', x))[[1]]))
  g <- flowCore::logicletGml2(parameters = "CD3", T = pr$t, W = pr$w, M = pr$m, A = pr$a)
  xv <- tr$inv(4, "CD3")
  expect_equal(co[3], as.numeric(flowCore::eval(g)(cbind(CD3 = xv))), tolerance = 1e-4)
})

test_that("gates on data compensated with --spillover carry that matrix", {
  # Gating-ML 2.0 section 7: a gate references a spectrumMatrix by its id and
  # names the compensated dimension by its fluorochrome, and no fluorochrome
  # may share a detector's name (7.2). "FCS" would tell the reader to apply each
  # file's own $SPILLOVER instead. Reads are built as read_fcs_resolved() builds
  # them: alter.names spelling in channel_names, $PnN in the keywords.
  sp <- matrix(c(1, 0.1, 0.02, 0.01,  0.2, 1, 0.05, 0,  0, 0.03, 1, 0,
                 0, 0, 0, 1), 4, byrow = TRUE,
               dimnames = list(NULL, c("FITC.A", "PE.A", "APC.A", "PerCP.A")))
  reads <- list(
    s1 = list(compensation = "compensated, supplied matrix (--spillover)",
              channel_names = c("FSC.A", "FITC.A", "PE.A", "BV421.A", "PerCP.A"),
              keywords = list(`$P1N` = "FSC-A", `$P2N` = "FITC-A", `$P3N` = "PE-A",
                              `$P4N` = "BV421-A", `$P5N` = "PerCP-A"),
              marker_cols = c(CD3 = 2L, CD4 = 3L, CD8 = 4L)))
  expect_identical(gating_ml_markers(reads),
                   c(CD3 = "FITC-A", CD4 = "PE-A", CD8 = "BV421-A"))
  expect_identical(gating_ml_compensation(reads), "supplied")
  reads$s1$compensation <- "compensated, file $SPILLOVER"
  expect_identical(gating_ml_compensation(reads), "FCS")
  reads$s1$compensation <- "compensated, supplied matrix (--spillover)"
  # APC-A is in no file, so the matrix the run applied is 3 x 3, its detectors
  # spelled as the files' $PnN though the matrix header used make.names
  sm <- gating_ml_spillover(reads, sp)
  expect_identical(colnames(sm), c("FITC-A", "PE-A", "PerCP-A"))
  expect_equal(unname(sm), unname(sp[c(1, 2, 4), c(1, 2, 4)]))

  tr <- make_transform("arcsinh", cofactor = 150)
  polys <- data.frame(label = "g", depth = 1L, marker_x = "CD3",
                      marker_y = "CD4", vertex = 1:4, x_transformed = c(1, 4, 4, 1),
                      y_transformed = c(1, 1, 3, 3), stringsAsFactors = FALSE)
  polys2 <- transform(polys, label = "h", marker_y = "CD8")
  p <- withr::local_tempfile(fileext = ".xml")
  suppressMessages(write_gating_ml(rbind(polys, polys2), p, transform = tr,
                                   channel_map = gating_ml_markers(reads),
                                   compensation_ref = "supplied", spillover = sm))
  x <- paste(readLines(p), collapse = "\n")
  expect_match(x, '<transforms:spectrumMatrix transforms:id="cyRAVEN_spillover">')
  # fluorochromes by marker, Comp- where the detector carries none; detectors by $PnN
  expect_match(x, paste0('<transforms:fluorochromes>\\s*',
                         '<data-type:fcs-dimension data-type:name="CD3"/>\\s*',
                         '<data-type:fcs-dimension data-type:name="CD4"/>\\s*',
                         '<data-type:fcs-dimension data-type:name="Comp-PerCP-A"/>'))
  expect_match(x, paste0('<transforms:detectors>\\s*',
                         '<data-type:fcs-dimension data-type:name="FITC-A"/>\\s*',
                         '<data-type:fcs-dimension data-type:name="PE-A"/>\\s*',
                         '<data-type:fcs-dimension data-type:name="PerCP-A"/>'))
  expect_identical(lengths(regmatches(x, gregexpr("<transforms:spectrum>", x)))[[1]], 3L)
  expect_match(x, 'transforms:coefficient transforms:value="0.1"')
  expect_false(grepl('compensation-ref="FCS"', x))
  expect_match(x, paste0('compensation-ref="cyRAVEN_spillover" gating:transformation-ref="T_CD3">\\s*',
                         '<data-type:fcs-dimension data-type:name="CD3"/>'))
  # CD8's detector is not in the matrix, so it was never compensated, and an
  # uncompensated dimension is named by its $PnN
  expect_match(x, paste0('compensation-ref="uncompensated" gating:transformation-ref="T_CD8">\\s*',
                         '<data-type:fcs-dimension data-type:name="BV421-A"/>'))
  expect_error(write_gating_ml(polys, p, transform = tr, compensation_ref = "supplied"),
               "needs the spillover matrix and channel_map")
})

test_that("a run's gates name each dimension by the files' $PnN", {
  # Gating-ML 2.0 section 4.1, and 4.2.2 for compensation-ref "FCS": a conformant
  # reader finds the column by $PnN; no FCS file holds a column named CD3.
  reads <- list(s1 = list(channel_names = c("FITC.A", "PE.A"),
                          keywords = list(`$P1N` = "FITC-A", `$P2N` = "PE-A"),
                          marker_cols = c(CD3 = 1L, CD4 = 2L)))
  tr <- make_transform("arcsinh", cofactor = 150)
  polys <- data.frame(label = "g", depth = 1L, marker_x = "CD3", marker_y = "CD4",
                      vertex = 1:4, x_transformed = c(1, 4, 4, 1),
                      y_transformed = c(1, 1, 3, 3), stringsAsFactors = FALSE)
  p <- withr::local_tempfile(fileext = ".xml")
  suppressMessages(write_gating_ml(polys, p, transform = tr,
                                   channel_map = gating_ml_markers(reads)))
  x <- paste(readLines(p), collapse = "\n")
  expect_match(x, paste0('compensation-ref="FCS" gating:transformation-ref="T_CD3">\\s*',
                         '<data-type:fcs-dimension data-type:name="FITC-A"/>'))
  expect_false(grepl('data-type:name="CD3"', x, fixed = TRUE))
  # a marker symbol equal to a detector name gets Comp- under a matrix
  sm <- matrix(c(1, 0.1, 0.2, 1), 2, dimnames = list(NULL, c("FITC-A", "PE-A")))
  reads$s1$marker_cols <- c(`FITC-A` = 1L, CD4 = 2L)
  polys$marker_x <- "FITC-A"
  suppressMessages(write_gating_ml(polys, p, transform = tr,
                                   channel_map = gating_ml_markers(reads),
                                   compensation_ref = "supplied", spillover = sm))
  x <- paste(readLines(p), collapse = "\n")
  expect_match(x, '<transforms:fluorochromes>\\s*<data-type:fcs-dimension data-type:name="Comp-FITC-A"/>')
  expect_match(x, 'transformation-ref="T_FITC-A">\\s*<data-type:fcs-dimension data-type:name="Comp-FITC-A"/>')
})

test_that("a marker missing from channel_map keeps its own name", {
  tr <- make_transform("arcsinh", cofactor = 150)
  polys <- data.frame(label = "g", depth = 1L, marker_x = "CD3", marker_y = "CD4",
                      vertex = 1:4, x_transformed = c(1, 4, 4, 1),
                      y_transformed = c(1, 1, 3, 3), stringsAsFactors = FALSE)
  p <- withr::local_tempfile(fileext = ".xml")
  suppressMessages(write_gating_ml(polys, p, transform = tr,
                                   channel_map = c(CD3 = "FITC-A")))
  x <- paste(readLines(p), collapse = "\n")
  expect_match(x, 'data-type:name="FITC-A"')
  expect_match(x, 'data-type:name="CD4"')
  expect_false(grepl('data-type:name="NA"', x, fixed = TRUE))
})
