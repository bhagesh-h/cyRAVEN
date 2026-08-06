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
  suppressMessages(write_gating_ml(polys, p, transform = tr, n_per_edge = 4L))

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
