# explore_matrix(rows =) must be an exact refactor, not an approximation.
#
# It exists so that --explore-cells-per-sample bounds PEAK MEMORY and not just
# the size of the result. Transforming every event and subsetting afterwards
# materialised the whole file first: with --max-events-per-file 0, a 4-million
# event acquisition is a ~770 MB matrix before one cell is discarded, and a run
# whose declared half completed was killed with SIGKILL partway through explore.
#
# The property that makes the change safe is that every transform is elementwise
# with parameters fixed at construction. If someone later adds a transform that
# derives a parameter from the data it is handed, these tests fail, which is the
# point.

make_rd <- function(n = 2000L, seed = 7L) {
  set.seed(seed)
  ex <- matrix(rexp(n * 6, 1 / 500), ncol = 6)
  colnames(ex) <- paste0("ch", 1:6)
  list(exprs = ex,
       marker_cols  = setNames(as.list(1:4), c("CD3", "CD4", "HLA-DR", "TCR-Vd1")),
       scatter_cols = setNames(as.list(5:6), c("FSC-A", "SSC-A")))
}
FEATS <- c("CD3", "CD4", "HLA-DR", "TCR-Vd1", "FSC-A", "SSC-A")

test_that("subset-then-transform equals transform-then-subset, exactly", {
  rd <- make_rd()
  tr <- make_transform("arcsinh", cofactor = 197.1)
  take <- sort(sample.int(nrow(rd$exprs), 300))

  old <- explore_matrix(rd, FEATS, tr)[take, , drop = FALSE]
  new <- explore_matrix(rd, FEATS, tr, rows = take)

  # identical(), not all.equal(): a tolerance would let a real change through.
  expect_identical(old, new)
  expect_identical(colnames(old), colnames(new))
})

test_that("it holds for the identity transform too", {
  rd <- make_rd()
  tr <- make_transform("none")
  take <- sort(sample.int(nrow(rd$exprs), 111))
  expect_identical(explore_matrix(rd, FEATS, tr)[take, , drop = FALSE],
                   explore_matrix(rd, FEATS, tr, rows = take))
})

test_that("rows = NULL keeps the old whole-file behaviour", {
  rd <- make_rd()
  tr <- make_transform("arcsinh", cofactor = 150)
  m <- explore_matrix(rd, FEATS, tr)
  expect_equal(nrow(m), nrow(rd$exprs))
  expect_identical(m, explore_matrix(rd, FEATS, tr, rows = seq_len(nrow(rd$exprs))))
})

test_that("a single retained row still returns a matrix, not a vector", {
  # vapply() collapses to a vector at n = 1; the guard in explore_matrix keeps
  # it a matrix. Without that, the caller's as.data.frame() produces one column
  # of six values instead of one row of six markers.
  rd <- make_rd()
  tr <- make_transform("arcsinh", cofactor = 150)
  m <- explore_matrix(rd, FEATS, tr, rows = 5L)
  expect_true(is.matrix(m))
  expect_equal(nrow(m), 1L)
  expect_equal(ncol(m), length(FEATS))
})

test_that("scatter-only and marker-only feature sets both subset correctly", {
  rd <- make_rd()
  tr <- make_transform("arcsinh", cofactor = 150)
  take <- sort(sample.int(nrow(rd$exprs), 50))
  for (f in list(c("CD3", "CD4"), c("FSC-A", "SSC-A"))) {
    expect_identical(explore_matrix(rd, f, tr)[take, , drop = FALSE],
                     explore_matrix(rd, f, tr, rows = take))
  }
})
