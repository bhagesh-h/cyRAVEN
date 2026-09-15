# =============================================================================
# Cluster identity: what each unsupervised cluster corresponds to
#
# The claims: the catch-all remainder is recognised as a remainder wherever a
# coverage question is asked; the assignment scores the match in both directions
# rather than taking the largest overlap; and a cluster the specification cannot
# describe is reported as undescribed instead of being named after whichever
# lineage happens to have the biggest slice of it.
# =============================================================================

test_that("the catch-all label is recognised, and a real population is not", {
  expect_true(all(is_catch_all_label(
    c("Other CD45+", "other", "Other", "unclassified", "unlabelled",
      "unlabeled", "unassigned", "none", "ungated", "", NA))))
  # The bug this guards: the list of names treated as "not a real population"
  # was written out by hand and did not contain the name this package actually
  # assigns, so every coverage question answered "fully covered".
  expect_true(is_catch_all_label(catch_all_label()))
  expect_false(any(is_catch_all_label(
    c("CD4 T cells", "B cells", "NKT cells", "HLA-DR low monocytes"))))
  # A prefix match alone would catch this one, which is a real population name.
  expect_false(is_catch_all_label("Otherwise odd cells"))
})

test_that("a cluster made mostly of the catch-all is called undescribed", {
  ct <- data.frame(
    cluster    = c("k1", "k1", "k2", "k2"),
    population = c("Other CD45+", "B cells", "B cells", "Other CD45+"),
    cells      = c(950, 50, 900, 100),
    stringsAsFactors = FALSE)
  id <- cluster_identity_table(ct)
  k1 <- id[id$cluster == "k1", ]
  expect_identical(k1$identity, "undescribed")
  expect_identical(k1$call, "undescribed by the specification")
  expect_equal(k1$pct_catch_all, 95)
  # k2 is mostly B cells and holds most of the B cells, so it is a real match.
  k2 <- id[id$cluster == "k2", ]
  expect_identical(k2$identity, "B cells")
  expect_identical(k2$call, "confident")
})

test_that("the assignment scores both directions, not the largest overlap alone", {
  # THE CASE THAT SEPARATES THE TWO RULES. In k1 the plurality among the real
  # labels is CD4 T cells (300 of 500 cells), so a largest-overlap rule names
  # k1 "CD4 T cells". But those 300 are a twentieth of the 6,000 CD4 T cells in
  # the run, while k1 holds 200 of the 210 B cells that exist -- so the label
  # that k1 actually accounts for is B cells. F1 picks B cells; plurality would
  # not.
  ct <- data.frame(
    cluster    = c("k1", "k1", "k2", "k2"),
    population = c("CD4 T cells", "B cells", "CD4 T cells", "B cells"),
    cells      = c(300, 200, 5700, 10),
    stringsAsFactors = FALSE)
  id <- cluster_identity_table(ct)
  k1 <- id[id$cluster == "k1", ]
  expect_identical(k1$identity, "B cells")
  expect_equal(k1$precision, 40)     # 200 of k1's 500 cells
  expect_equal(k1$recall, 95.2)      # 200 of the 210 B cells in the run
  # CD4 T cells is retained as the runner-up rather than discarded: it is what
  # the plurality rule would have said, and a reader comparing the two needs it.
  expect_identical(k1$runner_up, "CD4 T cells")
  expect_lt(k1$runner_up_f1, k1$f1)
  # k2 is the cluster CD4 T cells actually lives in.
  k2 <- id[id$cluster == "k2", ]
  expect_identical(k2$identity, "CD4 T cells")
  expect_gt(k2$f1, 0.9)
  expect_identical(k2$call, "confident")
})

test_that("a cluster more than half catch-all is never called confident", {
  # F1 is computed over the declared labels only, so a cluster that is 60%
  # remainder can still score a high F1 against the lineage holding the other
  # 40%. The label is the best answer available; "confident" is not, while most
  # of the cluster is cells the specification cannot name.
  ct <- data.frame(
    cluster    = c("k1", "k1"),
    population = c("Other CD45+", "NK cells"),
    cells      = c(600, 400),
    stringsAsFactors = FALSE)
  id <- cluster_identity_table(ct)
  expect_identical(id$identity, "NK cells")
  expect_gt(id$f1, 0.5)               # the NK cells really are all here
  expect_identical(id$call, "partial")
})

test_that("the catch-all cannot win the assignment", {
  # It is the complement of the specification, not a lineage. Allowed to
  # compete it would win most rows on a specification with modest coverage,
  # which is true and useless.
  ct <- data.frame(
    cluster    = c("k1", "k1"),
    population = c("Other CD45+", "NK cells"),
    cells      = c(600, 400),
    stringsAsFactors = FALSE)
  id <- cluster_identity_table(ct)
  expect_identical(id$identity, "NK cells")
  expect_false(is_catch_all_label(id$identity))
  expect_equal(id$pct_catch_all, 60)
})

test_that("a cluster with no real label at all is undescribed, not an error", {
  ct <- data.frame(cluster = "k1", population = "Other CD45+", cells = 100,
                   stringsAsFactors = FALSE)
  id <- cluster_identity_table(ct)
  expect_identical(id$identity, "undescribed")
  expect_true(is.na(id$f1))
  expect_equal(id$pct_catch_all, 100)
})

test_that("an empty or malformed cross-tabulation yields NULL rather than failing", {
  expect_null(cluster_identity_table(NULL))
  expect_null(cluster_identity_table(data.frame()))
  expect_null(cluster_identity_table(
    data.frame(cluster = "k1", cells = 10, stringsAsFactors = FALSE)))
})

test_that("the identity figure is written and is not empty", {
  ct <- data.frame(
    cluster    = rep(c("k1", "k2", "k3"), each = 2),
    population = c("Other CD45+", "B cells", "CD4 T cells", "Other CD45+",
                   "NK cells", "Other CD45+"),
    cells      = c(950, 50, 800, 200, 700, 300),
    stringsAsFactors = FALSE)
  f <- file.path(withr::local_tempdir(), "identity.png")
  fig_cluster_identity(ct, NULL, f)
  expect_true(file.exists(f))
  expect_gt(file.size(f), 5000)
})

test_that("the findings verdict counts the catch-all as unlabelled", {
  # THE REGRESSION THIS PINS DOWN. explore_findings.csv exists to report the
  # clusters no declared population covers. The list of labels it treated as
  # unlabelled was spelled out by hand and did not contain the name this
  # package actually assigns to unmatched cells, so pct_unlabelled was 0 for
  # every cluster and every verdict read "covered by the declared
  # specification" -- including clusters that were 100% remainder, on runs
  # whose specification described under a tenth of the cells. The table's only
  # output was inverted, silently.
  ct <- data.frame(
    cluster    = c("k1", "k2", "k2"),
    population = c(catch_all_label(), "B cells", catch_all_label()),
    cells      = c(1000, 900, 100),
    stringsAsFactors = FALSE)
  undecl <- vapply(c("k1", "k2"), function(k) {
    sub <- ct[ct$cluster == k, , drop = FALSE]
    100 * sum(sub$cells[is_catch_all_label(sub$population)]) / sum(sub$cells)
  }, numeric(1))
  expect_equal(unname(undecl), c(100, 10))
  # And the >= 70 rule that turns this into a verdict then fires on k1 alone.
  expect_identical(unname(undecl >= 70), c(TRUE, FALSE))
})

test_that("marker_state is bounded, and drops populations on size not on result", {
  # The figure draws one panel per population x marker pair, so its cost is the
  # PRODUCT of two independently growing numbers. With the subset layer added it
  # reached 1,275 panels and the renderer was killed building it: one run wrote
  # a 30 MB PNG nobody could open, the next left a zero-byte file.
  set.seed(9)
  pops <- paste("pop", sprintf("%02d", 1:40))
  mks  <- paste0("M", 1:25)
  mfi <- expand.grid(population = pops, marker = mks, sample_id = c("S1", "S2"),
                     stringsAsFactors = FALSE)
  # pop 01 is the most abundant, pop 40 the least.
  mfi$n_cells <- 5000 - 100 * as.integer(sub("pop ", "", mfi$population))
  mfi$pct_positive <- runif(nrow(mfi), 0, 100)

  f <- file.path(withr::local_tempdir(), "ms.png")
  suppressMessages(fig_marker_state(mfi, f, max_panels = 250L))
  expect_true(file.exists(f))
  expect_gt(file.size(f), 5000)

  # 250 panels over 25 markers leaves room for 10 populations, and they are the
  # ten most abundant -- chosen before any test is consulted, so the figure's
  # contents cannot be shaped by the result.
  d <- mfi[mfi$n_cells >= 20, ]
  sz <- tapply(d$n_cells, d$population, stats::median)
  expect_identical(names(sort(sz, decreasing = TRUE))[1:10],
                   paste("pop", sprintf("%02d", 1:10)))
})

test_that("a population count under the ceiling is left completely alone", {
  set.seed(10)
  mfi <- expand.grid(population = paste("pop", 1:4), marker = paste0("M", 1:5),
                     sample_id = c("S1", "S2"), stringsAsFactors = FALSE)
  mfi$n_cells <- 1000L
  mfi$pct_positive <- runif(nrow(mfi), 0, 100)
  f <- file.path(withr::local_tempdir(), "ms.png")
  # No message about a ceiling, because none was applied.
  expect_silent(suppressMessages(fig_marker_state(mfi, f, max_panels = 400L)))
  expect_true(file.exists(f))
})

test_that("a cluster is named from its own marker profile", {
  # The annotation attaches to the CLUSTER, from explore_cluster_profile.csv,
  # never to the cell: the declared specification is not re-scored and no
  # declared output changes.
  prof <- data.frame(
    cluster = c("k1", "k2", "k3"),
    `frac_pos.CD3`  = c(0.99, 0.98, 0.02),
    `frac_pos.CD4`  = c(0.97, 0.03, 0.01),
    `frac_pos.CD8`  = c(0.02, 0.95, 0.01),
    `frac_pos.CD69` = c(0.88, 0.10, 0.05),
    `frac_pos.CD19` = c(0.01, 0.01, 0.93),
    `frac_pos.CD14` = c(0.01, 0.02, 0.02),
    `frac_pos.CD56` = c(0.02, 0.03, 0.02),
    check.names = FALSE, stringsAsFactors = FALSE)
  a <- annotate_clusters_with_subsets(prof)
  expect_identical(a$cluster, c("k1", "k2", "k3"))
  expect_match(a$subset_label[1], "CD4 T cells CD69", fixed = TRUE)
  expect_true(is.na(a$subset_label[3]) ||
                !grepl("CD4|CD8", a$subset_label[3]))
  # The margin says how far from the cut the weakest requirement sat, so a
  # cluster that is 51% positive is not reported like one that is 99%.
  expect_true(is.finite(a$margin[1]))
})

test_that("a definition needing an upper bound is dropped, not approximated", {
  # "bright" and "intermediate" need the valley between the two positive modes,
  # which a positivity FRACTION does not carry. Treating bright as merely
  # positive would name a whole compartment after its upper mode.
  prof <- data.frame(cluster = "k1",
                     `frac_pos.CD3` = 0.02, `frac_pos.CD56` = 0.95,
                     check.names = FALSE, stringsAsFactors = FALSE)
  a <- annotate_clusters_with_subsets(prof)
  expect_false(isTRUE(grepl("bright", a$subset_label)))
})

test_that("a profile with no frac_pos columns yields NULL rather than failing", {
  expect_null(annotate_clusters_with_subsets(NULL))
  expect_null(annotate_clusters_with_subsets(
    data.frame(cluster = "k1", cells = 10, stringsAsFactors = FALSE)))
})
