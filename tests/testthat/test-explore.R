test_that("explore_features takes every eligible channel by default", {
  f <- cyRAVEN:::explore_features(c("CD3", "CD4", "CD45"), c("FSC-A", "SSC-A"))
  # The whole point of option 1: no lineage preference, scatter included.
  expect_setequal(f, c("CD3", "CD4", "CD45", "FSC-A", "SSC-A"))
})

test_that("explore_features honours an explicit list and exclusions", {
  expect_setequal(cyRAVEN:::explore_features(c("CD3", "CD4"), "SSC-A", prefer = c("CD3", "SSC-A")),
                  c("CD3", "SSC-A"))
  expect_setequal(cyRAVEN:::explore_features(c("CD3", "CD4"), "SSC-A", exclude = "SSC-A"),
                  c("CD3", "CD4"))
  # A requested channel the panel does not have is dropped, not invented.
  expect_setequal(cyRAVEN:::explore_features(c("CD3"), character(0), prefer = c("CD3", "CD99")),
                  "CD3")
})

test_that("two_mode_split separates a clear bimodal set and refuses a flat one", {
  sp <- cyRAVEN:::two_mode_split(c(0.1, 0.2, 0.15, 5.0, 5.2, 4.9))
  expect_false(is.null(sp))
  expect_setequal(sp$high, 4:6)
  expect_true(sp$cut > 0.2 && sp$cut < 4.9)

  expect_null(cyRAVEN:::two_mode_split(rep(1, 6)))   # no variance
  expect_null(cyRAVEN:::two_mode_split(c(1, 2)))     # too few clusters
})

test_that("explore_positivity uses each sample's own threshold", {
  # Two samples, the second stained twice as bright. A single pooled cut would
  # call the whole of sample B positive; per-sample cuts must not.
  X <- matrix(c(0, 1, 0, 1,   0, 2, 0, 2), ncol = 2,
              dimnames = list(NULL, c("CD3", "CD4")))
  cl <- c(1L, 1L, 2L, 2L)
  sid <- c("a", "a", "b", "b")
  thr <- list(a = c(CD3 = 0.5, CD4 = 0.5), b = c(CD3 = 0.5, CD4 = 1.5))

  pos <- cyRAVEN:::explore_positivity(X, cl, sid, thr)
  expect_equal(unname(pos["k1", "CD3"]), 0.5)
  expect_equal(unname(pos["k2", "CD4"]), 0.5)
  expect_equal(unname(attr(pos, "source")[["CD3"]]), "per_sample_threshold")
})

test_that("explore_positivity survives a feature with no threshold", {
  # Scatter channels never have a derived threshold, and `[[` on a NAMED NUMERIC
  # VECTOR raises "subscript out of bounds" for an absent name -- where the same
  # call on a list would return NULL. This is the shape thresholds_used.csv
  # produces, so it has to be the shape the test uses.
  X <- matrix(c(0, 1, 2, 3, 10, 20, 30, 40), ncol = 2,
              dimnames = list(NULL, c("CD3", "SSC-A")))
  thr <- list(a = c(CD3 = 0.5))          # named numeric, no SSC-A entry
  pos <- cyRAVEN:::explore_positivity(X, c(1L, 1L, 2L, 2L), rep("a", 4), thr)
  expect_true(all(is.finite(pos)))
  expect_equal(unname(attr(pos, "source")[["CD3"]]), "per_sample_threshold")
  # The feature with no threshold still gets a value, from the pooled median.
  expect_equal(unname(attr(pos, "source")[["SSC-A"]]), "pooled_median")
})

test_that("explore_positivity falls back to pooled medians without thresholds", {
  X <- matrix(c(0, 1, 2, 3), ncol = 1, dimnames = list(NULL, "CD3"))
  pos <- cyRAVEN:::explore_positivity(X, c(1L, 1L, 2L, 2L), rep("a", 4), list())
  expect_equal(unname(attr(pos, "source")[["CD3"]]), "pooled_median")
  expect_true(all(is.finite(pos)))
})

test_that("explore_phenotype names positives and negatives and admits uncertainty", {
  pos <- matrix(c(0.9, 0.05, 0.5), nrow = 1,
                dimnames = list("k1", c("CD19", "CD3", "CD4")))
  ph <- cyRAVEN:::explore_phenotype(pos)
  expect_true(grepl("CD19\\+", ph))
  expect_true(grepl("CD3-", ph))
  # A marker the cluster is mixed for is OMITTED rather than called negative.
  expect_false(grepl("CD4", ph))

  # Nothing resolvable must not silently produce an empty label.
  expect_equal(unname(cyRAVEN:::explore_phenotype(matrix(0.5, nrow = 1,
                        dimnames = list("k1", "CD3")))), "unresolved")
})

test_that("explore_qc_gate calls debris, dead and saturated", {
  feats <- c("CD45", "LiveDead", "CD3")
  pos <- rbind(k1 = c(0.9, 0.02, 0.8),   # leukocyte, alive        -> keep
               k2 = c(0.05, 0.02, 0.1),  # CD45 low                -> debris
               k3 = c(0.9, 0.95, 0.5),   # viability high          -> dead
               k4 = c(0.99, 0.99, 0.99)) # bright everywhere       -> saturated
  colnames(pos) <- feats
  med <- pos * 3
  g <- cyRAVEN:::explore_qc_gate(pos, med, sizes = c(100, 100, 100, 100),
                       leukocyte = "CD45", viability = "LiveDead",
                       have_thresholds = TRUE)
  calls <- setNames(g$call, g$cluster)
  expect_equal(unname(calls[["k1"]]), "keep")
  expect_equal(unname(calls[["k2"]]), "debris")
  expect_equal(unname(calls[["k3"]]), "dead")
  # Saturation is decided first: an aggregate is bright in the leukocyte marker
  # too, so calling it debris or dead would mislabel it.
  expect_equal(unname(calls[["k4"]]), "saturated")
  expect_true(all(nzchar(g$basis[g$call != "keep"])))
})

test_that("explore_qc_gate keeps the dead cluster when asked", {
  feats <- c("CD45", "LiveDead")
  pos <- rbind(k1 = c(0.9, 0.02), k2 = c(0.9, 0.95))
  colnames(pos) <- feats
  g <- cyRAVEN:::explore_qc_gate(pos, pos * 3, c(50, 50), leukocyte = "CD45",
                       viability = "LiveDead", have_thresholds = TRUE,
                       drop_dead = FALSE)
  expect_true(all(g$call == "keep"))
})

test_that("explore_write_spec emits a curatable YAML, never a silent adoption", {
  d <- withr::local_tempdir()
  pos <- rbind(k1 = c(0.9, 0.02), k2 = c(0.85, 0.9))
  colnames(pos) <- c("CD19", "CD3")
  prof <- data.frame(cluster = c("k1", "k2"), pct_of_gated = c(10, 5))
  nm <- cyRAVEN:::explore_write_spec(pos, prof, c("k1", "k2"), d)
  txt <- readLines(file.path(d, nm))
  expect_true(any(grepl("^populations:", txt)))
  # "above"/"below", the vocabulary score_populations() reads. This asserted
  # "pos"/"neg" until the emitter was fixed, which locked in a file that could
  # never be run: the header invites the reader to curate it and pass it to
  # --config, and doing so reported every population UNAVAILABLE with no error.
  expect_true(any(grepl('"CD19": above', txt, fixed = TRUE)))
  expect_true(any(grepl('"CD3": below', txt, fixed = TRUE)))
  # It must announce that it is a draft, or someone will run it as a spec.
  expect_true(any(grepl("SUGGESTIONS", txt)))

  # THE PROPERTY THAT MATTERS, rather than the spelling: what is written parses
  # as a specification and scores cells. A test on the string alone is what
  # allowed the unrunnable vocabulary to survive.
  y <- yaml::yaml.load_file(file.path(d, nm))
  expect_true(is.list(y$populations) && length(y$populations) >= 1L)
  tmat <- cbind(CD19 = c(5, 1, 5), CD3 = c(1, 5, 5))
  thr  <- c(CD19 = 3, CD3 = 3)
  got <- score_populations(tmat, thr, parent = rep(TRUE, 3), spec = y$populations)
  expect_true(length(got$masks) >= 1L)
  # k1 is CD19-above CD3-below: exactly the first cell.
  expect_equal(sum(got$masks[["k1"]]), 1L)
  expect_true(is.null(got$unavailable[["k1"]]))
})
