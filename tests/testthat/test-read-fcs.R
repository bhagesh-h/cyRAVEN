# =============================================================================
# Channel resolution
#
# The reading stage decides which columns become markers and which become the
# scatter gate. Both decisions are silent when they go wrong: a file whose
# markers do not resolve produces an empty population table several stages later,
# with a message that blames the specification rather than the panel.
# =============================================================================

#' A parameter block shaped like the one flowCore returns, without needing a
#' file on disk. `read_fcs_resolved()` is not called here because it reads FCS;
#' the classification logic it applies is what these tests exercise, through the
#' same expressions.
resolve_channels <- function(nm, desc) {
  sym <- trimws(sub(" : .*$", "", desc))
  sym[is.na(sym) | !nzchar(sym)] <- nm[is.na(sym) | !nzchar(sym)]

  is_area <- grepl("(^|[^A-Za-z])A$|\\.A$|-A$|Area$", nm) | grepl("- Area$", desc)
  scatter <- c(`FSC-A` = "FSC.*A(rea)?$", `FSC-H` = "FSC.*H(eight)?$",
               `SSC-A` = "SSC.*A(rea)?$", `SSC-H` = "SSC.*H(eight)?$")
  sc_cols <- integer(0)
  for (k in names(scatter)) {
    hit <- which(grepl(scatter[[k]], nm) | grepl(scatter[[k]], desc))
    if (length(hit)) sc_cols[k] <- hit[1]
  }
  fluor <- !grepl("^FSC|^SSC|^Time|Time.Stamp", nm, ignore.case = TRUE) &
           !grepl("^FSC|^SSC|^Time", sym, ignore.case = TRUE)
  height_only <- !any(is_area & fluor) && any(fluor)
  if (height_only) {
    is_h <- grepl("(^|[^A-Za-z])H$|\\.H$|-H$|Height$", nm) | grepl("- Height$", desc)
    if (any(is_h & fluor)) is_area <- is_area | is_h
  }
  for (k in c("FSC", "SSC"))
    if (!paste0(k, "-A") %in% names(sc_cols) && paste0(k, "-H") %in% names(sc_cols))
      sc_cols[paste0(k, "-A")] <- sc_cols[[paste0(k, "-H")]]
  keep <- which(is_area & fluor)
  mk <- setNames(keep, sym[keep])
  list(markers = mk[!duplicated(names(mk))], scatter = sc_cols,
       height_only = height_only)
}

test_that("an area acquisition resolves area channels and ignores height", {
  # The rule that keeps a marker recorded twice from being counted twice.
  r <- resolve_channels(
    nm   = c("FSC-A", "FSC-H", "SSC-A", "CD3-A", "CD3-H", "CD4-A", "Time"),
    desc = c("FSC-A", "FSC-H", "SSC-A", "CD3",   "CD3",   "CD4",   "Time"))
  expect_setequal(names(r$markers), c("CD3", "CD4"))
  expect_false(r$height_only)
  # The height duplicate of CD3 is not a second marker.
  expect_length(r$markers, 2L)
  expect_identical(unname(r$scatter[["FSC-A"]]), 1L)
})

test_that("a height-only acquisition resolves its markers rather than none", {
  # THE REGRESSION. Older instruments and some clinical archives record height
  # only. Applying the area rule unchanged resolved ZERO markers, and the run
  # then failed several stages later reporting every population UNAVAILABLE,
  # which points at the specification instead of the panel.
  r <- resolve_channels(
    nm   = c("FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL3-H", "FL4-H", "Time"),
    desc = c("FSC-Height", "SSC-Height", "CD15", "CD45", "CD14", "CD33", "Time"))
  expect_true(r$height_only)
  expect_setequal(names(r$markers), c("CD15", "CD45", "CD14", "CD33"))
})

test_that("the scatter gate falls back to height when no area channel exists", {
  # derive_scatter_gate() requires FSC-A and SSC-A by name and stops without
  # them, so a height-only file could not be gated at all.
  r <- resolve_channels(
    nm   = c("FSC-H", "SSC-H", "FL1-H", "Time"),
    desc = c("FSC-Height", "SSC-Height", "CD15", "Time"))
  expect_true(all(c("FSC-A", "SSC-A") %in% names(r$scatter)))
  # They point at the height columns, which is what makes the singlet guard
  # below necessary.
  expect_identical(unname(r$scatter[["FSC-A"]]), unname(r$scatter[["FSC-H"]]))
})

test_that("a single area channel suppresses the height fallback", {
  # A file carrying one area channel among otherwise height-only fluorescence is
  # NOT height-only, and the area rule selects that channel alone. This is the
  # state flowCore's GvHD frames are in before FL2-A is dropped, and the reason
  # the demonstration script drops it.
  r <- resolve_channels(
    nm   = c("FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL2-A", "FL4-H", "Time"),
    desc = c("FSC-Height", "SSC-Height", "CD15", "CD45", NA, "CD33", "Time"))
  expect_false(r$height_only)
  expect_length(r$markers, 1L)
  expect_identical(names(r$markers), "FL2-A")
})

test_that("scatter and time are never treated as markers on either path", {
  r <- resolve_channels(
    nm   = c("FSC-H", "SSC-H", "Time", "FL1-H"),
    desc = c("FSC-Height", "SSC-Height", "Time (51.20 sec.)", "CD15"))
  expect_identical(names(r$markers), "CD15")
})

test_that("the singlet gate refuses a ratio of a channel with itself", {
  # With the scatter fallback above, FSC-H and FSC-A name the same column. The
  # ratio would be identically 1: a MAD of zero, a band of zero width, and every
  # event discarded.
  ex <- cbind(FSC = c(100, 200, 300, 400), SSC = c(10, 20, 30, 40))
  sc <- c(`FSC-H` = 1L, `FSC-A` = 1L)
  out <- suppressMessages(derive_singlet_band(ex, sc, rep(TRUE, 4L)))
  expect_true(out$skipped)
  expect_true(all(out$mask))

  # Distinct columns still gate normally.
  ex2 <- cbind(`FSC-A` = c(100, 200, 300, 400), `FSC-H` = c(95, 190, 285, 380))
  sc2 <- c(`FSC-A` = 1L, `FSC-H` = 2L)
  out2 <- suppressMessages(derive_singlet_band(ex2, sc2, rep(TRUE, 4L)))
  expect_false(out2$skipped)
})
