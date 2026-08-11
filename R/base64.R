# SECTION 13a -- BASE64
# =============================================================================
#
# WHY THIS IS WRITTEN OUT RATHER THAN TAKEN FROM A PACKAGE. The run report is a
# single self-contained HTML file, which means every figure and every table has
# to travel inside it as a data URI. That needs a base64 encoder, and the ones
# on CRAN (base64enc, xfun) would each add a dependency to a package whose
# import list is deliberately short and whose container image is pinned. The
# algorithm is fifteen lines and has no edge cases beyond padding.
#
# WHY IT IS VECTORISED OVER THE WHOLE VECTOR. A per-byte loop over a 6 MB PNG is
# minutes in R. Reshaping to a 3-row matrix and indexing an alphabet once turns
# the same work into a handful of vector operations, and rawToChar() assembles
# the result without paste(collapse = "") over millions of elements.

#' Encode a raw vector as base64
#'
#' @param r A raw vector.
#' @return A single character string, padded with `=` as the standard requires.
#' @keywords internal
base64_encode_raw <- function(r) {
  if (!length(r)) return("")
  alphabet <- utf8ToInt(paste0("ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                               "abcdefghijklmnopqrstuvwxyz",
                               "0123456789+/"))
  n <- length(r)
  pad <- (3L - n %% 3L) %% 3L
  if (pad) r <- c(r, as.raw(rep(0L, pad)))
  m <- matrix(as.integer(r), nrow = 3L)
  i1 <- bitwShiftR(m[1L, ], 2L)
  i2 <- bitwOr(bitwShiftL(bitwAnd(m[1L, ], 3L),  4L), bitwShiftR(m[2L, ], 4L))
  i3 <- bitwOr(bitwShiftL(bitwAnd(m[2L, ], 15L), 2L), bitwShiftR(m[3L, ], 6L))
  i4 <- bitwAnd(m[3L, ], 63L)
  v <- alphabet[c(rbind(i1, i2, i3, i4)) + 1L]
  # The padding bytes added above encoded as 'A'; the standard requires '='.
  if (pad) v[seq.int(length(v) - pad + 1L, length(v))] <- utf8ToInt("=")
  rawToChar(as.raw(v))
}

#' Read a file and return it as a data URI
#'
#' @param path File path.
#' @param mime MIME type, defaulting to PNG.
#' @return `data:<mime>;base64,...`, or NA when the file cannot be read.
#' @keywords internal
file_data_uri <- function(path, mime = "image/png") {
  if (!file.exists(path)) return(NA_character_)
  r <- tryCatch(readBin(path, "raw", n = file.size(path)),
                error = function(e) NULL)
  if (is.null(r)) return(NA_character_)
  paste0("data:", mime, ";base64,", base64_encode_raw(r))
}
