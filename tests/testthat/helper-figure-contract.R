# The rendered-file contract every figure has to satisfy, shared by every test
# file that writes a PNG.
#
# WHY IT IS A HELPER. It lived at the top of test-figure-contract.R, which made it
# invisible to every other test file -- testthat evaluates each file in its own
# environment, so a second file testing figures had to either duplicate the header
# reader or assert something weaker. Duplicating it is how two files end up
# disagreeing about what "opaque" means.
#
# Reading the PNG header directly keeps this dependency-free: `png` is not
# installed in the pinned image and adding it to Suggests to run one assertion
# would grow the footprint for no analytical gain.

# IHDR is the first chunk and has a fixed layout, so the fields below sit at
# known offsets: 8-byte signature, 4-byte length, 4-byte type, then width,
# height, bit depth and colour type.
png_header <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  sig <- readBin(con, "raw", 8L)
  if (!identical(sig[1:4], as.raw(c(0x89, 0x50, 0x4E, 0x47))))
    stop("not a PNG: ", path)
  readBin(con, "raw", 8L)                       # chunk length + "IHDR"
  w  <- readBin(con, "integer", 1L, size = 4L, endian = "big")
  h  <- readBin(con, "integer", 1L, size = 4L, endian = "big")
  readBin(con, "raw", 1L)                       # bit depth
  ct <- as.integer(readBin(con, "raw", 1L))
  list(width = w, height = h, colour_type = ct)
}

# Colour type 6 is RGBA, which is how a figure acquires a transparent
# background; 2 is RGB, 3 is a palette, 0 is greyscale, all opaque.
OPAQUE_TYPES <- c(0L, 2L, 3L)
# A rendered dimension below this means a row or panel has collapsed rather
# than merely being small. The overlay bug produced exactly this shape.
FLOOR_PX <- 150L
# The ceiling safe_ggsave() enforces; past it the device silently fails.
CEIL_PX  <- 30000L

expect_figure_contract <- function(path, label) {
  expect_true(file.exists(path), info = paste(label, "was not written"))
  hdr <- png_header(path)
  expect_true(hdr$colour_type %in% OPAQUE_TYPES,
              info = paste0(label, ": colour type ", hdr$colour_type,
                            " (6 = RGBA, i.e. a transparent background)"))
  expect_gte(hdr$width,  FLOOR_PX)
  expect_gte(hdr$height, FLOOR_PX)
  expect_lte(hdr$width,  CEIL_PX)
  expect_lte(hdr$height, CEIL_PX)
  invisible(hdr)
}
