# SECTION 10b -- WRITING A LEARNED GATE OUT TO THE INSTRUMENT
# =============================================================================
#
# WHY A GATE FILE AND NOT A FIGURE. cluster_gate_strategy_*.png shows a
# cytometrist where to draw. That is enough to reproduce a gate by hand and not
# enough to sort on one, apply it to two hundred files, or hand it to a
# collaborator without also handing over this package. Gating-ML 2.0 is the ISAC
# exchange format for exactly this (Spidlen et al. 2015, Cytometry A 87:683) and
# is read by Cytobank, FlowRepository, and FlowJo through an ACS archive.
#
# WHY THE GATES ARE WRITTEN ON THE ANALYSIS SCALE. The polygons are fitted on
# the transformed scale, where their edges are straight. Gating-ML 2.0 carries
# that scale: its parametrised inverse hyperbolic sine and its logicle are the
# transforms cyRAVEN applies, written with other parameters, so a gate can be
# stated exactly where it was fitted (Spidlen et al. 2015; Gating-ML 2.0,
# section 6). asinh(x / c) is fasinh(x; T, M, A = 0) times M ln 10 with
# T = c sinh(M ln 10), and flowCore's logicle is Gating-ML's times M + A; both
# were checked against flowCore's asinhtGml2() and logicletGml2(). Vertex
# coordinates are divided by those factors and nothing is approximated.
#
# `space = "linear"` writes the earlier form instead, for software that does
# not read Gating-ML transforms: the region in the units the FCS file stores,
# where a straight edge on the analysis scale is a CURVE. Each edge is then
# subdivided before inversion, so the emitted boundary tracks the fitted one to
# well within the width of a histogram bin.
#
# NO NEW HARD DEPENDENCY. The XML is written directly. Building a GatingSet to
# hand to CytoML would pull in flowWorkspace for a document this file can emit in
# a hundred lines, and would make the export unavailable wherever those packages
# are not installed, which includes the container this package tests in.

#' Subdivide polygon edges and return them in linear units
#'
#' @param poly two-column matrix of vertices on the analysis scale
#' @param transform the transform object built by [make_transform()]
#' @param marker_x,marker_y marker names, needed for a per-marker transform
#' @param n_per_edge points inserted along each edge before inversion
#' @return two-column matrix of vertices in the units the FCS file stores
#' @export
polygon_to_linear <- function(poly, transform, marker_x, marker_y,
                              n_per_edge = 24L) {
  poly <- as.matrix(poly)
  if (nrow(poly) < 3L) return(NULL)
  if (is.null(transform$inv))
    stop("this transform has no inverse; cannot express the gate in linear units",
         call. = FALSE)
  n <- nrow(poly)
  k <- max(1L, as.integer(n_per_edge))
  t_seq <- seq(0, 1, length.out = k + 1L)[-(k + 1L)]
  xs <- numeric(0); ys <- numeric(0)
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    xs <- c(xs, poly[i, 1L] + t_seq * (poly[j, 1L] - poly[i, 1L]))
    ys <- c(ys, poly[i, 2L] + t_seq * (poly[j, 2L] - poly[i, 2L]))
  }
  cbind(transform$inv(xs, marker_x), transform$inv(ys, marker_y))
}

#' Escape text for inclusion in XML
#' @param x character vector
#' @return `x` with the five predefined entities replaced
#' @keywords internal
xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&apos;", x, fixed = TRUE)
}

#' Write learned gating strategies as a Gating-ML 2.0 document
#'
#' Each strategy becomes a chain of polygon gates, level *n* declared as the
#' child of level *n-1*, which is what makes the file a gating STRATEGY rather
#' than a bag of regions: software reading it applies them in the order they were
#' learned, on the cells the previous gate kept.
#'
#' COORDINATES ARE ON THE ANALYSIS SCALE by default: one Gating-ML
#' transformation per marker (`fasinh` for arcsinh, `logicle` for logicle)
#' and each vertex divided into its units, which states the fitted region
#' exactly. `space = "linear"` writes the region in the units the FCS file
#' stores instead, edges subdivided before inversion.
#'
#' COMPENSATION. A gate drawn on data compensated with the file's own
#' `$SPILLOVER` is declared `compensation-ref="FCS"`, which tells the reader to
#' apply that same matrix; data with no matrix in the file, already unmixed by
#' the acquisition software, are declared `uncompensated`. The earlier export
#' declared every gate `uncompensated`, which told readers to gate the stored,
#' uncompensated values for gates cyRAVEN had placed on compensated ones.
#' Data compensated with a matrix from `--spillover` (`"supplied"`) carry that
#' matrix in the document as a `spectrumMatrix`, and each gate dimension it
#' compensated references its identifier, as Gating-ML 2.0 section 7 defines:
#' `"FCS"` would tell the reader to apply each file's own `$SPILLOVER`, the matrix
#' the run replaced. A dimension whose detector the matrix does not hold was
#' never compensated and is declared `uncompensated`.
#'
#' READING "FCS" WITH flowCore. A reader applying an `"FCS"` gate fetches each
#' file's `$SPILLOVER` itself. flowCore 2.18 to 2.22, every release up to
#' Bioconductor 3.22, stops there on R 4.2 or later ("the condition has length >
#' 1": `getSpilloverFromFlowFrame()` tests `class(mySpill) != "matrix"`, and a
#' matrix has two classes since R 4.0); the fix (`!is(mySpill, "matrix")`,
#' upstream commit 903b5b3) is in flowCore 2.25 (Bioconductor devel). Checked
#' with flowUtils 1.63.0: on a patched 2.18 the three cases apply and count
#' exactly the events computed directly.
#'
#' DIMENSION NAMES. Gating-ML 2.0 names a dimension of an FCS file by its `$PnN`
#' (section 4.1), also under `"FCS"` (4.2.2), so `channel_map` should map each
#' marker symbol to the files' `$PnN` ([gating_ml_markers()]); without it the
#' marker symbols are written, which no FCS file holds as a column name. Under a
#' `spectrumMatrix` a compensated dimension is named by its fluorochrome
#' instead (7.3), here the marker symbol, prefixed `Comp-` where it would equal
#' a detector name, which section 7.2 forbids; a detector carrying no marker of
#' the run is `Comp-` and its `$PnN`.
#'
#' @param polygons a polygons table from [explain_external_labels()] or
#'   [explain_unmatched_clusters()]
#' @param path destination `.xml`
#' @param transform the transform object the run used
#' @param id_col column naming the strategy each row belongs to
#' @param x_col,y_col columns holding the vertex coordinates
#' @param n_per_edge edge subdivision passed to [polygon_to_linear()], for
#'   `space = "linear"`
#' @param channel_map optional named character vector, marker symbol to the
#'   files' `$PnN`; required with `compensation_ref = "supplied"`
#' @param space "transformed" (default where the transform has a Gating-ML
#'   form) or "linear"
#' @param compensation_ref "FCS" (default), "uncompensated" or "supplied"; see
#'   above
#' @param spillover for `compensation_ref = "supplied"`, the spillover matrix
#'   the files were compensated with, the files' `$PnN` as its column names
#'   (see [gating_ml_spillover()])
#' @param M decades for the Gating-ML asinh parameters. Any positive value
#'   states the same transform; 4.5 matches the logicle default.
#' @return `path`, invisibly, or NULL when there was nothing to write
#' @export
write_gating_ml <- function(polygons, path, transform, id_col = "label",
                            x_col = "x_transformed", y_col = "y_transformed",
                            n_per_edge = 24L, channel_map = NULL,
                            space = c("transformed", "linear"),
                            compensation_ref = c("FCS", "uncompensated", "supplied"),
                            spillover = NULL, M = 4.5) {
  space <- match.arg(space)
  compensation_ref <- match.arg(compensation_ref)
  if (identical(compensation_ref, "supplied") && (is.null(spillover) || is.null(channel_map)))
    stop("compensation_ref = \"supplied\" needs the spillover matrix and channel_map",
         call. = FALSE)
  gml <- if (identical(space, "transformed")) gating_ml_scale(transform, M) else NULL
  if (identical(space, "transformed") && is.null(gml)) space <- "linear"
  if (is.null(polygons) || !nrow(polygons)) return(invisible(NULL))
  need <- c(id_col, "depth", "marker_x", "marker_y", x_col, y_col)
  miss <- setdiff(need, names(polygons))
  if (length(miss))
    stop("polygon table is missing: ", paste(miss, collapse = ", "), call. = FALSE)

  pnn <- function(m) {
    v <- if (!is.null(channel_map)) unname(channel_map[m])
    if (is.null(v) || is.na(v)) m else v
  }
  # formatC pads to its digits; an attribute value is the number alone.
  num <- function(x, digits) trimws(formatC(x, format = "g", digits = digits))
  safe_id <- function(x) gsub("[^A-Za-z0-9_.-]", "_", as.character(x))

  # Under the supplied matrix, one fluorochrome per detector: the marker the
  # detector carries, or Comp- and the detector where that name is taken by a
  # detector or there is no marker.
  spill_id <- "cyRAVEN_spillover"
  det <- if (identical(compensation_ref, "supplied")) colnames(spillover) else character(0)
  det_of <- function(m) {
    i <- match(make.names(pnn(m)), make.names(det))
    if (is.na(i)) NA_character_ else det[i]
  }
  fluor_of_det <- function(d) {
    m <- names(channel_map)[match(make.names(d), make.names(channel_map))]
    if (is.na(m) || m %in% det) paste0("Comp-", d) else m
  }
  fluor <- vapply(det, fluor_of_det, "", USE.NAMES = FALSE)
  compensated <- function(m) identical(compensation_ref, "supplied") && !is.na(det_of(m))
  comp_of <- function(m) {
    if (!identical(compensation_ref, "supplied")) return(compensation_ref)
    if (compensated(m)) spill_id else "uncompensated"
  }
  chan <- function(m)
    xml_escape(if (compensated(m)) fluor_of_det(det_of(m)) else pnn(m))

  ns <- c("gating"    = "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
          "data-type" = "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes",
          "transforms" = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations")
  L <- c('<?xml version="1.0" encoding="UTF-8"?>',
         paste0('<gating:Gating-ML xmlns:gating="', ns[["gating"]], '"'),
         paste0('                  xmlns:data-type="', ns[["data-type"]], '"'),
         paste0('                  xmlns:transforms="', ns[["transforms"]], '"'),
         '                  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">',
         paste0('  <!-- written by cyRAVEN ',
                as.character(utils::packageVersion("cyRAVEN")),
                if (identical(space, "transformed"))
                  paste0('. Coordinates are on the ',
                         xml_escape(transform$label %||% "analysis"),
                         ' scale the gates were fitted on, through the',
                         ' transformations declared below. -->')
                else paste0('. Coordinates are in the linear units stored in the',
                            ' FCS file; polygon edges were subdivided before',
                            ' inversion from the ',
                            xml_escape(transform$label %||% "analysis"),
                            ' scale. -->')))

  tid <- function(m) paste0("T_", safe_id(m))
  if (identical(space, "transformed")) {
    mk <- unique(c(polygons$marker_x, polygons$marker_y))
    for (m in mk) {
      p <- gml(m)
      if (is.null(p)) next
      L <- c(L, paste0('  <transforms:transformation transforms:id="',
                       xml_escape(tid(m)), '">'),
             paste0('    <transforms:', p$element, ' ',
                    paste0('transforms:', names(p$params), '="',
                           num(unlist(p$params), 10),
                           '"', collapse = " "), '/>'),
             '  </transforms:transformation>')
    }
  }
  if (length(fluor)) {
    fcs_dims <- function(x) paste0('      <data-type:fcs-dimension data-type:name="',
                                   xml_escape(x), '"/>')
    L <- c(L, paste0('  <transforms:spectrumMatrix transforms:id="', spill_id, '">'),
           '    <transforms:fluorochromes>', fcs_dims(fluor), '    </transforms:fluorochromes>',
           '    <transforms:detectors>', fcs_dims(det),
           '    </transforms:detectors>')
    for (i in seq_len(nrow(spillover)))
      L <- c(L, '    <transforms:spectrum>',
             paste0('      <transforms:coefficient transforms:value="',
                    num(spillover[i, ], 10), '"/>'),
             '    </transforms:spectrum>')
    L <- c(L, '  </transforms:spectrumMatrix>')
  }
  dim_open <- function(m) paste0('    <gating:dimension gating:compensation-ref="',
                                 comp_of(m), '"',
                                 if (identical(space, "transformed"))
                                   paste0(' gating:transformation-ref="',
                                          xml_escape(tid(m)), '"') else "", '>')

  n_gates <- 0L
  for (id in unique(polygons[[id_col]])) {
    d <- polygons[polygons[[id_col]] == id, , drop = FALSE]
    prev <- NULL
    for (dp in sort(unique(d$depth))) {
      v <- d[d$depth == dp, , drop = FALSE]
      v <- v[order(v$vertex), , drop = FALSE]
      mx <- v$marker_x[1L]; my <- v$marker_y[1L]
      lin <- if (identical(space, "transformed")) {
        px <- gml(mx); py <- gml(my)
        if (is.null(px) || is.null(py)) NULL
        else cbind(v[[x_col]] / px$unit, v[[y_col]] / py$unit)
      } else tryCatch(polygon_to_linear(cbind(v[[x_col]], v[[y_col]]),
                                        transform, mx, my, n_per_edge),
                      error = function(e) NULL)
      if (is.null(lin)) next
      gid <- paste0(safe_id(id), "_L", dp)
      L <- c(L, paste0('  <gating:PolygonGate gating:id="', xml_escape(gid), '"',
                       if (!is.null(prev))
                         paste0(' gating:parent_id="', xml_escape(prev), '"')
                       else "", '>'),
             dim_open(mx),
             paste0('      <data-type:fcs-dimension data-type:name="', chan(mx), '"/>'),
             '    </gating:dimension>',
             dim_open(my),
             paste0('      <data-type:fcs-dimension data-type:name="', chan(my), '"/>'),
             '    </gating:dimension>')
      for (i in seq_len(nrow(lin)))
        L <- c(L, '    <gating:vertex>',
               paste0('      <gating:coordinate data-type:value="',
                      num(lin[i, 1L], 8), '"/>'),
               paste0('      <gating:coordinate data-type:value="',
                      num(lin[i, 2L], 8), '"/>'),
               '    </gating:vertex>')
      L <- c(L, '  </gating:PolygonGate>')
      prev <- gid
      n_gates <- n_gates + 1L
    }
  }
  if (!n_gates) return(invisible(NULL))
  L <- c(L, '</gating:Gating-ML>')
  writeLines(L, path)
  log_msg("wrote ", basename(path), " (", n_gates, " polygon gate(s) in ",
          length(unique(polygons[[id_col]])), " strategy/strategies, Gating-ML 2.0)")
  invisible(path)
}

#' The Gating-ML 2.0 form of a cyRAVEN transform
#'
#' @param transform the transform object from [make_transform()].
#' @param M decades for the asinh parameters.
#' @return a function of the marker name returning `list(element, params,
#'   unit)`, where a value on the analysis scale divided by `unit` is the
#'   Gating-ML transformed value; or NULL when the transform has no Gating-ML
#'   form here.
#' @keywords internal
gating_ml_scale <- function(transform, M = 4.5) {
  if (identical(transform$method, "arcsinh")) {
    cof <- transform$cofactor_of %||% function(m) transform$params$cofactor
    return(function(m) {
      c0 <- cof(m)
      list(element = "fasinh",
           params = list(T = c0 * sinh(M * log(10)), M = M, A = 0),
           unit = M * log(10))
    })
  }
  if (identical(transform$method, "logicle")) {
    return(function(m) {
      p <- transform$params[[m]]
      if (is.null(p)) return(NULL)
      # flowCore's logicle reaches M at T; Gating-ML's reaches 1 there, whatever A
      # is. Dividing by M + A was low by A/(M + A) once A > 0.
      list(element = "logicle",
           params = list(T = p$t, W = p$w, M = p$m, A = p$a),
           unit = p$m)
    })
  }
  NULL
}

#' The same polygons as a plain table in linear units
#'
#' The format anyone can read without an XML parser, and the one to use when
#' redrawing the gate by hand at the instrument. Carries both scales side by
#' side, so a vertex can be checked against the figure it came from.
#'
#' @inheritParams write_gating_ml
#' @return a data.frame, or NULL
#' @export
polygons_linear_table <- function(polygons, transform, id_col = "label",
                                  x_col = "x_transformed", y_col = "y_transformed",
                                  n_per_edge = 24L) {
  if (is.null(polygons) || !nrow(polygons)) return(NULL)
  rows <- list()
  for (id in unique(polygons[[id_col]])) {
    d <- polygons[polygons[[id_col]] == id, , drop = FALSE]
    for (dp in sort(unique(d$depth))) {
      v <- d[d$depth == dp, , drop = FALSE]
      v <- v[order(v$vertex), , drop = FALSE]
      mx <- v$marker_x[1L]; my <- v$marker_y[1L]
      lin <- tryCatch(polygon_to_linear(cbind(v[[x_col]], v[[y_col]]),
                                        transform, mx, my, n_per_edge),
                      error = function(e) NULL)
      if (is.null(lin)) next
      rows[[length(rows) + 1L]] <- data.frame(
        strategy = as.character(id), depth = dp, marker_x = mx, marker_y = my,
        vertex = seq_len(nrow(lin)),
        x_linear = round(lin[, 1L], 4), y_linear = round(lin[, 2L], 4),
        row.names = NULL, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' The compensation a Gating-ML gate should declare for this run
#'
#' "supplied" where the files were compensated with the `--spillover` matrix,
#' which the document must then carry ([gating_ml_spillover()]); "FCS" where
#' they were compensated with their own `$SPILLOVER`; "uncompensated" where none
#' had a matrix and the stored values were gated as they are.
#' @param reads the run's reads, each carrying `compensation`.
#' @return "supplied", "FCS" or "uncompensated".
#' @keywords internal
gating_ml_compensation <- function(reads) {
  st <- vapply(reads, function(r) as.character(r$compensation %||% NA_character_)[1], "")
  if (any(grepl("^compensated.*--spillover", st))) "supplied"
  else if (any(grepl("^compensated", st))) "FCS" else "uncompensated"
}

#' Each channel's `$PnN` as the file writes it
#'
#' `channel_names` carry the `alter.names` spelling (`FITC.A`); Gating-ML names a
#' dimension by the `$PnN` keyword itself (`FITC-A`).
#' @param r one read, carrying `channel_names` and `keywords`.
#' @return character vector parallel to `channel_names`.
#' @keywords internal
read_pnn <- function(r) {
  vapply(seq_along(r$channel_names), function(i)
    as.character(r$keywords[[paste0("$P", i, "N")]] %||% r$channel_names[i])[1], "")
}

#' The supplied spillover matrix as the run applied it
#'
#' [maybe_compensate()] inverts only the detectors a file holds, so the matrix
#' a reader must apply is the supplied one restricted to the detectors the
#' run's files hold. Its columns are renamed to the files' `$PnN`, since
#' [read_spillover_file()] also accepts a header in `make.names` spelling, which
#' no file holds.
#' @param reads the run's reads, each carrying `channel_names` and `keywords`.
#' @param spillover the matrix from [read_spillover_file()], or NULL.
#' @return the restricted matrix, or NULL.
#' @keywords internal
gating_ml_spillover <- function(reads, spillover) {
  if (is.null(spillover)) return(NULL)
  pn <- unique(unlist(lapply(reads, read_pnn)))
  i <- match(make.names(colnames(spillover)), make.names(pn))
  keep <- !is.na(i)
  out <- spillover[keep, keep, drop = FALSE]
  colnames(out) <- pn[i[keep]]
  out
}

#' Marker symbol to the files' `$PnN`, over the run's reads
#' @param reads the run's reads, each carrying `marker_cols` (indices into
#'   `channel_names`), `channel_names` and `keywords`.
#' @return named character vector, marker symbol to the first `$PnN` seen for
#'   it; the `channel_map` of [write_gating_ml()].
#' @keywords internal
gating_ml_markers <- function(reads) {
  mc <- do.call(c, unname(lapply(reads, function(r)
    stats::setNames(read_pnn(r)[r$marker_cols], names(r$marker_cols)))))
  mc[!duplicated(names(mc))]
}
