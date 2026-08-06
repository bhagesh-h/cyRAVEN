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
# WHY THE VERTICES ARE INVERSE-TRANSFORMED, AND WHY THAT NEEDS MORE OF THEM.
# The polygons are fitted on the analysis scale, where their edges are straight.
# A gate file consumed by other software has to describe the region in the units
# the FCS file stores, and the transform is non-linear, so a straight edge on the
# analysis scale is a CURVE in linear units. Inverse-transforming only the
# corners would emit a polygon with straight edges between them, which is a
# different region -- narrower or wider than the gate that was actually fitted
# and validated, by an amount that grows with how far the edge spans.
#
# So each edge is subdivided before inversion. With the default subdivision the
# emitted boundary tracks the fitted one to well within the width of a histogram
# bin, and the cost is vertices in a text file.
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
#' COORDINATES ARE LINEAR AND UNCOMPENSATED-REFERENCED. The gates are emitted in
#' the units the FCS file stores, so no transformation element is needed and
#' there is nothing for the reader to get wrong. `compensation-ref` is declared
#' as `uncompensated` because cyRAVEN applies the acquisition spillover matrix
#' from the file's own keyword block before gating; a reader that applies its own
#' compensation on top would be compensating twice.
#'
#' DIMENSION NAMES ARE MARKER SYMBOLS. cyRAVEN resolves channels through `$PnS`,
#' so a gate names CD3 rather than the detector. Software keyed on detector names
#' needs `channel_map` to translate.
#'
#' @param polygons a polygons table from [explain_external_labels()] or
#'   [explain_unmatched_clusters()]
#' @param path destination `.xml`
#' @param transform the transform object the run used
#' @param id_col column naming the strategy each row belongs to
#' @param x_col,y_col columns holding the vertex coordinates
#' @param n_per_edge edge subdivision passed to [polygon_to_linear()]
#' @param channel_map optional named character vector, marker symbol to channel
#' @return `path`, invisibly, or NULL when there was nothing to write
#' @export
write_gating_ml <- function(polygons, path, transform, id_col = "label",
                            x_col = "x_transformed", y_col = "y_transformed",
                            n_per_edge = 24L, channel_map = NULL) {
  if (is.null(polygons) || !nrow(polygons)) return(invisible(NULL))
  need <- c(id_col, "depth", "marker_x", "marker_y", x_col, y_col)
  miss <- setdiff(need, names(polygons))
  if (length(miss))
    stop("polygon table is missing: ", paste(miss, collapse = ", "), call. = FALSE)

  chan <- function(m) xml_escape(unname(channel_map[m]) %||% m)
  safe_id <- function(x) gsub("[^A-Za-z0-9_.-]", "_", as.character(x))

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
                '. Coordinates are in the linear units stored in the FCS file;',
                ' polygon edges were subdivided before inversion from the ',
                xml_escape(transform$label %||% "analysis"), ' scale. -->'))

  n_gates <- 0L
  for (id in unique(polygons[[id_col]])) {
    d <- polygons[polygons[[id_col]] == id, , drop = FALSE]
    prev <- NULL
    for (dp in sort(unique(d$depth))) {
      v <- d[d$depth == dp, , drop = FALSE]
      v <- v[order(v$vertex), , drop = FALSE]
      mx <- v$marker_x[1L]; my <- v$marker_y[1L]
      lin <- tryCatch(polygon_to_linear(cbind(v[[x_col]], v[[y_col]]),
                                        transform, mx, my, n_per_edge),
                      error = function(e) NULL)
      if (is.null(lin)) next
      gid <- paste0(safe_id(id), "_L", dp)
      L <- c(L, paste0('  <gating:PolygonGate gating:id="', xml_escape(gid), '"',
                       if (!is.null(prev))
                         paste0(' gating:parent_id="', xml_escape(prev), '"')
                       else "", '>'),
             '    <gating:dimension gating:compensation-ref="uncompensated">',
             paste0('      <data-type:fcs-dimension data-type:name="', chan(mx), '"/>'),
             '    </gating:dimension>',
             '    <gating:dimension gating:compensation-ref="uncompensated">',
             paste0('      <data-type:fcs-dimension data-type:name="', chan(my), '"/>'),
             '    </gating:dimension>')
      for (i in seq_len(nrow(lin)))
        L <- c(L, '    <gating:vertex>',
               paste0('      <gating:coordinate data-type:value="',
                      formatC(lin[i, 1L], format = "g", digits = 8), '"/>'),
               paste0('      <gating:coordinate data-type:value="',
                      formatC(lin[i, 2L], format = "g", digits = 8), '"/>'),
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
