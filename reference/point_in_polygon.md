# Which points fall inside a polygon

WHY A WINDING-NUMBER TEST AND NOT A RAY CAST: both are short, but the
winding rule gives the same answer for a vertex-on-boundary point
regardless of which direction the ray was cast, so a cell sitting
exactly on a gate edge does not change classification when the polygon
is rotated. Vectorised over points.

## Usage

``` r
point_in_polygon(xy, poly)
```

## Arguments

- xy:

  n x 2 matrix of points

- poly:

  m x 2 matrix of polygon vertices in order

## Value

logical vector
