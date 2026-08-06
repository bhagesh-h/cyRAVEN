# Vertices of the polygon defined by an intersection of half-planes

WHY IT ENUMERATES PAIRS: with k planes there are at most choose(k, 2)
candidate corners, 66 at the default k = 12, and each is one 2x2 solve.
A general half-space enumeration algorithm would be faster
asymptotically and slower here, and would be one more thing that can be
subtly wrong.

## Usage

``` r
halfplane_polygon(W, b, tol = 1e-09)
```

## Arguments

- W:

  2 x k matrix of normals

- b:

  length-k offsets

- tol:

  feasibility tolerance

## Value

matrix of vertices in draw order, or NULL

## Details

Returns NULL when the region is empty or degenerate, which is the honest
answer: an optimiser that produced no feasible region has not produced a
gate.
