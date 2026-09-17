# Subdivide polygon edges and return them in linear units

Subdivide polygon edges and return them in linear units

## Usage

``` r
polygon_to_linear(poly, transform, marker_x, marker_y, n_per_edge = 24L)
```

## Arguments

- poly:

  two-column matrix of vertices on the analysis scale

- transform:

  the transform object built by
  [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md)

- marker_x, marker_y:

  marker names, needed for a per-marker transform

- n_per_edge:

  points inserted along each edge before inversion

## Value

two-column matrix of vertices in the units the FCS file stores
