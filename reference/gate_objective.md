# Log-likelihood and gradient of a soft convex gate

WHAT: class-weighted binary cross-entropy between the relaxed gate
membership `prod_j sigmoid(s * (w_j . x + b_j))` and the label, plus a
tightness penalty.

## Usage

``` r
gate_objective(
  theta,
  X,
  y,
  wt,
  fixed_W,
  centroid,
  s = 40,
  lambda = 0,
  eps = 0.001
)
```

## Arguments

- theta:

  packed parameters: one angle per free plane, then every offset

- X:

  n x 2 matrix of marker values, scaled to the unit square

- y:

  0/1 target indicator

- wt:

  per-cell weight

- fixed_W:

  2 x 4 matrix of the frozen (PC-aligned) unit normals

- centroid:

  length-2 centroid of the target cells

- s:

  sigmoid steepness

- lambda:

  tightness penalty strength

- eps:

  membership clamp

## Value

objective value with a `gradient` attribute

## Details

WHY THE PENALTY EXISTS AND WHAT IT DOES. Cross-entropy alone has no
opinion about where a plane sits once it has separated the classes, so
planes drift outward and the gate ends up loose – correct on the cells
it was shown, and admitting anything that happens to lie beyond them.
The penalty is the mean distance from each plane to the centroid of the
target cells, which pulls every boundary inward. It is a shrinkage prior
on gate size, and its strength `lambda` is not guessed:
learn_convex_gate() searches it against held-out F1, so the data chooses
how tight the gate should be.

WHY THE MEMBERSHIP IS CLAMPED. Deep inside the gate p is numerically 1,
and the cross-entropy gradient for a misplaced non-target there is
proportional to p/(1-p), which overflows. Rescaling p onto
`[eps, 1-eps]` bounds that ratio by 1/eps while staying smooth, so a
single badly-placed cell cannot dominate the step.
