# Solve Max-P regionalization problem

Maximize the number of regions such that each region satisfies a minimum
threshold constraint on a spatial extensive attribute.

## Usage

``` r
rust_max_p(
  attrs,
  threshold_var,
  threshold,
  adj_i,
  adj_j,
  n_iterations,
  n_sa_iterations,
  cooling_rate,
  tabu_length,
  seed,
  homogeneous,
  compact,
  compact_weight,
  compact_metric,
  centroids_x,
  centroids_y,
  areas,
  moments
)
```

## Arguments

- attrs:

  Attribute matrix (n x p) for computing within-region dissimilarity

- threshold_var:

  Values of the threshold variable (e.g., population)

- threshold:

  Minimum sum required per region

- adj_i:

  Row indices of adjacency (0-based)

- adj_j:

  Column indices of adjacency (0-based)

- n_iterations:

  Number of construction phase iterations

- n_sa_iterations:

  Number of simulated annealing iterations

- cooling_rate:

  SA cooling rate (e.g., 0.99)

- tabu_length:

  Tabu list length for SA

- seed:

  Random seed

- compact:

  Whether to optimize for compactness

- compact_weight:

  Weight for compactness vs dissimilarity (0-1)

- compact_metric:

  Compactness metric to use, either "centroid" or "nmi" (normalized
  moment of inertia)

- centroids_x:

  X coordinates of unit centroids (for compactness)

- centroids_y:

  Y coordinates of unit centroids (for compactness)

- areas:

  Areas of units (for compactness)

- moments:

  Second areal moments of units (for compactness)

## Value

List with labels (1-based), n_regions, objective, and compactness
