# Solve MCLP (Maximum Coverage Location Problem)

Solve MCLP (Maximum Coverage Location Problem)

## Usage

``` r
rust_mclp(cost_matrix, weights, service_radius, n_facilities, fixed_facilities)
```

## Arguments

- cost_matrix:

  Cost/distance matrix (demand x facilities)

- weights:

  Demand weights

- service_radius:

  Maximum service distance

- n_facilities:

  Number of facilities to locate

## Value

List with selected facilities and coverage
