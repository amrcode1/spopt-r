# Solve P-Median facility location problem

Solve P-Median facility location problem

## Usage

``` r
rust_p_median(cost_matrix, weights, n_facilities, fixed_facilities)
```

## Arguments

- cost_matrix:

  Cost/distance matrix (demand x facilities)

- weights:

  Demand weights

- n_facilities:

  Number of facilities to locate (p)

- fixed_facilities:

  Optional indices of pre-selected facilities (1-based)

## Value

List with selected facilities and assignments
