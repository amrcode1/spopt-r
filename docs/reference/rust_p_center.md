# Solve P-Center facility location problem

Solve P-Center facility location problem

## Usage

``` r
rust_p_center(cost_matrix, n_facilities, method, fixed_facilities)
```

## Arguments

- cost_matrix:

  Cost/distance matrix (demand x facilities)

- n_facilities:

  Number of facilities to locate

- method:

  Algorithm method: "binary_search" (default) or "mip"

## Value

List with selected facilities, assignments, and max distance
