# Solve Capacitated Vehicle Routing Problem (CVRP)

Find minimum-cost routes for multiple vehicles, each with a capacity
limit, to serve all customers from a depot. Uses Clarke-Wright savings
heuristic with 2-opt, relocate, and swap improvement.

## Usage

``` r
rust_vrp(
  cost_matrix,
  depot,
  demands,
  capacity,
  max_vehicles,
  method,
  service_times,
  max_route_time,
  balance_time,
  earliest,
  latest
)
```

## Arguments

- cost_matrix:

  Square cost/distance matrix (n x n)

- depot:

  Depot index (0-based)

- demands:

  Demand at each location (depot demand should be 0)

- capacity:

  Vehicle capacity

- max_vehicles:

  Maximum number of vehicles (NULL for unlimited)

- method:

  Algorithm: "savings" (construction only) or "2-opt" (with improvement)

- service_times:

  Optional service time at each stop (NULL for zero)

- max_route_time:

  Optional maximum total time per route (NULL for unlimited)

- balance_time:

  Whether to run route-time balancing phase

- earliest:

  Optional earliest arrival times at each stop

- latest:

  Optional latest arrival times at each stop

## Value

List with vehicle assignments, visit orders, costs, and route details
