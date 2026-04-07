# Least-cost corridor routing on a raster grid

Least-cost corridor routing on a raster grid

## Usage

``` r
rust_corridor(
  values,
  n_rows,
  n_cols,
  cell_width,
  cell_height,
  origin_cell,
  dest_cell,
  neighbours,
  method
)
```

## Arguments

- values:

  Flattened cell values (row-major), NaN = impassable

- n_rows:

  Number of rows in the raster

- n_cols:

  Number of columns in the raster

- cell_width:

  Cell width in CRS units

- cell_height:

  Cell height in CRS units

- origin_cell:

  0-based cell index of origin

- dest_cell:

  0-based cell index of destination

- neighbours:

  Cell connectivity: 4, 8, or 16

- method:

  Routing algorithm: "dijkstra", "bidirectional", or "astar"

## Value

List with path_cells, total_cost, solve_time_ms, graph_build_time_ms,
n_edges
