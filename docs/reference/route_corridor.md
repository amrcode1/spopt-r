# Least-Cost Corridor Routing

Find the minimum-cost path between an origin and destination across a
raster friction surface. The cost surface must be in a projected CRS
with equal-area cells (not EPSG:4326). Cell values represent traversal
friction (cost per unit distance). Higher values = more expensive. NA
cells are impassable.

## Usage

``` r
route_corridor(
  cost_surface,
  origin,
  destination,
  neighbours = 8L,
  method = c("dijkstra", "bidirectional", "astar"),
  resolution_factor = 1L
)
```

## Arguments

- cost_surface:

  A terra SpatRaster (single band) or an `spopt_corridor_graph` object
  created by
  [`corridor_graph`](https://walker-data.com/spopt-r/reference/corridor_graph.md).
  When a graph object is supplied, `neighbours` and `resolution_factor`
  are fixed at graph build time and cannot be overridden.

- origin:

  sf/sfc POINT (single feature) or numeric vector c(x, y) in the CRS of
  cost_surface. Must fall on a non-NA cell.

- destination:

  sf/sfc POINT (single feature) or numeric vector c(x, y). Must fall on
  a non-NA cell.

- neighbours:

  Integer. Cell connectivity: 4, 8 (default), or 16. 4 = cardinal only.
  8 = cardinal + diagonal. 16 = adds knight's-move. Cannot be overridden
  when `cost_surface` is a corridor graph (fixed at
  [`corridor_graph`](https://walker-data.com/spopt-r/reference/corridor_graph.md)
  build time).

- method:

  Character. Routing algorithm:

  - "dijkstra" (default): standard Dijkstra's shortest path.

  - "bidirectional": bidirectional Dijkstra, ~2x faster.

  - "astar": A\* with Euclidean heuristic, fastest for distant pairs.

- resolution_factor:

  Integer, default 1L. Values \> 1 aggregate the cost surface before
  routing (e.g., 2L halves resolution). Cannot be overridden when
  `cost_surface` is a corridor graph (fixed at
  [`corridor_graph`](https://walker-data.com/spopt-r/reference/corridor_graph.md)
  build time).

## Value

An sf LINESTRING object (single feature) with columns:

- `total_cost`: accumulated traversal cost (friction \* distance units)

- `n_cells`: number of cells in the path

- `straight_line_dist`: Euclidean distance origin to destination (CRS
  units)

- `path_dist`: actual path length (CRS units)

- `sinuosity`: path_dist / straight_line_dist (NA if origin ==
  destination)

The returned linestring starts at the user-supplied origin coordinates,
passes through the cell centers along the optimal path, and ends at the
user-supplied destination coordinates.

The "spopt" attribute contains metadata including `total_cost`,
`n_cells`, `method`, `neighbours`, `solve_time`, `graph_build_time`,
`cell_indices`, and grid dimensions.

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
library(sf)

# Build a simple friction surface (projected CRS)
r <- rast(nrows = 500, ncols = 500, xmin = 0, xmax = 500000,
          ymin = 0, ymax = 500000, crs = "EPSG:32614")
values(r) <- runif(ncell(r), 0.5, 2.0)

origin <- st_sfc(st_point(c(50000, 50000)), crs = 32614)
dest   <- st_sfc(st_point(c(450000, 450000)), crs = 32614)

path <- route_corridor(r, origin, dest)
plot(r)
plot(st_geometry(path), add = TRUE, col = "red", lwd = 2)

# Graph caching for multiple OD pairs
g <- corridor_graph(r, neighbours = 8L)
path1 <- route_corridor(g, origin, dest, method = "astar")
} # }
```
