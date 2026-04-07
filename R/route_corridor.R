# ---------------------------------------------------------------------------
# Internal helpers (shared by route_corridor and corridor_graph)
# ---------------------------------------------------------------------------

#' Validate and prepare a cost surface for corridor routing
#' @noRd
prepare_cost_surface <- function(cost_surface, neighbours, resolution_factor) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "Package 'terra' is required for corridor routing. ",
      "Install with install.packages('terra').",
      call. = FALSE
    )
  }

  if (!inherits(cost_surface, "SpatRaster")) {
    stop("cost_surface must be a terra SpatRaster", call. = FALSE)
  }
  if (terra::nlyr(cost_surface) != 1) {
    stop("cost_surface must have exactly 1 layer", call. = FALSE)
  }

  if (terra::crs(cost_surface) == "") {
    stop("cost_surface must have a CRS", call. = FALSE)
  }
  if (terra::is.lonlat(cost_surface)) {
    stop(
      "cost_surface must be in a projected CRS, not geographic coordinates ",
      "(e.g., EPSG:4326). Use terra::project() to reproject.",
      call. = FALSE
    )
  }

  vals_check <- terra::values(cost_surface, mat = FALSE)
  non_na <- vals_check[!is.na(vals_check)]
  if (any(!is.finite(non_na)) || any(non_na <= 0)) {
    stop(
      "All non-NA cell values in cost_surface must be finite and strictly positive (> 0).",
      call. = FALSE
    )
  }

  neighbours <- as.integer(neighbours)
  if (!neighbours %in% c(4L, 8L, 16L)) {
    stop("`neighbours` must be 4, 8, or 16", call. = FALSE)
  }

  resolution_factor <- as.integer(resolution_factor)
  if (is.na(resolution_factor) || resolution_factor < 1L) {
    stop("`resolution_factor` must be an integer >= 1", call. = FALSE)
  }
  if (resolution_factor > 1L) {
    cost_surface <- terra::aggregate(
      cost_surface,
      fact = resolution_factor,
      fun = "mean",
      na.rm = FALSE
    )
  }

  cost_surface
}

#' Normalize a point to numeric c(x, y) with CRS reconciliation
#' @noRd
normalize_corridor_point <- function(pt, name, raster_crs) {
  if (inherits(pt, "sf")) {
    if (nrow(pt) != 1) {
      stop(sprintf("`%s` must be a single point feature", name), call. = FALSE)
    }
    geom_type <- as.character(sf::st_geometry_type(pt, by_geometry = TRUE))
    if (geom_type != "POINT") {
      stop(sprintf("`%s` must be a POINT geometry, got %s", name, geom_type), call. = FALSE)
    }
    if (is.na(sf::st_crs(pt))) {
      stop(sprintf(
        "`%s` has no CRS. Supply coordinates in the CRS of cost_surface or set the CRS on the sf object.",
        name
      ), call. = FALSE)
    }
    if (sf::st_crs(pt) != raster_crs) {
      pt <- sf::st_transform(pt, raster_crs)
    }
    coords <- sf::st_coordinates(pt)
    c(coords[1, "X"], coords[1, "Y"])

  } else if (inherits(pt, "sfc")) {
    if (length(pt) != 1) {
      stop(sprintf("`%s` must be a single point geometry", name), call. = FALSE)
    }
    geom_type <- as.character(sf::st_geometry_type(pt, by_geometry = TRUE))
    if (geom_type != "POINT") {
      stop(sprintf("`%s` must be a POINT geometry, got %s", name, geom_type), call. = FALSE)
    }
    if (is.na(sf::st_crs(pt))) {
      stop(sprintf(
        "`%s` has no CRS. Supply coordinates in the CRS of cost_surface or set the CRS on the sf object.",
        name
      ), call. = FALSE)
    }
    if (sf::st_crs(pt) != raster_crs) {
      pt <- sf::st_transform(pt, raster_crs)
    }
    coords <- sf::st_coordinates(pt)
    c(coords[1, "X"], coords[1, "Y"])

  } else if (is.numeric(pt) && length(pt) == 2L) {
    if (any(!is.finite(pt))) {
      stop(sprintf("`%s` coordinates must be finite", name), call. = FALSE)
    }
    pt

  } else {
    stop(
      sprintf("`%s` must be an sf/sfc POINT or a numeric c(x, y) vector", name),
      call. = FALSE
    )
  }
}

#' Validate cell indices for origin/destination
#' @noRd
validate_corridor_cells <- function(cost_surface, origin_xy, dest_xy) {
  origin_cell <- terra::cellFromXY(cost_surface, matrix(origin_xy, ncol = 2))
  dest_cell   <- terra::cellFromXY(cost_surface, matrix(dest_xy, ncol = 2))

  if (is.na(origin_cell)) {
    stop("origin falls outside the cost_surface extent", call. = FALSE)
  }
  if (is.na(dest_cell)) {
    stop("destination falls outside the cost_surface extent", call. = FALSE)
  }

  origin_val <- terra::extract(cost_surface, origin_cell)[1, 1]
  dest_val   <- terra::extract(cost_surface, dest_cell)[1, 1]

  if (is.na(origin_val)) {
    stop(
      "Origin falls on an NA cell (impassable). ",
      "Ensure the point falls on a valid cell, or adjust the cost surface.",
      call. = FALSE
    )
  }
  if (is.na(dest_val)) {
    stop(
      "Destination falls on an NA cell (impassable). ",
      "Ensure the point falls on a valid cell, or adjust the cost surface.",
      call. = FALSE
    )
  }

  list(origin_cell = origin_cell, dest_cell = dest_cell)
}

#' Build sf LINESTRING from Rust corridor result
#' @noRd
corridor_result_to_sf <- function(result, cost_surface, origin_xy, dest_xy,
                                  method, neighbours) {
  raster_crs <- sf::st_crs(terra::crs(cost_surface))

  # Convert 0-based cell indices to 1-based
  cell_indices_1 <- result$path_cells + 1L

  # Get cell center coordinates
  path_xy <- terra::xyFromCell(cost_surface, cell_indices_1)

  # Get cell centers for origin and dest (for metadata)
  origin_cell <- terra::cellFromXY(cost_surface, matrix(origin_xy, ncol = 2))
  dest_cell   <- terra::cellFromXY(cost_surface, matrix(dest_xy, ncol = 2))
  origin_cell_center <- terra::xyFromCell(cost_surface, as.integer(origin_cell))
  dest_cell_center   <- terra::xyFromCell(cost_surface, as.integer(dest_cell))

  # Compute straight-line distance
  straight_line_dist <- sqrt(sum((origin_xy - dest_xy)^2))
  is_degenerate <- straight_line_dist == 0

  if (is_degenerate) {
    coords <- rbind(origin_xy, origin_xy)
    colnames(coords) <- c("x", "y")
  } else {
    coords <- rbind(
      matrix(origin_xy, ncol = 2, dimnames = list(NULL, c("x", "y"))),
      path_xy,
      matrix(dest_xy, ncol = 2, dimnames = list(NULL, c("x", "y")))
    )
  }

  # Build sf LINESTRING
  line <- sf::st_linestring(coords)
  sfc  <- sf::st_sfc(line, crs = raster_crs)
  output <- sf::st_sf(geometry = sfc)

  # Compute metrics
  path_dist <- if (is_degenerate) 0 else as.numeric(sf::st_length(output))
  sinuosity <- if (straight_line_dist > 0) path_dist / straight_line_dist else NA_real_

  n_rows <- terra::nrow(cost_surface)
  n_cols <- terra::ncol(cost_surface)
  res    <- terra::res(cost_surface)

  output$total_cost         <- result$total_cost
  output$n_cells            <- length(result$path_cells)
  output$straight_line_dist <- straight_line_dist
  output$path_dist          <- path_dist
  output$sinuosity          <- sinuosity

  attr(output, "spopt") <- list(
    total_cost       = result$total_cost,
    n_cells          = length(result$path_cells),
    method           = method,
    neighbours       = neighbours,
    n_cells_surface  = as.integer(n_rows) * as.integer(n_cols),
    n_edges_graph    = result$n_edges,
    solve_time       = result$solve_time_ms / 1000,
    graph_build_time = result$graph_build_time_ms / 1000,
    cell_indices     = cell_indices_1,
    origin_cell_center = as.numeric(origin_cell_center),
    dest_cell_center   = as.numeric(dest_cell_center),
    raster_dims      = c(n_rows, n_cols),
    cell_size        = c(res[1], res[2])
  )

  class(output) <- c("spopt_corridor", class(output))
  output
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Least-Cost Corridor Routing
#'
#' Find the minimum-cost path between an origin and destination across a
#' raster friction surface. The cost surface must be in a projected CRS
#' with equal-area cells (not EPSG:4326). Cell values represent traversal
#' friction (cost per unit distance). Higher values = more expensive.
#' NA cells are impassable.
#'
#' @param cost_surface A terra SpatRaster (single band) or an
#'   \code{spopt_corridor_graph} object created by \code{\link{corridor_graph}}.
#'   When a graph object is supplied, \code{neighbours} and
#'   \code{resolution_factor} are fixed at graph build time and cannot be
#'   overridden.
#' @param origin sf/sfc POINT (single feature) or numeric vector c(x, y)
#'   in the CRS of cost_surface. Must fall on a non-NA cell.
#' @param destination sf/sfc POINT (single feature) or numeric vector c(x, y).
#'   Must fall on a non-NA cell.
#' @param neighbours Integer. Cell connectivity: 4, 8 (default), or 16.
#'   4 = cardinal only. 8 = cardinal + diagonal. 16 = adds knight's-move.
#'   Cannot be overridden when \code{cost_surface} is a corridor graph
#'   (fixed at \code{\link{corridor_graph}} build time).
#' @param method Character. Routing algorithm:
#'   - "dijkstra" (default): standard Dijkstra's shortest path.
#'   - "bidirectional": bidirectional Dijkstra, ~2x faster.
#'   - "astar": A* with Euclidean heuristic, fastest for distant pairs.
#' @param resolution_factor Integer, default 1L. Values > 1 aggregate the
#'   cost surface before routing (e.g., 2L halves resolution). Cannot be
#'   overridden when \code{cost_surface} is a corridor graph (fixed at
#'   \code{\link{corridor_graph}} build time).
#'
#' @return An sf LINESTRING object (single feature) with columns:
#'   \itemize{
#'     \item \code{total_cost}: accumulated traversal cost (friction * distance units)
#'     \item \code{n_cells}: number of cells in the path
#'     \item \code{straight_line_dist}: Euclidean distance origin to destination (CRS units)
#'     \item \code{path_dist}: actual path length (CRS units)
#'     \item \code{sinuosity}: path_dist / straight_line_dist (NA if origin == destination)
#'   }
#'
#'   The returned linestring starts at the user-supplied origin coordinates,
#'   passes through the cell centers along the optimal path, and ends at
#'   the user-supplied destination coordinates.
#'
#'   The "spopt" attribute contains metadata including \code{total_cost},
#'   \code{n_cells}, \code{method}, \code{neighbours}, \code{solve_time},
#'   \code{graph_build_time}, \code{cell_indices}, and grid dimensions.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' # Build a simple friction surface (projected CRS)
#' r <- rast(nrows = 500, ncols = 500, xmin = 0, xmax = 500000,
#'           ymin = 0, ymax = 500000, crs = "EPSG:32614")
#' values(r) <- runif(ncell(r), 0.5, 2.0)
#'
#' origin <- st_sfc(st_point(c(50000, 50000)), crs = 32614)
#' dest   <- st_sfc(st_point(c(450000, 450000)), crs = 32614)
#'
#' path <- route_corridor(r, origin, dest)
#' plot(r)
#' plot(st_geometry(path), add = TRUE, col = "red", lwd = 2)
#'
#' # Graph caching for multiple OD pairs
#' g <- corridor_graph(r, neighbours = 8L)
#' path1 <- route_corridor(g, origin, dest, method = "astar")
#' }
#'
#' @export
route_corridor <- function(cost_surface,
                           origin,
                           destination,
                           neighbours = 8L,
                           method = c("dijkstra", "bidirectional", "astar"),
                           resolution_factor = 1L) {

  method <- match.arg(method)

  # --- dispatch: cached graph vs raster ---
  if (inherits(cost_surface, "spopt_corridor_graph")) {
    return(.route_corridor_cached(
      cost_surface, origin, destination, method,
      call_args = match.call()
    ))
  }

  # --- raster path ---
  cost_surface <- prepare_cost_surface(cost_surface, neighbours, resolution_factor)

  neighbours <- as.integer(neighbours)
  raster_crs <- sf::st_crs(terra::crs(cost_surface))

  origin_xy <- normalize_corridor_point(origin, "origin", raster_crs)
  dest_xy   <- normalize_corridor_point(destination, "destination", raster_crs)

  cells <- validate_corridor_cells(cost_surface, origin_xy, dest_xy)

  # Extract raster data
  values <- terra::values(cost_surface, mat = FALSE)
  values[is.na(values)] <- NaN
  n_rows <- terra::nrow(cost_surface)
  n_cols <- terra::ncol(cost_surface)
  res    <- terra::res(cost_surface)

  origin_cell_0 <- as.integer(cells$origin_cell - 1L)
  dest_cell_0   <- as.integer(cells$dest_cell - 1L)

  result <- rust_corridor(
    values, as.integer(n_rows), as.integer(n_cols),
    res[1], res[2],
    origin_cell_0, dest_cell_0,
    neighbours, method
  )

  corridor_result_to_sf(result, cost_surface, origin_xy, dest_xy, method, neighbours)
}

# ---------------------------------------------------------------------------
# Cached graph path
# ---------------------------------------------------------------------------

#' @noRd
.route_corridor_cached <- function(graph, origin, destination, method, call_args) {
  # Reject conflicting overrides
  explicit_args <- names(call_args)
  if ("neighbours" %in% explicit_args) {
    stop(
      "`neighbours` cannot be overridden on a cached graph. ",
      "This is fixed at corridor_graph() build time.",
      call. = FALSE
    )
  }
  if ("resolution_factor" %in% explicit_args) {
    stop(
      "`resolution_factor` cannot be overridden on a cached graph. ",
      "This is fixed at corridor_graph() build time.",
      call. = FALSE
    )
  }

  # Pointer validity
  ptr_valid <- tryCatch(
    { rust_corridor_graph_info(graph$ptr); TRUE },
    error = function(e) FALSE
  )
  if (!ptr_valid) {
    stop(
      "Graph pointer is invalid (possibly deserialized). ",
      "Rebuild with corridor_graph().",
      call. = FALSE
    )
  }

  meta <- attr(graph, "spopt")
  raster_crs <- meta$crs
  neighbours <- meta$neighbours

  origin_xy <- normalize_corridor_point(origin, "origin", raster_crs)
  dest_xy   <- normalize_corridor_point(destination, "destination", raster_crs)

  cells <- validate_corridor_cells(graph$cost_surface, origin_xy, dest_xy)

  origin_cell_0 <- as.integer(cells$origin_cell - 1L)
  dest_cell_0   <- as.integer(cells$dest_cell - 1L)

  result <- rust_corridor_solve_cached(
    graph$ptr, origin_cell_0, dest_cell_0, method
  )

  corridor_result_to_sf(
    result, graph$cost_surface, origin_xy, dest_xy, method, neighbours
  )
}

# ---------------------------------------------------------------------------
# Graph constructor
# ---------------------------------------------------------------------------

#' Build a Corridor Graph for Cached Routing
#'
#' Pre-build the routing graph from a cost surface so that multiple
#' \code{\link{route_corridor}} calls can skip graph construction.
#' The returned object is a snapshot of the cost surface at build time;
#' subsequent edits to the raster do not affect the graph.
#'
#' The object retains both the CSR graph (in Rust) and an independent copy
#' of the raster (for coordinate mapping). The printed \code{graph_storage}
#' reflects only the CSR arrays, not the raster copy.
#'
#' @param cost_surface A terra SpatRaster (single band). Same requirements
#'   as \code{\link{route_corridor}}.
#' @param neighbours Integer. Cell connectivity: 4, 8 (default), or 16.
#'   Fixed at build time.
#' @param resolution_factor Integer, default 1L. If > 1, the surface is
#'   aggregated before graph construction. Fixed at build time.
#'
#' @return An opaque \code{spopt_corridor_graph} object for use with
#'   \code{\link{route_corridor}}.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(nrows = 500, ncols = 500, xmin = 0, xmax = 500000,
#'           ymin = 0, ymax = 500000, crs = "EPSG:32614")
#' values(r) <- runif(ncell(r), 0.5, 2.0)
#'
#' g <- corridor_graph(r, neighbours = 8L)
#' print(g)
#' path <- route_corridor(g, c(50000, 50000), c(450000, 450000))
#' }
#'
#' @export
corridor_graph <- function(cost_surface, neighbours = 8L, resolution_factor = 1L) {
  cost_surface <- prepare_cost_surface(cost_surface, neighbours, resolution_factor)

  neighbours <- as.integer(neighbours)

  values <- terra::values(cost_surface, mat = FALSE)
  values[is.na(values)] <- NaN
  n_rows <- terra::nrow(cost_surface)
  n_cols <- terra::ncol(cost_surface)
  res    <- terra::res(cost_surface)

  ptr <- rust_corridor_build_graph(
    values, as.integer(n_rows), as.integer(n_cols),
    res[1], res[2], neighbours
  )

  info <- rust_corridor_graph_info(ptr)

  surface_snapshot <- terra::deepcopy(cost_surface)

  structure(
    list(
      ptr = ptr,
      cost_surface = surface_snapshot
    ),
    spopt = list(
      n_rows          = info$n_rows,
      n_cols          = info$n_cols,
      cell_size       = c(info$cell_width, info$cell_height),
      crs             = sf::st_crs(terra::crs(cost_surface)),
      neighbours      = info$neighbours,
      n_edges         = info$n_edges,
      graph_build_time = info$build_time_ms / 1000,
      n_cells_surface = info$n_rows * info$n_cols,
      graph_storage_mb = round(info$memory_bytes / 1e6, 1)
    ),
    class = "spopt_corridor_graph"
  )
}
