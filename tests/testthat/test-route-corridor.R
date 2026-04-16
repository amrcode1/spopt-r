skip_if_not_installed("terra")
skip_if_not_installed("sf")

library(terra)
library(sf)

# make_test_raster() is defined in helper-raster.R

# ---------------------------------------------------------------------------
# 1. Uniform surface -- path approximates straight line
# ---------------------------------------------------------------------------
test_that("uniform surface produces near-straight-line path", {
  r <- make_test_raster(50, 50, vals = 1)
  origin <- c(500, 500)
  dest   <- c(4500, 4500)

  path <- route_corridor(r, origin, dest)

  expect_s3_class(path, "spopt_corridor")
  expect_s3_class(path, "sf")
  expect_true("total_cost" %in% names(path))
  expect_true("n_cells" %in% names(path))
  expect_true("sinuosity" %in% names(path))

  # On a uniform grid, sinuosity should be close to 1
  expect_lt(path$sinuosity, 1.5)
})

# ---------------------------------------------------------------------------
# 2. Barrier avoidance (NA band)
# ---------------------------------------------------------------------------
test_that("path routes around NA barrier", {
  r <- make_test_raster(50, 50, vals = 1)
  # Place vertical NA barrier in the middle (cols 24-26), but leave a gap at rows 1-3
  cells <- cellFromRowCol(r, rep(4:50, each = 3), rep(24:26, times = 47))
  r[cells] <- NA

  origin <- c(500, 2500)
  dest   <- c(4500, 2500)

  path <- route_corridor(r, origin, dest)

  # Cost should exceed straight-line distance since path must detour
  expect_gt(path$total_cost, path$straight_line_dist)
  expect_gt(path$sinuosity, 1.0)
})

# ---------------------------------------------------------------------------
# 3. High-cost avoidance
# ---------------------------------------------------------------------------
test_that("path avoids high-cost band", {
  r <- make_test_raster(50, 50, vals = 1)
  # High-cost vertical band in the middle
  cells <- cellFromRowCol(r, rep(1:50, each = 3), rep(24:26, 50))
  r[cells] <- 100

  origin <- c(500, 2500)
  dest   <- c(4500, 2500)

  path <- route_corridor(r, origin, dest)
  # The path should detour, so sinuosity > 1

  expect_gt(path$sinuosity, 1.0)
})

# ---------------------------------------------------------------------------
# 4. Preference corridor (low-friction band)
# ---------------------------------------------------------------------------
test_that("path follows low-friction corridor", {
  r <- make_test_raster(50, 50, vals = 10)
  # Low-friction corridor along row 25
  cells <- cellFromRowCol(r, rep(25, 50), 1:50)
  r[cells] <- 0.3

  origin <- c(100, 2500)
  dest   <- c(4900, 2500)

  path <- route_corridor(r, origin, dest)

  # Path should be cheaper than going through the expensive zone
  meta <- attr(path, "spopt")
  expect_true(meta$total_cost < 10 * path$straight_line_dist)
})

# ---------------------------------------------------------------------------
# 5. CRS validation -- geographic CRS rejected
# ---------------------------------------------------------------------------
test_that("geographic CRS is rejected", {
  r <- rast(
    nrows = 10, ncols = 10,
    xmin = -90, xmax = -80, ymin = 30, ymax = 40,
    crs = "EPSG:4326"
  )
  values(r) <- 1

  expect_error(
    route_corridor(r, c(-85, 35), c(-82, 37)),
    "projected CRS"
  )
})

# ---------------------------------------------------------------------------
# 6. Cell value validation -- zero values
# ---------------------------------------------------------------------------
test_that("zero cell values are rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  r[5] <- 0

  expect_error(
    route_corridor(r, c(100, 100), c(900, 900)),
    "positive"
  )
})

# ---------------------------------------------------------------------------
# 7. Origin on NA cell
# ---------------------------------------------------------------------------
test_that("origin on NA cell is rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  r[cellFromXY(r, matrix(c(500, 500), ncol = 2))] <- NA

  expect_error(
    route_corridor(r, c(500, 500), c(900, 900)),
    "impassable|NA"
  )
})

# ---------------------------------------------------------------------------
# 8. Method equivalence
# ---------------------------------------------------------------------------
test_that("all three methods produce the same total cost", {
  r <- make_test_raster(30, 30, vals = 1)
  # Add some variation
  set.seed(42)
  values(r) <- runif(ncell(r), 0.5, 2.0)

  origin <- c(200, 200)
  dest   <- c(2800, 2800)

  res_dijk  <- route_corridor(r, origin, dest, method = "dijkstra")
  res_bidir <- route_corridor(r, origin, dest, method = "bidirectional")
  res_astar <- route_corridor(r, origin, dest, method = "astar")

  expect_equal(res_dijk$total_cost, res_bidir$total_cost, tolerance = 1e-6)
  expect_equal(res_dijk$total_cost, res_astar$total_cost, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 9. Resolution factor
# ---------------------------------------------------------------------------
test_that("resolution_factor produces valid but coarser path", {
  r <- make_test_raster(50, 50, vals = 1)
  set.seed(42)
  values(r) <- runif(ncell(r), 0.5, 2.0)

  origin <- c(500, 500)
  dest   <- c(4500, 4500)

  path_fine   <- route_corridor(r, origin, dest, resolution_factor = 1L)
  path_coarse <- route_corridor(r, origin, dest, resolution_factor = 2L)

  # Coarse path should have fewer cells
  expect_lt(path_coarse$n_cells, path_fine$n_cells)
  # Both should be valid
  expect_s3_class(path_coarse, "spopt_corridor")
})

# ---------------------------------------------------------------------------
# 10. Endpoint geometry
# ---------------------------------------------------------------------------
test_that("linestring starts and ends at user-supplied coordinates", {
  r <- make_test_raster(50, 50, vals = 1)
  origin <- c(123.4, 567.8)
  dest   <- c(4321.0, 3456.7)

  path <- route_corridor(r, origin, dest)
  coords <- sf::st_coordinates(path)

  # First point = origin
  expect_equal(as.numeric(coords[1, "X"]), origin[1], tolerance = 1e-10)
  expect_equal(as.numeric(coords[1, "Y"]), origin[2], tolerance = 1e-10)

  # Last point = destination
  n <- nrow(coords)
  expect_equal(as.numeric(coords[n, "X"]), dest[1], tolerance = 1e-10)
  expect_equal(as.numeric(coords[n, "Y"]), dest[2], tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 11. Neighbour ordering: 16 <= 8 <= 4 cost
# ---------------------------------------------------------------------------
test_that("more connectivity produces lower or equal cost", {
  r <- make_test_raster(30, 30, vals = 1)
  set.seed(42)
  values(r) <- runif(ncell(r), 0.5, 2.0)

  origin <- c(200, 200)
  dest   <- c(2800, 2800)

  cost_4  <- route_corridor(r, origin, dest, neighbours = 4L)$total_cost
  cost_8  <- route_corridor(r, origin, dest, neighbours = 8L)$total_cost
  cost_16 <- route_corridor(r, origin, dest, neighbours = 16L)$total_cost

  expect_lte(cost_16, cost_8 + 1e-6)
  expect_lte(cost_8, cost_4 + 1e-6)
})

# ---------------------------------------------------------------------------
# 12. CRS transform of sf input
# ---------------------------------------------------------------------------
test_that("sf point in different CRS is auto-transformed", {
  r <- make_test_raster(50, 50, vals = 1, crs = "EPSG:32618")

  # Create origin and destination as sf points in the raster CRS
  origin_utm <- st_sfc(st_point(c(500, 500)), crs = 32618)
  dest_utm   <- st_sfc(st_point(c(4500, 4500)), crs = 32618)

  # Transform origin to lonlat; route_corridor should auto-transform back
  origin_ll <- st_transform(origin_utm, 4326)

  path <- route_corridor(r, origin_ll, dest_utm)
  expect_s3_class(path, "spopt_corridor")
})

# ---------------------------------------------------------------------------
# 13. CRS-less sf rejection
# ---------------------------------------------------------------------------
test_that("sf point without CRS is rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  pt_no_crs <- st_sfc(st_point(c(500, 500)))

  expect_error(
    route_corridor(r, pt_no_crs, c(900, 900)),
    "no CRS"
  )
})

# ---------------------------------------------------------------------------
# 14. Non-POINT rejection
# ---------------------------------------------------------------------------
test_that("non-POINT geometry is rejected", {
  r <- make_test_raster(10, 10, vals = 1)

  # Multipoint
  mp <- st_sfc(st_multipoint(matrix(c(500, 500, 900, 900), ncol = 2, byrow = TRUE)),
               crs = 32618)
  expect_error(
    route_corridor(r, mp, c(900, 900)),
    "POINT"
  )

  # Multi-row sf
  pts <- st_as_sf(data.frame(x = c(100, 200), y = c(100, 200)),
                  coords = c("x", "y"), crs = 32618)
  expect_error(
    route_corridor(r, pts, c(900, 900)),
    "single point"
  )
})

# ---------------------------------------------------------------------------
# 15. No path (complete NA barrier)
# ---------------------------------------------------------------------------
test_that("no-path scenario raises error", {
  r <- make_test_raster(20, 20, vals = 1)
  # Complete horizontal NA barrier
  cells <- cellFromRowCol(r, rep(10, 20), 1:20)
  r[cells] <- NA

  expect_error(
    route_corridor(r, c(500, 500), c(500, 1500)),
    "No path"
  )
})

# ---------------------------------------------------------------------------
# 16. Origin == destination (degenerate case)
# ---------------------------------------------------------------------------
test_that("origin == destination returns valid result with NA sinuosity", {
  r <- make_test_raster(10, 10, vals = 1)
  pt <- c(500, 500)

  path <- route_corridor(r, pt, pt)

  expect_s3_class(path, "spopt_corridor")
  expect_true(is.na(path$sinuosity))
  expect_equal(path$total_cost, 0, tolerance = 1e-10)
  expect_equal(path$straight_line_dist, 0, tolerance = 1e-10)
  expect_equal(path$path_dist, 0, tolerance = 1e-10)
})

# ===========================================================================
# Graph caching tests
# ===========================================================================

# ---------------------------------------------------------------------------
# 17. Semantic equivalence: direct vs cached
# ---------------------------------------------------------------------------
test_that("cached graph produces identical results to direct path", {
  set.seed(42)
  r <- make_test_raster(50, 50, vals = runif(2500, 0.5, 2.0))
  o <- c(500, 500)
  d <- c(4500, 4500)

  direct <- route_corridor(r, o, d, method = "astar")
  g      <- corridor_graph(r, neighbours = 8L)
  cached <- route_corridor(g, o, d, method = "astar")

  expect_equal(cached$total_cost, direct$total_cost, tolerance = 1e-10)
  expect_equal(
    attr(cached, "spopt")$cell_indices,
    attr(direct, "spopt")$cell_indices
  )
  expect_equal(sf::st_coordinates(cached), sf::st_coordinates(direct))
})

# ---------------------------------------------------------------------------
# 18. Multiple OD pairs on one graph
# ---------------------------------------------------------------------------
test_that("multiple OD pairs work on a single cached graph", {
  r <- make_test_raster(50, 50, vals = 1)
  g <- corridor_graph(r, neighbours = 8L)

  p1 <- route_corridor(g, c(500, 500), c(4500, 4500), method = "dijkstra")
  p2 <- route_corridor(g, c(500, 2500), c(4500, 2500), method = "dijkstra")
  p3 <- route_corridor(g, c(2500, 500), c(2500, 4500), method = "dijkstra")

  expect_s3_class(p1, "spopt_corridor")
  expect_s3_class(p2, "spopt_corridor")
  expect_s3_class(p3, "spopt_corridor")

  # graph_build_time should be 0 for cached routes
  expect_equal(attr(p1, "spopt")$graph_build_time, 0)
})

# ---------------------------------------------------------------------------
# 19. Method equivalence on cached graph
# ---------------------------------------------------------------------------
test_that("all 3 methods produce same cost on cached graph", {
  set.seed(99)
  r <- make_test_raster(30, 30, vals = runif(900, 0.5, 2.0))
  g <- corridor_graph(r, neighbours = 8L)

  o <- c(200, 200)
  d <- c(2800, 2800)

  c1 <- route_corridor(g, o, d, method = "dijkstra")$total_cost
  c2 <- route_corridor(g, o, d, method = "bidirectional")$total_cost
  c3 <- route_corridor(g, o, d, method = "astar")$total_cost

  expect_equal(c1, c2, tolerance = 1e-6)
  expect_equal(c1, c3, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 20. Resolution factor with graph
# ---------------------------------------------------------------------------
test_that("corridor_graph works with resolution_factor", {
  r <- make_test_raster(50, 50, vals = 1)
  g <- corridor_graph(r, neighbours = 8L, resolution_factor = 2L)

  path <- route_corridor(g, c(500, 500), c(4500, 4500))
  expect_s3_class(path, "spopt_corridor")
  # Coarser grid means fewer cells
  expect_lt(path$n_cells, 50)
})

# ---------------------------------------------------------------------------
# 21. Print method for corridor graph
# ---------------------------------------------------------------------------
test_that("print.spopt_corridor_graph produces expected output", {
  r <- make_test_raster(50, 50, vals = 1)
  g <- corridor_graph(r, neighbours = 8L)

  out <- capture.output(print(g))
  expect_true(any(grepl("Corridor graph", out)))
  expect_true(any(grepl("50 x 50", out)))
  expect_true(any(grepl("Neighbours: 8", out)))
})

# ---------------------------------------------------------------------------
# 22. Pointer invalidation
# ---------------------------------------------------------------------------
test_that("invalidated pointer gives informative error", {
  r <- make_test_raster(10, 10, vals = 1)
  g <- corridor_graph(r, neighbours = 8L)
  g$ptr <- NULL

  expect_error(
    route_corridor(g, c(500, 500), c(900, 900)),
    "invalid|deserialized"
  )
})

# ---------------------------------------------------------------------------
# 23. Graph stores neighbours
# ---------------------------------------------------------------------------
test_that("graph metadata stores the correct neighbours", {
  r <- make_test_raster(10, 10, vals = 1)
  g4  <- corridor_graph(r, neighbours = 4L)
  g16 <- corridor_graph(r, neighbours = 16L)

  expect_equal(attr(g4, "spopt")$neighbours, 4L)
  expect_equal(attr(g16, "spopt")$neighbours, 16L)
})

# ---------------------------------------------------------------------------
# 24. Override rejection
# ---------------------------------------------------------------------------
test_that("overriding neighbours on cached graph errors", {
  r <- make_test_raster(10, 10, vals = 1)
  g <- corridor_graph(r, neighbours = 8L)

  expect_error(
    route_corridor(g, c(500, 500), c(900, 900), neighbours = 4L),
    "cannot be overridden"
  )
  expect_error(
    route_corridor(g, c(500, 500), c(900, 900), resolution_factor = 2L),
    "cannot be overridden"
  )
})

# ---------------------------------------------------------------------------
# 25. NA cell on cached graph
# ---------------------------------------------------------------------------
test_that("NA cell check works on cached graph", {
  r <- make_test_raster(10, 10, vals = 1)
  r[cellFromXY(r, matrix(c(500, 500), ncol = 2))] <- NA
  g <- corridor_graph(r, neighbours = 8L)

  expect_error(
    route_corridor(g, c(500, 500), c(900, 900)),
    "impassable|NA"
  )
})

# ===========================================================================
# Waypoints (ordered multi-city routing)
# ===========================================================================

# ---------------------------------------------------------------------------
# W1. Regression: waypoints=NULL retains legacy schema and equivalent output
# ---------------------------------------------------------------------------
test_that("waypoints=NULL preserves legacy schema and semantics", {
  set.seed(17)
  r <- make_test_raster(40, 40, vals = runif(1600, 0.5, 2.0))
  o <- c(500, 500); d <- c(3500, 3500)

  path <- route_corridor(r, o, d)

  expect_s3_class(path, "spopt_corridor")
  expect_s3_class(path, "sf")
  expect_equal(nrow(path), 1L)
  expect_setequal(
    setdiff(names(path), c("geometry")),
    c("total_cost", "n_cells", "straight_line_dist", "path_dist", "sinuosity")
  )
  meta <- attr(path, "spopt")
  # Waypoint fields must be absent in the no-waypoint case
  expect_null(meta$n_waypoints_input)
  expect_null(meta$segment_costs)
  expect_null(meta$waypoints_input_xy)
})

# ---------------------------------------------------------------------------
# W2. Equivalence to manual chaining
# ---------------------------------------------------------------------------
test_that("waypoint total_cost equals manual chain of two route_corridor calls", {
  set.seed(42)
  r <- make_test_raster(50, 50, vals = runif(2500, 0.5, 2.0))
  o <- c(500, 500); w <- c(2500, 2500); d <- c(4500, 4500)

  leg1 <- route_corridor(r, o, w)
  leg2 <- route_corridor(r, w, d)
  manual_total <- leg1$total_cost + leg2$total_cost

  via <- route_corridor(r, o, d, waypoints = list(w))

  expect_equal(via$total_cost, manual_total, tolerance = 1e-9)
  expect_equal(attr(via, "spopt")$n_waypoints_input, 1L)
  expect_equal(attr(via, "spopt")$n_waypoints_effective, 1L)
  expect_equal(attr(via, "spopt")$n_segments_effective, 2L)
  expect_equal(length(attr(via, "spopt")$segment_costs), 2L)
  expect_equal(sum(attr(via, "spopt")$segment_costs),
               via$total_cost, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# W3. No backtrack spike at off-cell-center waypoint
# ---------------------------------------------------------------------------
test_that("off-cell-center waypoint does not produce a backtrack vertex", {
  r <- make_test_raster(50, 50, vals = 1)
  # Cell size is 100; waypoint at an off-center coord
  via <- route_corridor(r, c(250, 250), c(4750, 4750),
                        waypoints = list(c(2537, 2523)))
  coords <- sf::st_coordinates(via)[, c("X", "Y")]
  # Assert no v_i where v_{i-1} == v_{i+1} and v_{i-1} != v_i
  backtracks <- 0L
  for (i in 2:(nrow(coords) - 1)) {
    if (all(coords[i - 1L, ] == coords[i + 1L, ]) &&
        !all(coords[i - 1L, ] == coords[i, ])) {
      backtracks <- backtracks + 1L
    }
  }
  expect_equal(backtracks, 0L)
})

# ---------------------------------------------------------------------------
# W4. Two waypoints: segment_costs sum to total
# ---------------------------------------------------------------------------
test_that("two waypoints: segment_costs sum equals total_cost", {
  set.seed(3)
  r <- make_test_raster(40, 40, vals = runif(1600, 0.5, 2.0))
  via <- route_corridor(r, c(250, 250), c(3750, 3750),
                        waypoints = list(c(1250, 1750), c(2750, 1500)))
  meta <- attr(via, "spopt")
  expect_equal(meta$n_segments_effective, 3L)
  expect_equal(length(meta$segment_costs), 3L)
  expect_equal(sum(meta$segment_costs), via$total_cost, tolerance = 1e-9)
  expect_equal(length(meta$segment_solve_times), 3L)
})

# ---------------------------------------------------------------------------
# W5. Raster input and cached-graph input agree
# ---------------------------------------------------------------------------
test_that("raster and cached-graph inputs produce matching waypoint routes", {
  set.seed(5)
  r <- make_test_raster(40, 40, vals = runif(1600, 0.5, 2.0))
  g <- corridor_graph(r, neighbours = 8L)
  wp <- list(c(1250, 1750), c(2750, 1500))

  via_r <- route_corridor(r, c(250, 250), c(3750, 3750), waypoints = wp)
  via_g <- route_corridor(g, c(250, 250), c(3750, 3750), waypoints = wp)

  expect_equal(via_g$total_cost, via_r$total_cost, tolerance = 1e-10)
  expect_equal(attr(via_g, "spopt")$cell_indices,
               attr(via_r, "spopt")$cell_indices)
  expect_equal(sf::st_coordinates(via_g),
               sf::st_coordinates(via_r))
})

# ---------------------------------------------------------------------------
# W6. graph_build_time recorded once; solve times per-segment
# ---------------------------------------------------------------------------
test_that("graph_build_time is a single positive scalar for transient graph", {
  set.seed(7)
  r <- make_test_raster(40, 40, vals = runif(1600, 0.5, 2.0))
  via <- route_corridor(r, c(250, 250), c(3750, 3750),
                        waypoints = list(c(1500, 2500)))
  meta <- attr(via, "spopt")

  expect_true(is.numeric(meta$graph_build_time))
  expect_equal(length(meta$graph_build_time), 1L)
  expect_gte(meta$graph_build_time, 0)
  expect_equal(length(meta$segment_solve_times), meta$n_segments_effective)
  expect_true(all(meta$segment_solve_times >= 0))
})

# ---------------------------------------------------------------------------
# W7. Waypoint on NA cell -> error naming waypoint index
# ---------------------------------------------------------------------------
test_that("waypoint on NA cell errors with waypoint index", {
  r <- make_test_raster(20, 20, vals = 1)
  r[cellFromXY(r, matrix(c(1000, 1000), ncol = 2))] <- NA
  expect_error(
    route_corridor(r, c(250, 250), c(1750, 1750),
                   waypoints = list(c(1000, 1000))),
    "waypoint 1"
  )
})

# ---------------------------------------------------------------------------
# W8. Waypoint outside raster extent -> error
# ---------------------------------------------------------------------------
test_that("waypoint outside raster extent errors with waypoint index", {
  r <- make_test_raster(20, 20, vals = 1)
  expect_error(
    route_corridor(r, c(250, 250), c(1750, 1750),
                   waypoints = list(c(10000, 10000))),
    "waypoint 1"
  )
})

# ---------------------------------------------------------------------------
# W9. Exact-duplicate waypoint -> warning, still valid
# ---------------------------------------------------------------------------
test_that("exact duplicate waypoint emits warning and dedups", {
  r <- make_test_raster(20, 20, vals = 1)
  wp <- c(1000, 1000)
  expect_warning(
    via <- route_corridor(r, c(250, 250), c(1750, 1750),
                          waypoints = list(wp, wp)),
    "identical"
  )
  meta <- attr(via, "spopt")
  expect_equal(meta$n_waypoints_input, 2L)
  # One waypoint dropped as duplicate -> effective is 1
  expect_equal(meta$n_waypoints_effective, 1L)
})

# ---------------------------------------------------------------------------
# W10. Same-cell consecutive points -> warning
# ---------------------------------------------------------------------------
test_that("close-but-distinct waypoints snapping to same cell emit warning", {
  r <- make_test_raster(20, 20, vals = 1)
  # Two coords unambiguously within the same 100x100 cell (1000..1100 x 1000..1100)
  expect_warning(
    route_corridor(r, c(250, 250), c(1750, 1750),
                   waypoints = list(c(1050, 1050), c(1080, 1060))),
    "same raster cell"
  )
})

# ---------------------------------------------------------------------------
# W11. output = "segments" returns N+1 rows with correct columns
# ---------------------------------------------------------------------------
test_that("output='segments' returns N+1 rows with leg columns", {
  r <- make_test_raster(30, 30, vals = 1)
  via <- route_corridor(r, c(250, 250), c(2750, 2750),
                        waypoints = list(c(1250, 1750)),
                        output = "segments")
  expect_s3_class(via, "spopt_corridor_segments")
  expect_s3_class(via, "sf")
  expect_equal(nrow(via), 2L)
  expected_cols <- c("segment", "from_label", "to_label",
                     "from_x", "from_y", "to_x", "to_y",
                     "total_cost", "n_cells", "path_dist",
                     "straight_line_dist", "sinuosity")
  expect_true(all(expected_cols %in% names(via)))
  expect_equal(via$from_label, c("origin", "waypoint 1"))
  expect_equal(via$to_label,   c("waypoint 1", "destination"))
})

# ---------------------------------------------------------------------------
# W12. output = "combined" with waypoints returns single-row sf
# ---------------------------------------------------------------------------
test_that("output='combined' with waypoints returns 1-row sf", {
  r <- make_test_raster(30, 30, vals = 1)
  via <- route_corridor(r, c(250, 250), c(2750, 2750),
                        waypoints = list(c(1250, 1750)),
                        output = "combined")
  expect_s3_class(via, "spopt_corridor")
  expect_equal(nrow(via), 1L)
})

# ---------------------------------------------------------------------------
# W13. Legacy spopt keys preserved in combined-mode waypoint output
# ---------------------------------------------------------------------------
test_that("existing spopt attribute keys are preserved with waypoints", {
  r <- make_test_raster(30, 30, vals = 1)
  via <- route_corridor(r, c(250, 250), c(2750, 2750),
                        waypoints = list(c(1250, 1750)))
  meta <- attr(via, "spopt")
  legacy_keys <- c("total_cost", "n_cells", "method", "neighbours",
                   "n_cells_surface", "n_edges_graph", "solve_time",
                   "graph_build_time", "cell_indices",
                   "origin_cell_center", "dest_cell_center",
                   "raster_dims", "cell_size")
  expect_true(all(legacy_keys %in% names(meta)))
})

# ---------------------------------------------------------------------------
# W14. Exact-supplied and snapped waypoint coords both stored
# ---------------------------------------------------------------------------
test_that("waypoints_input_xy preserves exact coords, waypoint_cell_xy is snapped", {
  r <- make_test_raster(30, 30, vals = 1)
  # Cell size = 100; off-center coord will snap to a different point
  wp <- c(1537, 1523)
  via <- route_corridor(r, c(250, 250), c(2750, 2750),
                        waypoints = list(wp))
  meta <- attr(via, "spopt")
  expect_equal(as.numeric(meta$waypoints_input_xy[1, ]), wp)
  expect_equal(dim(meta$waypoints_input_xy), dim(meta$waypoint_cell_xy))
  # Exact should differ from snapped for an off-center waypoint
  expect_false(isTRUE(all.equal(as.numeric(meta$waypoints_input_xy[1, ]),
                                 as.numeric(meta$waypoint_cell_xy[1, ]))))
})

# ---------------------------------------------------------------------------
# W15. leg_* fields behave correctly
# ---------------------------------------------------------------------------
test_that("leg_straight_line_dist and leg_sinuosity populate with waypoints", {
  r <- make_test_raster(30, 30, vals = 1)
  via <- route_corridor(r, c(250, 250), c(2750, 2750),
                        waypoints = list(c(1000, 2250)))

  expect_true("leg_straight_line_dist" %in% names(via))
  expect_true("leg_sinuosity" %in% names(via))
  # Leg sum should be >= direct (detour through a waypoint lengthens the baseline)
  expect_gte(via$leg_straight_line_dist, via$straight_line_dist - 1e-8)
  # path_dist >= leg_straight_line_dist (routed path can't be shorter than its required baseline)
  expect_gte(via$path_dist, via$leg_straight_line_dist - 1e-8)
})

# ---------------------------------------------------------------------------
# W16. resolution_factor with waypoints -> no double-aggregation
# ---------------------------------------------------------------------------
test_that("resolution_factor>1 with waypoints doesn't double-aggregate", {
  set.seed(11)
  r <- make_test_raster(60, 60, vals = runif(3600, 0.5, 2.0))
  g <- corridor_graph(r, neighbours = 8L, resolution_factor = 2L)
  n_edges_ref <- attr(g, "spopt")$n_edges

  via <- route_corridor(r, c(500, 500), c(5500, 5500),
                        waypoints = list(c(3000, 2000)),
                        resolution_factor = 2L)
  expect_equal(attr(via, "spopt")$n_edges_graph, n_edges_ref)
})

# ---------------------------------------------------------------------------
# W17. sf POINTs input for waypoints works
# ---------------------------------------------------------------------------
test_that("sf/sfc POINT input for waypoints works with auto-transform", {
  r <- make_test_raster(30, 30, vals = 1, crs = "EPSG:32618")

  wp_sfc <- sf::st_sfc(sf::st_point(c(1250, 1750)), crs = 32618)
  via <- route_corridor(r, c(250, 250), c(2750, 2750), waypoints = wp_sfc)
  expect_s3_class(via, "spopt_corridor")
  expect_equal(attr(via, "spopt")$n_waypoints_input, 1L)

  # Multi-feature sf
  wp_sf <- sf::st_as_sf(data.frame(x = c(1000, 2000), y = c(1500, 2000)),
                        coords = c("x", "y"), crs = 32618)
  via2 <- route_corridor(r, c(250, 250), c(2750, 2750), waypoints = wp_sf)
  expect_equal(attr(via2, "spopt")$n_waypoints_input, 2L)
})

# ---------------------------------------------------------------------------
# W18. Destination is preserved when a waypoint snaps to destination's cell
# ---------------------------------------------------------------------------
# Regression: same-cell elision must not drop the user's destination.
test_that("waypoint sharing destination's cell does not overwrite destination", {
  r <- make_test_raster(10, 10, vals = 1)  # cell size 100, extent 0..1000
  # Cells (905,905) and (950,950) both fall in the same cell (col 10, row 1)
  expect_warning(
    via <- route_corridor(r, c(50, 50), c(950, 950),
                          waypoints = list(c(905, 905))),
    "same raster cell"
  )
  coords <- sf::st_coordinates(via)
  n <- nrow(coords)
  expect_equal(as.numeric(coords[1, "X"]), 50, tolerance = 1e-10)
  expect_equal(as.numeric(coords[1, "Y"]), 50, tolerance = 1e-10)
  expect_equal(as.numeric(coords[n, "X"]), 950, tolerance = 1e-10)
  expect_equal(as.numeric(coords[n, "Y"]), 950, tolerance = 1e-10)
  # Waypoint is elided, so effective is 0
  meta <- attr(via, "spopt")
  expect_equal(meta$n_waypoints_input, 1L)
  expect_equal(meta$n_waypoints_effective, 0L)
})

# ---------------------------------------------------------------------------
# W19. graph_build_time reflects real construction time for single-segment paths
# ---------------------------------------------------------------------------
# Regression: segments-mode on a raster without waypoints still does a Rust
# graph build internally; that time must not be silently zeroed.
test_that("segments-mode on raster reports non-zero graph_build_time", {
  set.seed(23)
  r <- make_test_raster(40, 40, vals = runif(1600, 0.5, 2.0))
  via <- route_corridor(r, c(250, 250), c(3750, 3750),
                        waypoints = list(c(1500, 2500)),
                        output = "segments")
  meta <- attr(via, "spopt")
  expect_true(meta$graph_build_time > 0)
})

test_that("collapsed-to-one-segment route reports real graph_build_time", {
  r <- make_test_raster(10, 10, vals = 1)
  # waypoint snaps to destination cell -> effective 0 waypoints, single segment
  suppressWarnings(
    via <- route_corridor(r, c(50, 50), c(950, 950),
                          waypoints = list(c(905, 905)))
  )
  expect_true(attr(via, "spopt")$graph_build_time > 0)
})

# ---------------------------------------------------------------------------
# W20. Segments-mode aggregate n_cells dedups shared join cells
# ---------------------------------------------------------------------------
# Regression: aggregate n_cells in segments mode must equal the deduped cell
# count from combined mode (shared join cells not double-counted).
test_that("segments-mode n_cells matches combined-mode n_cells", {
  set.seed(31)
  r <- make_test_raster(30, 30, vals = runif(900, 0.5, 2.0))
  wp <- list(c(1250, 1750))

  via_c <- route_corridor(r, c(250, 250), c(2750, 2750),
                          waypoints = wp, output = "combined")
  via_s <- route_corridor(r, c(250, 250), c(2750, 2750),
                          waypoints = wp, output = "segments")
  expect_equal(attr(via_s, "spopt")$n_cells, attr(via_c, "spopt")$n_cells)
})

# ---------------------------------------------------------------------------
# W21b. Segments-mode spopt attr carries the same aggregate metadata keys
#       as combined-mode (cell_indices, origin_cell_center, dest_cell_center)
# ---------------------------------------------------------------------------
test_that("segments-mode spopt attr matches combined-mode aggregate metadata", {
  set.seed(41)
  r <- make_test_raster(30, 30, vals = runif(900, 0.5, 2.0))
  wp <- list(c(1250, 1750))

  via_c <- route_corridor(r, c(250, 250), c(2750, 2750),
                          waypoints = wp, output = "combined")
  via_s <- route_corridor(r, c(250, 250), c(2750, 2750),
                          waypoints = wp, output = "segments")
  mc <- attr(via_c, "spopt")
  ms <- attr(via_s, "spopt")

  # The three keys codex flagged as missing from segments mode
  expect_equal(ms$cell_indices, mc$cell_indices)
  expect_equal(ms$origin_cell_center, mc$origin_cell_center)
  expect_equal(ms$dest_cell_center, mc$dest_cell_center)
  # leg_* metrics should also match combined-mode values
  expect_equal(ms$leg_straight_line_dist, mc$leg_straight_line_dist,
               tolerance = 1e-9)
  expect_equal(ms$leg_sinuosity, mc$leg_sinuosity, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# W21. waypoints_input_xy preserves the full supplied list after elision
# ---------------------------------------------------------------------------
# Regression: even when duplicates or same-cell waypoints are elided,
# waypoints_input_xy must have n_waypoints_input rows (original user list).
test_that("waypoints_input_xy retains all supplied points after elision", {
  r <- make_test_raster(20, 20, vals = 1)
  wp <- c(1000, 1000)
  suppressWarnings(
    via <- route_corridor(r, c(250, 250), c(1750, 1750),
                          waypoints = list(wp, wp, c(1250, 1500)))
  )
  meta <- attr(via, "spopt")
  # 3 supplied, but two exact duplicates are elided to 1 effective waypoint
  expect_equal(meta$n_waypoints_input, 3L)
  expect_equal(nrow(meta$waypoints_input_xy), 3L)
  # waypoint_cell_xy is the effective count, not the input count
  expect_equal(nrow(meta$waypoint_cell_xy), meta$n_waypoints_effective)
})
