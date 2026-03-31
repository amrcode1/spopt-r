skip_if_not_installed("terra")
skip_if_not_installed("sf")

library(terra)
library(sf)

# Helper: create a projected raster with uniform values
make_test_raster <- function(nrow = 50, ncol = 50, vals = 1,
                             crs = "EPSG:32618") {
  r <- rast(
    nrows = nrow, ncols = ncol,
    xmin = 0, xmax = ncol * 100,
    ymin = 0, ymax = nrow * 100,
    crs = crs
  )
  values(r) <- vals
  r
}

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
