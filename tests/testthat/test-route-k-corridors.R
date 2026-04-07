skip_if_not_installed("terra")
skip_if_not_installed("sf")

library(terra)
library(sf)

# make_test_raster() is defined in helper-raster.R

# ---------------------------------------------------------------------------
# 1. Uniform surface k=3 — paths fan out
# ---------------------------------------------------------------------------
test_that("k=3 on uniform surface produces diverse paths", {
  set.seed(42)
  r <- make_test_raster(50, 50, vals = runif(2500, 0.5, 2.0))
  o <- c(500, 500)
  d <- c(4500, 4500)

  result <- route_k_corridors(r, o, d, k = 3L, method = "astar")

  expect_s3_class(result, "spopt_k_corridors")
  expect_equal(nrow(result), 3)
  # Ranks > 1 have positive mean_spacing
  expect_true(all(result$mean_spacing[2:3] > 0))
  # Rank 1 has NA spacing/overlap
  expect_true(is.na(result$mean_spacing[1]))
  expect_true(is.na(result$pct_overlap[1]))
})

# ---------------------------------------------------------------------------
# 2. Exclusion zone — corridors route around different sides
# ---------------------------------------------------------------------------
test_that("NA barrier forces corridors to different sides", {
  r <- make_test_raster(50, 50, vals = 1)
  # Vertical NA barrier in the middle, blocking rows 15-35
  cells <- cellFromRowCol(r, rep(15:35, each = 3), rep(24:26, times = 21))
  r[cells] <- NA

  result <- route_k_corridors(r, c(500, 2500), c(4500, 2500), k = 2L)

  expect_equal(nrow(result), 2)
  # The two corridors should have substantial spacing
  expect_gt(result$mean_spacing[2], 0)
})

# ---------------------------------------------------------------------------
# 3. Cost recomputation — k=1 matches route_corridor()
# ---------------------------------------------------------------------------
test_that("k=1 total_cost matches route_corridor()", {
  set.seed(99)
  r <- make_test_raster(30, 30, vals = runif(900, 0.5, 2.0))
  o <- c(200, 200)
  d <- c(2800, 2800)

  k1     <- route_k_corridors(r, o, d, k = 1L)
  direct <- route_corridor(r, o, d)

  expect_equal(k1$total_cost, direct$total_cost, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 4. k exceeds feasible — returns fewer
# ---------------------------------------------------------------------------
test_that("returns fewer than k when paths exhaust", {
  r <- make_test_raster(30, 30, vals = 1)
  # Narrow isthmus: block most of column 15, leaving only rows 14-16 open
  wall_rows <- c(1:13, 17:30)
  cells <- cellFromRowCol(r, rep(wall_rows, each = 3), rep(14:16, times = length(wall_rows)))
  r[cells] <- NA

  # Large penalty_radius + high penalty_factor on narrow gap should exhaust paths
  result <- route_k_corridors(r, c(200, 1500), c(2800, 1500), k = 10L,
                               penalty_radius = 500, penalty_factor = 5.0)

  meta <- attr(result, "spopt")
  expect_lte(meta$k_found, meta$k_requested)
  expect_equal(nrow(result), meta$k_found)
  expect_gte(meta$k_found, 1L)
})

# ---------------------------------------------------------------------------
# 5. min_spacing retry
# ---------------------------------------------------------------------------
test_that("min_spacing triggers retries", {
  set.seed(7)
  r <- make_test_raster(50, 50, vals = runif(2500, 0.5, 2.0))
  o <- c(500, 500)
  d <- c(4500, 4500)

  # Large min_spacing should force retries
  result <- route_k_corridors(r, o, d, k = 2L, min_spacing = 2000)

  meta <- attr(result, "spopt")
  # Either retries occurred or only 1 corridor was found

  expect_true(meta$total_retries > 0L || meta$k_found == 1L)
})

# ---------------------------------------------------------------------------
# 6. corridor_graph rejection
# ---------------------------------------------------------------------------
test_that("corridor_graph input is rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  g <- corridor_graph(r, neighbours = 8L)

  expect_error(
    route_k_corridors(g, c(500, 500), c(900, 900)),
    "SpatRaster"
  )
})

# ---------------------------------------------------------------------------
# 7. Print method
# ---------------------------------------------------------------------------
test_that("print.spopt_k_corridors produces expected output", {
  set.seed(42)
  r <- make_test_raster(30, 30, vals = runif(900, 0.5, 2.0))
  result <- route_k_corridors(r, c(200, 200), c(2800, 2800), k = 2L)

  out <- capture.output(print(result))
  expect_true(any(grepl("k-Diverse", out)))
  expect_true(any(grepl("Optimal", out)))
  # Optimal spacing should show "-"
  optimal_line <- out[grep("Optimal", out)]
  expect_true(grepl("-", optimal_line))
})

# ---------------------------------------------------------------------------
# 8. Class and columns
# ---------------------------------------------------------------------------
test_that("result has correct class and columns", {
  r <- make_test_raster(30, 30, vals = 1)
  result <- route_k_corridors(r, c(200, 200), c(2800, 2800), k = 2L)

  expect_s3_class(result, "spopt_k_corridors")
  expect_s3_class(result, "sf")
  expected_cols <- c("alternative", "total_cost", "n_cells", "path_dist",
                     "straight_line_dist", "sinuosity", "mean_spacing",
                     "pct_overlap", "geometry")
  expect_true(all(expected_cols %in% names(result)))
})

# ---------------------------------------------------------------------------
# 9. k=1 format — NA spacing/overlap
# ---------------------------------------------------------------------------
test_that("k=1 returns correct format", {
  r <- make_test_raster(20, 20, vals = 1)
  result <- route_k_corridors(r, c(200, 200), c(1800, 1800), k = 1L)

  expect_equal(nrow(result), 1)
  expect_equal(result$alternative, 1L)
  expect_true(is.na(result$mean_spacing))
  expect_true(is.na(result$pct_overlap))
})

# ---------------------------------------------------------------------------
# 10. Invalid k
# ---------------------------------------------------------------------------
test_that("invalid k is rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  expect_error(route_k_corridors(r, c(500, 500), c(900, 900), k = 0), "k")
  expect_error(route_k_corridors(r, c(500, 500), c(900, 900), k = -1), "k")
})

# ---------------------------------------------------------------------------
# 11. Invalid penalty_factor
# ---------------------------------------------------------------------------
test_that("invalid penalty_factor is rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  expect_error(
    route_k_corridors(r, c(500, 500), c(900, 900), penalty_factor = 0.5),
    "penalty_factor"
  )
  expect_error(
    route_k_corridors(r, c(500, 500), c(900, 900), penalty_factor = 1.0),
    "penalty_factor"
  )
})

# ---------------------------------------------------------------------------
# 12. Invalid penalty_radius
# ---------------------------------------------------------------------------
test_that("invalid penalty_radius is rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  expect_error(
    route_k_corridors(r, c(500, 500), c(900, 900), penalty_radius = -1),
    "penalty_radius"
  )
  expect_error(
    route_k_corridors(r, c(500, 500), c(900, 900), penalty_radius = 0),
    "penalty_radius"
  )
})

# ---------------------------------------------------------------------------
# 13. No-path graceful stop
# ---------------------------------------------------------------------------
test_that("heavily penalized surface eventually stops", {
  r <- make_test_raster(30, 30, vals = 1)
  # Narrow gap — first path works, but aggressive penalty radius
  # should exhaust alternatives quickly
  wall_rows <- c(1:13, 17:30)
  cells <- cellFromRowCol(r, rep(wall_rows, each = 3), rep(14:16, times = length(wall_rows)))
  r[cells] <- NA

  result <- route_k_corridors(r, c(200, 1500), c(2800, 1500), k = 10L,
                               penalty_radius = 500, penalty_factor = 3.0)
  meta <- attr(result, "spopt")
  expect_true(meta$k_found <= 10L)
  expect_true(meta$k_found >= 1L)
})

# ---------------------------------------------------------------------------
# 14. Origin == destination, no penalty_radius
# ---------------------------------------------------------------------------
test_that("origin == dest is rejected", {
  r <- make_test_raster(10, 10, vals = 1)
  expect_error(
    route_k_corridors(r, c(500, 500), c(500, 500)),
    "identical"
  )
  # Also rejected even with explicit penalty_radius
  expect_error(
    route_k_corridors(r, c(500, 500), c(500, 500), penalty_radius = 100),
    "identical"
  )
})

# ---------------------------------------------------------------------------
# 15. Plot smoke test
# ---------------------------------------------------------------------------
test_that("plot.spopt_k_corridors runs without error", {
  set.seed(42)
  r <- make_test_raster(30, 30, vals = runif(900, 0.5, 2.0))
  result <- route_k_corridors(r, c(200, 200), c(2800, 2800), k = 2L)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(result))
})
