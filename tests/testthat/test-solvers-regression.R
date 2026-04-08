# Backend-agnostic regression tests for facility location solvers.
#
# These pin known-good objective values from the Rust backend (captured
# April 2026) and assert structural invariants that any correct backend
# must satisfy. Exact facility selections are NOT pinned because
# tie-breaking may differ across backends.
#
# Tolerances vary by solver: 1e-1 for p_median/mclp, 1e-2 for
# p_center/p_dispersion, ~1 for cflp. Loose enough to survive
# backend changes, tight enough to catch model formulation errors.

skip_if_not_installed("sf")
library(sf)

# ---------------------------------------------------------------------------
# Shared fixture helper
# ---------------------------------------------------------------------------
make_solver_data <- function(n_demand, n_fac, seed = 123) {
  set.seed(seed)
  demand <- st_as_sf(
    data.frame(
      x = runif(n_demand, 0, 10),
      y = runif(n_demand, 0, 10),
      pop = rpois(n_demand, 50) + 1L
    ),
    coords = c("x", "y")
  )
  facilities <- st_as_sf(
    data.frame(
      x = runif(n_fac, 0, 10),
      y = runif(n_fac, 0, 10)
    ),
    coords = c("x", "y")
  )
  list(demand = demand, facilities = facilities)
}

# ---------------------------------------------------------------------------
# P-Median
# Pinned: objective = 3498.296111
# ---------------------------------------------------------------------------
test_that("p_median: pinned objective and invariants", {
  d <- make_solver_data(30, 10, seed = 101)

  result <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
  meta <- attr(result, "spopt")

  # Pinned objective
  expect_equal(meta$objective, 3498.296, tolerance = 1e-1)

  # Structural invariants
  expect_equal(meta$n_selected, 3)
  expect_equal(sum(result$facilities$.selected), 3)
  expect_true(all(result$demand$.facility %in% which(result$facilities$.selected)))
  expect_equal(nrow(result$demand), 30)
  expect_true(meta$mean_distance > 0)
})

# ---------------------------------------------------------------------------
# LSCP
# Pinned: objective = 2, n_selected = 2, coverage_pct = 100
# ---------------------------------------------------------------------------
test_that("lscp: pinned objective and coverage", {
  d <- make_solver_data(20, 8, seed = 202)

  result <- lscp(d$demand, d$facilities, service_radius = 5.0)
  meta <- attr(result, "spopt")

  # Pinned values
  expect_equal(meta$objective, 2)
  expect_equal(meta$n_selected, 2)
  expect_equal(meta$coverage_pct, 100)
  expect_equal(meta$uncoverable_demand, 0)

  # All demand should be covered
  expect_true(all(result$demand$.covered))
})

# ---------------------------------------------------------------------------
# MCLP
# Pinned: objective = 1667, coverage_pct = 66.0
# ---------------------------------------------------------------------------
test_that("mclp: pinned objective and coverage", {
  d <- make_solver_data(50, 10, seed = 303)

  result <- mclp(d$demand, d$facilities, service_radius = 3.0,
                 n_facilities = 3, weight_col = "pop")
  meta <- attr(result, "spopt")

  # Pinned values
  expect_equal(meta$objective, 1667, tolerance = 1e-1)
  expect_equal(meta$n_selected, 3)
  expect_equal(meta$coverage_pct, 66.0, tolerance = 0.5)
  expect_equal(meta$covered_weight, 1667, tolerance = 1e-1)
  expect_equal(meta$total_weight, 2527)
})

# ---------------------------------------------------------------------------
# P-Center (both methods)
# Pinned: max_distance = 4.555352 (both methods must agree)
# ---------------------------------------------------------------------------
test_that("p_center: pinned max_distance and method equivalence", {
  d <- make_solver_data(30, 10, seed = 404)

  result_bs <- p_center(d$demand, d$facilities, n_facilities = 3,
                        method = "binary_search")
  result_mip <- p_center(d$demand, d$facilities, n_facilities = 3,
                         method = "mip")
  meta_bs <- attr(result_bs, "spopt")
  meta_mip <- attr(result_mip, "spopt")

  # Pinned value
  expect_equal(meta_bs$max_distance, 4.5554, tolerance = 1e-2)

  # Method equivalence
  expect_equal(meta_bs$max_distance, meta_mip$max_distance, tolerance = 1e-4)

  # Structural invariants
  expect_equal(meta_bs$n_selected, 3)
  expect_equal(meta_mip$n_selected, 3)
  expect_true(all(result_bs$demand$.facility %in% which(result_bs$facilities$.selected)))
  expect_true(all(result_mip$demand$.facility %in% which(result_mip$facilities$.selected)))

  # Verify max_distance is actually the max assignment distance
  cost_mat <- as.matrix(sf::st_distance(d$demand, d$facilities))
  assigned_dists_bs <- sapply(seq_len(nrow(d$demand)), function(i) {
    cost_mat[i, result_bs$demand$.facility[i]]
  })
  expect_equal(max(assigned_dists_bs), meta_bs$max_distance, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# P-Dispersion
# Pinned: min_distance = 5.336305
# ---------------------------------------------------------------------------
test_that("p_dispersion: pinned min_distance and pairwise invariant", {
  set.seed(505)
  fac <- st_as_sf(
    data.frame(x = runif(15, 0, 10), y = runif(15, 0, 10)),
    coords = c("x", "y")
  )

  result <- p_dispersion(fac, n_facilities = 4)
  meta <- attr(result, "spopt")

  # Pinned value
  expect_equal(meta$min_distance, 5.3363, tolerance = 1e-2)
  expect_equal(meta$n_selected, 4)

  # Mathematical invariant: all pairwise distances >= min_distance
  selected_idx <- which(result$.selected)
  dist_mat <- as.matrix(sf::st_distance(result[selected_idx, ]))
  pairwise <- dist_mat[upper.tri(dist_mat)]
  expect_true(all(pairwise >= meta$min_distance - 1e-4))
})

# ---------------------------------------------------------------------------
# CFLP (fixed number of facilities)
# Pinned: objective = 4650.079, n_split = 0
# ---------------------------------------------------------------------------
test_that("cflp: pinned objective, allocation parity, and capacity invariants", {
  d <- make_solver_data(40, 10, seed = 606)
  cap_val <- max(sum(d$demand$pop) / 3, max(d$demand$pop))
  d$facilities$cap <- rep(cap_val, 10)

  result <- cflp(d$demand, d$facilities, n_facilities = 4,
                 weight_col = "pop", capacity_col = "cap")
  meta <- attr(result, "spopt")

  # Pinned values
  expect_equal(meta$objective, 4650.08, tolerance = 1)
  expect_equal(meta$n_selected, 4)
  expect_equal(meta$n_split_demand, 0)

  # Allocation matrix: must exist, correct dimensions, rows sum to 1
  alloc <- meta$allocation_matrix
  expect_true(!is.null(alloc))
  expect_equal(nrow(alloc), 40)
  expect_equal(ncol(alloc), 10)
  expect_equal(as.numeric(rowSums(alloc)), rep(1.0, 40), tolerance = 1e-6)

  # Allocation only to selected facilities
  selected <- which(result$facilities$.selected)
  non_selected <- setdiff(1:10, selected)
  expect_true(all(alloc[, non_selected] < 1e-6))

  # Capacity constraints: weighted allocation <= capacity for each selected facility
  weights <- d$demand$pop
  for (j in selected) {
    allocated_weight <- sum(weights * alloc[, j])
    expect_lte(allocated_weight, cap_val + 1e-6)
  }

  # Utilizations in [0, 1]
  selected_fac <- result$facilities[result$facilities$.selected, ]
  expect_true(all(selected_fac$.utilization >= -1e-6))
  expect_true(all(selected_fac$.utilization <= 1 + 1e-6))

  # Primary assignment consistent with allocation (only check rows
  # where a single facility has dominant allocation > 0.99)
  for (i in seq_len(40)) {
    max_alloc <- max(alloc[i, ])
    if (max_alloc > 0.99) {
      primary <- result$demand$.facility[i]
      expect_equal(primary, which.max(alloc[i, ]))
    }
  }
})

# ---------------------------------------------------------------------------
# CFLP (cost-based, no fixed n)
# Pinned: objective = 3409.637, n_selected = 7
# ---------------------------------------------------------------------------
test_that("cflp with facility costs: pinned objective and allocation", {
  set.seed(707)
  d <- make_solver_data(30, 8, seed = 707)
  d$facilities$cap <- rep(sum(d$demand$pop), 8)
  d$facilities$fcost <- runif(8, 10, 50)

  result <- cflp(d$demand, d$facilities, n_facilities = 0,
                 weight_col = "pop", capacity_col = "cap",
                 facility_cost_col = "fcost")
  meta <- attr(result, "spopt")

  # Pinned values
  expect_equal(meta$objective, 3409.64, tolerance = 1)
  expect_equal(meta$n_selected, 7)

  # Allocation rows sum to 1
  alloc <- meta$allocation_matrix
  expect_equal(as.numeric(rowSums(alloc)), rep(1.0, 30), tolerance = 1e-6)

  # All demands assigned
  expect_true(all(!is.na(result$demand$.facility)))
})
