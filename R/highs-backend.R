# Internal solver functions using the CRAN highs package.
# These replicate the MIP models from the Rust HiGHS backend.
# Each function takes the same arguments and returns the same list
# structure as the corresponding rust_* wrapper.
#
# All constraint matrices are built with vectorized triplet construction
# (no nested R loops). The pattern for assignment models (p-median,
# p-center, cflp) uses shared index-generation helpers.

# ---------------------------------------------------------------------------
# Status check: stop with informative error if solver did not find optimal
# ---------------------------------------------------------------------------
.check_highs_status <- function(res, solver_name) {
  if (res$status_message != "Optimal") {
    stop(sprintf("%s solver returned non-optimal status: %s", solver_name,
                 res$status_message), call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# LSCP: Location Set Covering Problem
# Minimize number of facilities to cover all demand within service_radius.
# ---------------------------------------------------------------------------
.solve_lscp <- function(cost_matrix, service_radius) {
  n_demand <- nrow(cost_matrix)
  n_fac <- ncol(cost_matrix)

  L <- rep(1.0, n_fac)
  lower <- rep(0, n_fac)
  upper <- rep(1, n_fac)
  types <- rep("I", n_fac)

  coverage <- cost_matrix <= service_radius
  coverable <- which(rowSums(coverage) > 0)
  uncoverable <- n_demand - length(coverable)

  if (length(coverable) == 0) {
    return(list(
      selected = integer(0), n_selected = 0L, objective = 0,
      covered_demand = 0L, total_demand = n_demand,
      coverage_pct = 0, status = "Optimal", uncoverable_demand = uncoverable
    ))
  }

  A <- coverage[coverable, , drop = FALSE]
  storage.mode(A) <- "double"
  A <- Matrix::Matrix(A, sparse = TRUE)
  lhs <- rep(1, nrow(A))
  rhs <- rep(Inf, nrow(A))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "LSCP")

  sol <- res$primal_solution
  selected <- which(sol > 0.5)

  covered <- rowSums(cost_matrix[, selected, drop = FALSE] <= service_radius) > 0
  covered_demand <- sum(covered)
  coverage_pct <- (covered_demand / n_demand) * 100

  list(
    selected = as.integer(selected),
    n_selected = length(selected),
    objective = res$info$objective_function_value,
    covered_demand = as.integer(covered_demand),
    total_demand = as.integer(n_demand),
    coverage_pct = coverage_pct,
    status = res$status_message,
    uncoverable_demand = as.integer(uncoverable)
  )
}

# ---------------------------------------------------------------------------
# P-Median: minimize total weighted distance with exactly p facilities.
# ---------------------------------------------------------------------------
.solve_p_median <- function(cost_matrix, weights, n_facilities, fixed_facilities) {
  n_d <- nrow(cost_matrix)
  n_f <- ncol(cost_matrix)
  p <- n_facilities

  if (is.null(fixed_facilities)) fixed_facilities <- integer(0)

  # Variable layout: y[1..n_f], x[n_f+1 .. n_f+n_d*n_f] (row-major x[i,j])
  n_vars <- n_f + n_d * n_f

  # Objective: 0 for y, weights[i]*cost[i,j] for x (vectorized)
  L <- c(rep(0, n_f), as.vector(t(cost_matrix * weights)))

  # Bounds + fixed facilities
  lower <- rep(0, n_vars)
  upper <- rep(1, n_vars)
  lower[fixed_facilities] <- 1

  types <- c(rep("I", n_f), rep("C", n_d * n_f))

  # Sparse constraint matrix (fully vectorized triplet construction)
  A <- .build_assignment_constraints(n_d, n_f)

  lhs <- c(p, rep(1, n_d), rep(-Inf, n_d * n_f))
  rhs <- c(p, rep(1, n_d), rep(0, n_d * n_f))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "P-Median")

  sol <- res$primal_solution
  selected <- which(sol[seq_len(n_f)] > 0.5)

  # Extract assignments: reshape x portion into matrix, take column of max per row
  x_mat <- matrix(sol[(n_f + 1L):n_vars], nrow = n_d, ncol = n_f, byrow = TRUE)
  assignments <- max.col(x_mat, ties.method = "first")

  # Mean weighted distance (vectorized)
  assigned_costs <- cost_matrix[cbind(seq_len(n_d), assignments)]
  mean_distance <- sum(weights * assigned_costs) / sum(weights)

  list(
    selected = as.integer(selected),
    assignments = as.integer(assignments),
    n_selected = length(selected),
    objective = res$info$objective_function_value,
    mean_distance = mean_distance
  )
}

# ---------------------------------------------------------------------------
# Shared helper: build the standard assignment constraint matrix.
#
# Variable layout: y[1..n_f], x[n_f+1 .. n_f+n_d*n_f] (row-major)
# Constraints:
#   Row 1:              sum_j y[j] = p
#   Rows 2..n_d+1:      sum_j x[i,j] = 1 for each i
#   Rows n_d+2..end:     x[i,j] - y[j] <= 0 for each i,j
#
# Returns a sparse Matrix of dimensions (1 + n_d + n_d*n_f) x (n_f + n_d*n_f)
# ---------------------------------------------------------------------------
.build_assignment_constraints <- function(n_d, n_f) {
  n_vars <- n_f + n_d * n_f
  n_con <- 1L + n_d + n_d * n_f

  # Constraint 1: sum_j y[j] = p  (n_f nonzeros)
  c1_i <- rep(1L, n_f)
  c1_j <- seq_len(n_f)
  c1_x <- rep(1, n_f)

  # Constraint 2: sum_j x[i,j] = 1  (n_d * n_f nonzeros)
  # Row index: 1 + i, repeated n_f times per demand
  c2_i <- rep(1L + seq_len(n_d), each = n_f)
  # Column index: n_f + (i-1)*n_f + j
  c2_j <- n_f + seq_len(n_d * n_f)
  c2_x <- rep(1, n_d * n_f)

  # Constraint 3: x[i,j] - y[j] <= 0  (2 * n_d * n_f nonzeros)
  # Row index for constraint (i,j): 1 + n_d + (i-1)*n_f + j
  con3_rows <- 1L + n_d + seq_len(n_d * n_f)

  # x[i,j] coefficient (+1): same column as c2_j
  c3a_i <- con3_rows
  c3a_j <- n_f + seq_len(n_d * n_f)
  c3a_x <- rep(1, n_d * n_f)

  # -y[j] coefficient (-1): column j cycles 1..n_f for each demand
  c3b_i <- con3_rows
  c3b_j <- rep(seq_len(n_f), times = n_d)
  c3b_x <- rep(-1, n_d * n_f)

  Matrix::sparseMatrix(
    i = c(c1_i, c2_i, c3a_i, c3b_i),
    j = c(c1_j, c2_j, c3a_j, c3b_j),
    x = c(c1_x, c2_x, c3a_x, c3b_x),
    dims = c(n_con, n_vars)
  )
}
