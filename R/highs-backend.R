# Internal solver functions using the CRAN highs package.
# These replicate the MIP models from the Rust HiGHS backend.
# Each function takes the same arguments and returns the same list
# structure as the corresponding rust_* wrapper.

# ---------------------------------------------------------------------------
# LSCP: Location Set Covering Problem
# Minimize number of facilities to cover all demand within service_radius.
# ---------------------------------------------------------------------------
.solve_lscp <- function(cost_matrix, service_radius) {
  n_demand <- nrow(cost_matrix)
  n_fac <- ncol(cost_matrix)

  # Variables: y[j] binary (1..n_fac)
  L <- rep(1.0, n_fac)
  lower <- rep(0, n_fac)
  upper <- rep(1, n_fac)
  types <- rep("I", n_fac)

  # Coverage matrix: a[i,j] = 1 if cost_matrix[i,j] <= service_radius
  coverage <- cost_matrix <= service_radius

  # Constraints: for each coverable demand i, sum_j a[i,j]*y[j] >= 1
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

  sol <- res$primal_solution
  selected <- which(sol > 0.5)

  # Coverage: demand i is covered if any selected facility is within radius
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
  n_demand <- nrow(cost_matrix)
  n_fac <- ncol(cost_matrix)
  p <- n_facilities

  if (is.null(fixed_facilities)) fixed_facilities <- integer(0)

  # Variables: y[j] (1..n_fac), x[i,j] (n_fac+1 .. n_fac+n_demand*n_fac)
  n_vars <- n_fac + n_demand * n_fac

  # Objective: 0 for y, weights[i]*cost[i,j] for x
  L <- numeric(n_vars)
  for (i in seq_len(n_demand)) {
    for (j in seq_len(n_fac)) {
      L[n_fac + (i - 1L) * n_fac + j] <- weights[i] * cost_matrix[i, j]
    }
  }

  # Bounds
  lower <- rep(0, n_vars)
  upper <- rep(1, n_vars)

  # Fixed facilities: force y[j] = 1
  for (j in fixed_facilities) {
    lower[j] <- 1
  }

  # Types: y = integer, x = continuous
  types <- c(rep("I", n_fac), rep("C", n_demand * n_fac))

  # Constraints
  # 1: sum_j y[j] = p                           (1 row)
  # 2: sum_j x[i,j] = 1 for each i              (n_demand rows)
  # 3: x[i,j] - y[j] <= 0 for each i,j          (n_demand*n_fac rows)
  n_con <- 1L + n_demand + n_demand * n_fac

  # Pre-count nonzeros
  n_nz <- n_fac + n_demand * n_fac + 2L * n_demand * n_fac
  tri_i <- integer(n_nz)
  tri_j <- integer(n_nz)
  tri_x <- numeric(n_nz)
  k <- 0L

  # Constraint 1: sum_j y[j] = p
  for (j in seq_len(n_fac)) {
    k <- k + 1L
    tri_i[k] <- 1L; tri_j[k] <- j; tri_x[k] <- 1.0
  }

  # Constraint 2: sum_j x[i,j] = 1 for each i
  for (i in seq_len(n_demand)) {
    row <- 1L + i
    for (j in seq_len(n_fac)) {
      k <- k + 1L
      tri_i[k] <- row
      tri_j[k] <- n_fac + (i - 1L) * n_fac + j
      tri_x[k] <- 1.0
    }
  }

  # Constraint 3: x[i,j] - y[j] <= 0
  for (i in seq_len(n_demand)) {
    for (j in seq_len(n_fac)) {
      con_row <- 1L + n_demand + (i - 1L) * n_fac + j
      x_col <- n_fac + (i - 1L) * n_fac + j
      # x[i,j] coefficient
      k <- k + 1L
      tri_i[k] <- con_row; tri_j[k] <- x_col; tri_x[k] <- 1.0
      # -y[j] coefficient
      k <- k + 1L
      tri_i[k] <- con_row; tri_j[k] <- j; tri_x[k] <- -1.0
    }
  }

  A <- Matrix::sparseMatrix(i = tri_i, j = tri_j, x = tri_x,
                            dims = c(n_con, n_vars))

  lhs <- c(p, rep(1, n_demand), rep(-Inf, n_demand * n_fac))
  rhs <- c(p, rep(1, n_demand), rep(0, n_demand * n_fac))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  sol <- res$primal_solution

  # Extract selected facilities (1-based)
  y_sol <- sol[seq_len(n_fac)]
  selected <- which(y_sol > 0.5)

  # Extract assignments: for each demand, facility with max x value
  assignments <- integer(n_demand)
  for (i in seq_len(n_demand)) {
    x_vals <- sol[n_fac + (i - 1L) * n_fac + seq_len(n_fac)]
    assignments[i] <- which.max(x_vals)
  }

  # Mean weighted distance
  total_weighted_dist <- sum(vapply(seq_len(n_demand), function(i) {
    weights[i] * cost_matrix[i, assignments[i]]
  }, numeric(1)))
  total_weight <- sum(weights)
  mean_distance <- total_weighted_dist / total_weight

  list(
    selected = as.integer(selected),
    assignments = as.integer(assignments),
    n_selected = length(selected),
    objective = res$info$objective_function_value,
    mean_distance = mean_distance
  )
}
