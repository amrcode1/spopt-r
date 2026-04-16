#' Vehicle Routing Problem (VRP)
#'
#' Solves the Capacitated Vehicle Routing Problem (CVRP): find minimum-cost
#' routes for a fleet of vehicles, each with a capacity limit, to serve all
#' customers from a central depot. Uses Clarke-Wright savings heuristic for
#' construction with 2-opt, relocate, and swap local search improvement.
#'
#' @param locations An sf object representing all locations (depot + customers).
#' @param depot Integer. Row index of the depot in `locations` (1-based). Defaults to 1.
#' @param demand_col Character. Column name in `locations` containing the demand
#'   at each stop. The depot's demand is ignored (set to 0 internally).
#' @param vehicle_capacity Numeric. Maximum capacity per vehicle.
#' @param n_vehicles Integer or NULL. Maximum number of vehicles. If NULL
#'   (default), uses as many as needed to satisfy capacity constraints.
#' @param cost_matrix Optional. Pre-computed square cost/distance matrix (n x n).
#'   If NULL, computed from `locations` using `distance_metric`.
#' @param distance_metric Distance metric when computing from geometry:
#'   "euclidean" (default) or "manhattan".
#' @param method Algorithm: "2-opt" (default, savings + local search)
#'   or "savings" (Clarke-Wright construction only).
#' @param service_time Optional service/dwell time at each stop. Supply either
#'   a numeric vector of length `nrow(locations)` or a column name in
#'   `locations`. The depot's service time is forced to 0.
#' @param max_route_time Optional maximum total time per route (travel time +
#'   service time). Routes that would exceed this limit are split.
#' @param balance Optional balancing mode. Currently supports `"time"` to
#'   minimize the longest route's total time after cost optimization. This
#'   runs a post-optimization phase that may slightly increase total cost
#'   (bounded to a small percentage). Ignored when `method = "savings"`.
#' @param earliest Optional earliest arrival/service time at each stop. Supply
#'   either a numeric vector of length `nrow(locations)` or a column name in
#'   `locations`. Must be supplied together with `latest`.
#' @param latest Optional latest arrival/service time at each stop. Must be
#'   supplied together with `earliest`. A vehicle may arrive early and wait,
#'   but cannot begin service after this time.
#' @return An sf object (the input `locations`) with added columns:
#'
#'   \itemize{
#'     \item `.vehicle`: Vehicle assignment (1, 2, ...). Depot is 0.
#'     \item `.visit_order`: Visit sequence within each vehicle's route.
#'   }
#'   Metadata is stored in the "spopt" attribute, including `n_vehicles`,
#'   `total_cost`, `total_time`, per-vehicle `vehicle_costs`, `vehicle_times`,
#'   `vehicle_loads`, and `vehicle_stops`.
#'
#' @details
#' The CVRP extends the TSP to multiple vehicles with capacity constraints.
#' The solver uses a two-phase approach:
#'
#' 1. **Construction** (Clarke-Wright savings): Starts with one route per
#'    customer, then iteratively merges routes that produce the greatest
#'    distance savings while respecting capacity limits.
#' 2. **Improvement** (if `method = "2-opt"`): Applies intra-route 2-opt
#'    (reversing subtours), inter-route relocate (moving a customer between
#'    routes), and inter-route swap (exchanging customers between routes).
#'
#' @section Units and cost matrices:
#' The `cost_matrix` can contain travel times, distances, or any generalized
#' cost. The solver minimizes the sum of matrix values along each route.
#' When using time-based features (`service_time`, `max_route_time`,
#' `earliest`/`latest`, `balance = "time"`), the cost matrix should be in
#' time-compatible units (e.g., minutes). Otherwise the time-based outputs
#' and constraints will be dimensionally inconsistent.
#'
#' In the per-vehicle summary, **Cost** is the raw matrix objective (travel
#' only), while **Time** adds service duration and any waiting induced by
#' time windows. When no service time or windows are used, Time equals Cost
#' and is not shown separately.
#'
#' @section Use Cases:
#' \itemize{
#'   \item **Oilfield logistics**: Route vacuum trucks to well pads with
#'     fluid volume constraints
#'   \item **Delivery routing**: Multiple delivery vehicles with weight/volume limits
#'   \item **Service dispatch**: Assign and sequence jobs across a fleet
#'     of technicians or crews
#'   \item **Waste collection**: Route collection vehicles with load capacity
#' }
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Depot + 20 customers with demands
#' locations <- st_as_sf(
#'   data.frame(
#'     id = 1:21, x = runif(21), y = runif(21),
#'     demand = c(0, rpois(20, 10))
#'   ),
#'   coords = c("x", "y")
#' )
#'
#' result <- route_vrp(locations, depot = 1, demand_col = "demand", vehicle_capacity = 40)
#'
#' # How many vehicles needed?
#' attr(result, "spopt")$n_vehicles
#'
#' # Per-vehicle costs
#' attr(result, "spopt")$vehicle_costs
#' }
#'
#' @seealso [route_tsp()] for single-vehicle routing
#'
#' @export
route_vrp <- function(locations,
                depot = 1L,
                demand_col,
                vehicle_capacity,
                n_vehicles = NULL,
                cost_matrix = NULL,
                distance_metric = "euclidean",
                method = "2-opt",
                service_time = NULL,
                max_route_time = NULL,
                balance = NULL,
                earliest = NULL,
                latest = NULL) {
  resolve_route_vector <- function(value, arg_name, default = NULL) {
    if (is.null(value)) {
      return(default)
    }

    if (is.character(value) && length(value) == 1L) {
      if (!value %in% names(locations)) {
        stop(sprintf("Column '%s' supplied to `%s` was not found in `locations`",
                     value, arg_name), call. = FALSE)
      }
      value <- locations[[value]]
    }

    value <- as.numeric(value)
    if (length(value) != n) {
      stop(sprintf("`%s` must have length %d, got %d", arg_name, n, length(value)),
           call. = FALSE)
    }
    if (anyNA(value)) {
      stop(sprintf("`%s` must not contain NA values", arg_name), call. = FALSE)
    }

    value
  }

  # Input validation
  if (!inherits(locations, "sf")) {
    stop("`locations` must be an sf object", call. = FALSE)
  }

  n <- nrow(locations)

  if (n < 3) {
    stop("VRP requires at least 3 locations (1 depot + 2 customers)", call. = FALSE)
  }

  if (!demand_col %in% names(locations)) {
    stop(paste0("Demand column '", demand_col, "' not found in locations"), call. = FALSE)
  }

  depot <- as.integer(depot)
  if (depot < 1 || depot > n) {
    stop(sprintf("`depot` must be between 1 and %d", n), call. = FALSE)
  }

  method <- match.arg(method, c("2-opt", "savings"))

  demands <- as.numeric(locations[[demand_col]])
  demands[depot] <- 0  # depot has no demand

  if (vehicle_capacity <= 0) {
    stop("`vehicle_capacity` must be positive", call. = FALSE)
  }

  if (any(demands > vehicle_capacity, na.rm = TRUE)) {
    too_big <- which(demands > vehicle_capacity)
    stop(sprintf(
      "Demand at location(s) %s exceeds vehicle capacity (%.1f > %.1f)",
      paste(too_big, collapse = ", "),
      max(demands[too_big]),
      vehicle_capacity
    ), call. = FALSE)
  }

  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(locations, locations, type = distance_metric)
  }

  if (nrow(cost_matrix) != n || ncol(cost_matrix) != n) {
    stop(sprintf(
      "cost_matrix must be %d x %d, got %d x %d",
      n, n, nrow(cost_matrix), ncol(cost_matrix)
    ), call. = FALSE)
  }

  if (any(is.na(cost_matrix))) {
    n_na <- sum(is.na(cost_matrix))
    warning(sprintf(
      "cost_matrix contains %d NA values. Replacing with large value.", n_na
    ))
    max_cost <- max(cost_matrix, na.rm = TRUE)
    cost_matrix[is.na(cost_matrix)] <- max_cost * 100
  }

  # Resolve service_time
  service_vec <- resolve_route_vector(service_time, "service_time", default = NULL)
  if (!is.null(service_vec)) {
    if (any(service_vec < 0)) {
      stop("`service_time` must contain non-negative values", call. = FALSE)
    }
    service_vec[depot] <- 0  # depot has no service time
  }

  # Validate max_route_time
  if (!is.null(max_route_time)) {
    max_route_time <- as.numeric(max_route_time)
    if (max_route_time <= 0 || !is.finite(max_route_time)) {
      stop("`max_route_time` must be positive and finite", call. = FALSE)
    }

    # Pre-check: each customer must be individually reachable within time limit
    svc <- if (is.null(service_vec)) rep(0, n) else service_vec
    for (i in seq_len(n)) {
      if (i == depot) next
      round_trip <- cost_matrix[depot, i] + cost_matrix[i, depot] + svc[i]
      if (round_trip > max_route_time) {
        stop(sprintf(
          "Customer %d is unreachable within max_route_time=%.1f (depot->customer->depot + service = %.1f)",
          i, max_route_time, round_trip
        ), call. = FALSE)
      }
    }
  }

  # Resolve time windows
  if (xor(is.null(earliest), is.null(latest))) {
    stop("`earliest` and `latest` must be supplied together", call. = FALSE)
  }

  earliest_vec <- resolve_route_vector(earliest, "earliest", default = NULL)
  latest_vec <- resolve_route_vector(latest, "latest", default = NULL)

  if (!is.null(earliest_vec) && any(earliest_vec > latest_vec)) {
    stop("All `earliest` values must be less than or equal to `latest`", call. = FALSE)
  }

  use_windows <- !is.null(earliest_vec)

  # Pre-check: each customer individually feasible with windows
  if (use_windows) {
    svc <- if (is.null(service_vec)) rep(0, n) else service_vec
    depot_depart <- earliest_vec[depot]
    for (i in seq_len(n)) {
      if (i == depot) next
      arrive <- depot_depart + cost_matrix[depot, i]
      svc_start <- max(arrive, earliest_vec[i])
      if (svc_start > latest_vec[i] + 1e-10) {
        stop(sprintf(
          "Customer %d is infeasible: earliest service start (%.1f) exceeds latest window (%.1f)",
          i, svc_start, latest_vec[i]
        ), call. = FALSE)
      }
      if (!is.null(max_route_time)) {
        depart <- svc_start + svc[i]
        return_time <- depart + cost_matrix[i, depot]
        route_time <- return_time - depot_depart
        if (route_time > max_route_time + 1e-10) {
          stop(sprintf(
            "Customer %d is unreachable within max_route_time=%.1f (route time including waiting = %.1f)",
            i, max_route_time, route_time
          ), call. = FALSE)
        }
      }
    }
  }

  # Resolve balance mode
  balance_time <- FALSE
  if (!is.null(balance)) {
    balance <- match.arg(balance, c("time"))
    if (use_windows) {
      warning("balance is ignored when time windows are active", call. = FALSE)
    } else {
      balance_time <- TRUE
    }
  }

  if (!is.null(n_vehicles)) {
    n_vehicles <- as.integer(n_vehicles)
    if (n_vehicles < 1L) {
      stop("`n_vehicles` must be a positive integer", call. = FALSE)
    }
  }
  max_v <- n_vehicles

  start_time <- Sys.time()

  # Quick lower-bound feasibility check: ceiling(total_demand/capacity)
  # is the absolute minimum vehicles needed. More complex cases are
  # handled by the Rust solver; we validate the result afterward.
  if (!is.null(max_v)) {
    total_demand <- sum(demands)
    lb <- as.integer(ceiling(total_demand / vehicle_capacity))
    if (max_v < lb) {
      stop(sprintf(
        "Cannot satisfy n_vehicles=%d: total demand (%.1f) requires at least %d vehicles with capacity %.1f.",
        max_v, total_demand, lb, vehicle_capacity
      ), call. = FALSE)
    }
  }

  result <- rust_vrp(
    cost_matrix,
    depot - 1L,
    demands,
    vehicle_capacity,
    max_v,
    method,
    service_vec,
    max_route_time,
    balance_time,
    earliest_vec,
    latest_vec
  )

  end_time <- Sys.time()

  # Post-check: verify the solver actually achieved the requested fleet size
  # and didn't silently drop any customers
  if (!is.null(max_v) && result$n_vehicles > max_v) {
    if (!is.null(max_route_time)) {
      stop(sprintf(
        "Cannot satisfy n_vehicles=%d with max_route_time=%.1f: solver needed %d vehicles. Increase max_route_time, n_vehicles, or both.",
        max_v, max_route_time, result$n_vehicles
      ), call. = FALSE)
    } else {
      stop(sprintf(
        "Cannot satisfy n_vehicles=%d: solver needed %d vehicles to respect vehicle_capacity=%.1f. Increase vehicle_capacity or n_vehicles.",
        max_v, result$n_vehicles, vehicle_capacity
      ), call. = FALSE)
    }
  }

  # Check no customers were silently dropped (vehicle==0 means unassigned)
  unassigned <- which(result$vehicle == 0)
  # Depot is expected to be 0; filter it out
  unassigned <- setdiff(unassigned, depot)
  if (length(unassigned) > 0) {
    stop(sprintf(
      "Solver failed to assign %d customer(s) to vehicles. This indicates an infeasible n_vehicles=%d with vehicle_capacity=%.1f.",
      length(unassigned), max_v, vehicle_capacity
    ), call. = FALSE)
  }

  # Build output
  output <- locations
  # result$vehicle is 0-indexed for depot, 1-based for vehicles
  output$.vehicle <- result$vehicle
  output$.vehicle[depot] <- 0L  # mark depot
  output$.visit_order <- result$visit_order
  output$.visit_order[depot] <- 0L

  # Add arrival/departure columns when time windows are active
  if (use_windows) {
    output$.arrival_time <- result$arrival_times
    output$.arrival_time[depot] <- NA_real_
    output$.departure_time <- result$departure_times
    output$.departure_time[depot] <- NA_real_
  }

  metadata <- list(
    algorithm = "vrp",
    method = method,
    n_locations = n,
    depot = depot,
    n_vehicles = result$n_vehicles,
    total_cost = result$total_cost,
    total_time = result$total_time,
    initial_cost = result$initial_cost,
    improvement_pct = result$improvement_pct,
    iterations = result$iterations,
    vehicle_costs = result$vehicle_costs,
    vehicle_times = result$vehicle_times,
    vehicle_loads = result$vehicle_loads,
    vehicle_stops = result$vehicle_stops,
    vehicle_capacity = vehicle_capacity,
    max_route_time = max_route_time,
    has_service_time = !is.null(service_vec),
    has_time_windows = use_windows,
    balance = if (!is.null(balance) && method != "savings" && !use_windows) balance else NULL,
    balance_iterations = result$balance_iterations,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_vrp", "spopt_route", class(output))

  output
}
