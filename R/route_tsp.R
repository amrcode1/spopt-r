#' Traveling Salesman Problem (TSP)
#'
#' Solves fixed-start routing problems over an `sf` layer. By default the route
#' is closed (start and end at the same location), but open routes and fixed
#' start/end paths are also supported. Optional time windows can be supplied in
#' the same units as the travel-time matrix.
#'
#' @param locations An sf object representing locations to visit.
#' @param depot Integer. Backward-compatible alias for `start` when `start` is
#'   not supplied. Defaults to 1.
#' @param start Integer. Row index of the route start in `locations` (1-based).
#'   If omitted, `depot` is used.
#' @param end Integer or NULL. Row index of the route end in `locations` (1-based).
#'   If omitted, defaults to `start` for a closed route. Set to NULL explicitly
#'   for an open route that may end at any stop.
#' @param cost_matrix Optional. Pre-computed square distance/cost matrix (n x n).
#'   If NULL, computed from `locations` using `distance_metric`.
#' @param distance_metric Distance metric when computing from geometry:
#'   "euclidean" (default) or "manhattan".
#' @param method Algorithm: "2-opt" (default, nearest-neighbor + local search)
#'   or "nn" (nearest-neighbor only).
#' @param earliest Optional earliest arrival/service times. Supply either a
#'   numeric vector of length `nrow(locations)` or a column name in `locations`.
#' @param latest Optional latest arrival/service times. Must be supplied
#'   together with `earliest`.
#' @param service_time Optional service duration at each stop. Supply either a
#'   numeric vector of length `nrow(locations)` or a column name in `locations`.
#'
#' @return An sf object (the input `locations`) with added columns:
#'
#'   \itemize{
#'     \item `.visit_order`: Visit sequence (1 = route start)
#'     \item `.tour_position`: Position in the returned tour/path
#'     \item `.arrival_time`: Arrival/service start time at each visited stop
#'     \item `.departure_time`: Departure time after service
#'   }
#'   Metadata is stored in the "spopt" attribute, including `total_cost`,
#'   `nn_cost`, `improvement_pct`, `tour`, `start`, `end`, `route_type`,
#'   and `solve_time`.
#'
#' @details
#' Supported route variants:
#'
#' 1. **Closed tour**: `start = 1, end = 1` (default)
#' 2. **Open route**: `start = 1, end = NULL`
#' 3. **Fixed path**: `start = 1, end = 5`
#'
#' Time windows use the same units as `cost_matrix`. When windows are supplied,
#' the solver constructs a feasible route and only accepts local-search moves
#' that preserve feasibility.
#'
#' @seealso [route_vrp()] for multi-vehicle routing, [distance_matrix()] for
#'   computing cost matrices
#'
#' @export
route_tsp <- function(locations,
                depot = 1L,
                start = NULL,
                end = NULL,
                cost_matrix = NULL,
                distance_metric = "euclidean",
                method = "2-opt",
                earliest = NULL,
                latest = NULL,
                service_time = NULL) {
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

  if (!inherits(locations, "sf")) {
    stop("`locations` must be an sf object", call. = FALSE)
  }

  n <- nrow(locations)
  if (n < 3) {
    stop("TSP requires at least 3 locations", call. = FALSE)
  }

  if (is.null(start)) {
    start <- depot
  }
  start <- as.integer(start)
  if (start < 1 || start > n) {
    stop(sprintf("`start` must be between 1 and %d", n), call. = FALSE)
  }

  if (missing(end)) {
    end <- start
  } else if (!is.null(end)) {
    end <- as.integer(end)
  }
  if (!is.null(end) && (end < 1 || end > n)) {
    stop(sprintf("`end` must be between 1 and %d", n), call. = FALSE)
  }

  method <- match.arg(method, c("2-opt", "nn"))

  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(locations, locations, type = distance_metric)
  }

  if (nrow(cost_matrix) != n || ncol(cost_matrix) != n) {
    stop(sprintf(
      "cost_matrix must be %d x %d (matching locations), got %d x %d",
      n, n, nrow(cost_matrix), ncol(cost_matrix)
    ), call. = FALSE)
  }

  if (any(is.na(cost_matrix))) {
    n_na <- sum(is.na(cost_matrix))
    warning(sprintf(
      "cost_matrix contains %d NA values. Replacing with large value.",
      n_na
    ))
    max_cost <- max(cost_matrix, na.rm = TRUE)
    cost_matrix[is.na(cost_matrix)] <- max_cost * 100
  }

  if (xor(is.null(earliest), is.null(latest))) {
    stop("`earliest` and `latest` must be supplied together", call. = FALSE)
  }

  earliest_vec <- resolve_route_vector(earliest, "earliest", default = NULL)
  latest_vec <- resolve_route_vector(latest, "latest", default = NULL)
  service_vec <- resolve_route_vector(service_time, "service_time", default = NULL)

  if (!is.null(service_vec) && any(service_vec < 0)) {
    stop("`service_time` must contain non-negative values", call. = FALSE)
  }
  if (!is.null(earliest_vec) && any(earliest_vec > latest_vec)) {
    stop("All `earliest` values must be less than or equal to `latest`", call. = FALSE)
  }

  start_time <- Sys.time()
  result <- rust_tsp(
    cost_matrix,
    start - 1L,
    method,
    earliest_vec,
    latest_vec,
    service_vec
  )
  end_time <- Sys.time()

  output <- locations
  output$.visit_order <- NA_integer_
  output$.arrival_time <- NA_real_
  output$.departure_time <- NA_real_

  tour <- result$tour
  for (i in seq_along(tour)) {
    idx <- tour[i]
    if (is.na(output$.visit_order[idx])) {
      output$.visit_order[idx] <- i
      output$.arrival_time[idx] <- result$arrival_time[i]
      output$.departure_time[idx] <- result$departure_time[i]
    }
  }

  output$.tour_position <- match(seq_len(n), tour)

  metadata <- list(
    algorithm = "tsp",
    method = method,
    n_locations = n,
    depot = if (!missing(depot) && start == depot && !is.null(end) && start == end) depot else NULL,
    start = start,
    end = end,
    route_type = if (is.null(end)) "open" else if (start == end) "closed" else "path",
    has_time_windows = !is.null(earliest_vec) || !is.null(service_vec),
    tour = tour,
    total_cost = result$total_cost,
    nn_cost = result$nn_cost,
    improvement_pct = result$improvement_pct,
    iterations = result$iterations,
    arrival_time = result$arrival_time,
    departure_time = result$departure_time,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_tsp", "spopt_route", class(output))

  output
}
