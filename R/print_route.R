#' @export
print.spopt_tsp <- function(x, ...) {
  if (!inherits(x, "sf")) return(NextMethod())
  meta <- attr(x, "spopt")
  cat(sprintf("TSP route: %d locations, %s tour\n",
              meta$n_locations, meta$route_type))
  cat(sprintf("  Method: %s | Cost: %.2f | Improvement: %.1f%% over NN\n",
              meta$method, meta$total_cost, meta$improvement_pct))
  if (isTRUE(meta$has_time_windows)) {
    cat("  Time windows: active\n")
  }
  cat(sprintf("  Solve time: %.3fs\n", meta$solve_time))
  invisible(x)
}

#' @export
print.spopt_vrp <- function(x, ...) {
  if (!inherits(x, "sf")) return(NextMethod())
  meta <- attr(x, "spopt")
  cat(sprintf("VRP routes: %d locations, %d vehicles (depot: %d)\n",
              meta$n_locations, meta$n_vehicles, meta$depot))
  cat(sprintf("  Method: %s | Total cost: %.1f | Improvement: %.1f%%\n",
              meta$method, meta$total_cost, meta$improvement_pct))

  constraints <- character(0)
  if (!is.null(meta$vehicle_capacity)) {
    constraints <- c(constraints, sprintf("Capacity: %.0f", meta$vehicle_capacity))
  }
  if (!is.null(meta$max_route_time)) {
    constraints <- c(constraints, sprintf("Max route time: %.0f", meta$max_route_time))
  }
  if (isTRUE(meta$has_time_windows)) {
    constraints <- c(constraints, "Time windows: active")
  }
  if (!is.null(meta$balance)) {
    bal_str <- sprintf("Balance: %s (%d move%s)",
                       meta$balance, meta$balance_iterations,
                       if (meta$balance_iterations == 1) "" else "s")
    constraints <- c(constraints, bal_str)
  }
  if (length(constraints) > 0) {
    cat(sprintf("  %s\n", paste(constraints, collapse = " | ")))
  }

  cat(sprintf("  Solve time: %.3fs\n", meta$solve_time))
  invisible(x)
}

#' @export
summary.spopt_tsp <- function(object, ...) {
  meta <- attr(object, "spopt")
  print.spopt_tsp(object)

  cat("\nTour sequence:\n")
  tour <- meta$tour
  cat(sprintf("  %s\n", paste(tour, collapse = " -> ")))

  if (isTRUE(meta$has_time_windows) &&
      !is.null(meta$arrival_time) &&
      length(meta$arrival_time) > 0) {
    cat("\nSchedule:\n")
    # Use tour order from metadata, exclude depot return (last element of closed tour)
    n_sched <- length(tour)
    if (meta$route_type == "closed" && n_sched > 1 && tour[n_sched] == tour[1]) {
      n_sched <- n_sched - 1
    }
    schedule <- data.frame(
      Stop = tour[seq_len(n_sched)],
      Arrival = round(meta$arrival_time[seq_len(n_sched)], 2),
      Departure = round(meta$departure_time[seq_len(n_sched)], 2)
    )
    print(schedule, row.names = FALSE)
  }

  invisible(object)
}

#' @export
summary.spopt_vrp <- function(object, ...) {
  meta <- attr(object, "spopt")
  print.spopt_vrp(object)

  show_time <- !isTRUE(all.equal(meta$vehicle_times, meta$vehicle_costs))

  cat("\nPer-vehicle summary:\n")
  tbl <- data.frame(
    Vehicle = seq_len(meta$n_vehicles),
    Stops = meta$vehicle_stops,
    Load = meta$vehicle_loads,
    Cost = meta$vehicle_costs
  )
  if (show_time) {
    tbl$Time <- meta$vehicle_times
  }
  print(tbl, row.names = FALSE)

  if (show_time) {
    cat("\n  Cost = matrix objective (travel only); Time = travel + service + waiting\n")
  }

  invisible(object)
}

#' @export
print.spopt_corridor <- function(x, ...) {
  meta <- attr(x, "spopt")
  n_wp <- meta$n_waypoints_input
  has_wp <- !is.null(n_wp) && n_wp > 0L

  header <- if (has_wp) "Least-cost corridor (via waypoints)" else "Least-cost corridor"
  cat(header, "\n", sep = "")
  cat("  Method:", meta$method, "\n")
  cat("  Total cost:", round(meta$total_cost, 2), "\n")
  cat("  Path distance:", round(x$path_dist, 0), "\n")
  cat("  Cells traversed:", meta$n_cells, "\n")
  cat("  Sinuosity:", round(x$sinuosity, 3), "\n")
  if (has_wp && !is.null(x$leg_sinuosity)) {
    cat("  Leg sinuosity:", round(x$leg_sinuosity, 3), "\n")
  }
  cat("  Solve time:", round(meta$solve_time, 3), "s\n")

  if (has_wp) {
    n_eff <- meta$n_waypoints_effective %||% n_wp
    if (!is.null(n_eff) && n_eff != n_wp) {
      cat(sprintf("  Waypoints: %d supplied, %d effective (duplicates/same-cell elided)\n",
                  n_wp, n_eff))
    } else {
      cat(sprintf("  Waypoints: %d\n", n_wp))
    }
    seg_costs <- meta$segment_costs
    seg_dists <- meta$segment_path_dists
    if (length(seg_costs) > 0L) {
      cat(sprintf("  %-10s  %12s  %12s\n", "Segment", "Cost", "Distance"))
      for (i in seq_along(seg_costs)) {
        cat(sprintf("  %-10s  %12s  %12.0f\n",
                    sprintf("%d", i),
                    format(round(seg_costs[i], 0), big.mark = ","),
                    seg_dists[i]))
      }
    }
  }

  invisible(x)
}

#' @export
print.spopt_corridor_segments <- function(x, ...) {
  meta <- attr(x, "spopt")
  cat("Least-cost corridor: per-segment output\n")
  cat("  Method:", meta$method, "\n")
  cat(sprintf("  Segments: %d (waypoints: %d supplied, %d effective)\n",
              meta$n_segments_effective,
              meta$n_waypoints_input,
              meta$n_waypoints_effective))
  cat("  Total cost:", round(meta$total_cost, 2), "\n")
  cat("  Solve time:", round(meta$solve_time, 3), "s",
      sprintf("(graph build: %.3fs)\n", meta$graph_build_time))

  cat(sprintf("\n  %-4s  %-14s  %-14s  %12s  %12s  %9s\n",
              "#", "From", "To", "Cost", "Distance", "Sinuosity"))
  for (i in seq_len(nrow(x))) {
    row <- x[i, , drop = FALSE]
    sinu_str <- if (is.na(row$sinuosity)) "-" else sprintf("%.3f", row$sinuosity)
    cat(sprintf("  %-4d  %-14s  %-14s  %12s  %12.0f  %9s\n",
                row$segment,
                substr(row$from_label, 1, 14),
                substr(row$to_label, 1, 14),
                format(round(row$total_cost, 0), big.mark = ","),
                row$path_dist,
                sinu_str))
  }
  invisible(x)
}

# Null-coalescing helper (local utility for print method)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.spopt_corridor_graph <- function(x, ...) {
  meta <- attr(x, "spopt")
  cat("Corridor graph\n")
  cat(sprintf("  Grid: %d x %d (%s cells)\n",
      meta$n_rows, meta$n_cols,
      format(meta$n_cells_surface, big.mark = ",")))
  cat(sprintf("  Cell size: %.1f x %.1f\n",
      meta$cell_size[1], meta$cell_size[2]))
  cat(sprintf("  Neighbours: %d (%s edges)\n",
      meta$neighbours,
      format(meta$n_edges, big.mark = ",")))
  cat(sprintf("  Build time: %.3fs | Graph storage: ~%.1f MB\n",
      meta$graph_build_time, meta$graph_storage_mb))
  invisible(x)
}

#' @export
print.spopt_k_corridors <- function(x, ...) {
  meta <- attr(x, "spopt")
  cat("k-Diverse Corridor Routing (spopt)\n")
  cat(sprintf("  Corridors found: %d of %d requested\n",
      meta$k_found, meta$k_requested))
  cat(sprintf("  Penalty: %.1fx within %.1f of each prior path\n",
      meta$penalty_factor, meta$penalty_radius))
  routing_time <- meta$total_solve_time + meta$total_graph_build_time
  cat(sprintf("  Routing time: %.3fs (solve: %.3fs, graph build: %.3fs)\n\n",
      routing_time, meta$total_solve_time, meta$total_graph_build_time))

  cat(sprintf("  %-15s  %10s  %10s  %9s  %10s  %7s\n",
      "", "Cost", "Distance", "Sinuosity", "Spacing", "Overlap"))
  for (i in seq_len(nrow(x))) {
    row <- x[i, , drop = FALSE]
    label <- if (row$alternative == 1L) "Optimal" else sprintf("Alternative %d", row$alternative - 1L)
    spacing_str <- if (is.na(row$mean_spacing)) "-" else sprintf("%.1f", row$mean_spacing)
    overlap_str <- if (is.na(row$pct_overlap)) "-" else sprintf("%.1f%%", row$pct_overlap * 100)
    cat(sprintf("  %-15s  %10s  %10.0f  %9.3f  %10s  %7s\n",
        label,
        format(round(row$total_cost, 0), big.mark = ","),
        row$path_dist,
        row$sinuosity,
        spacing_str,
        overlap_str))
  }
  invisible(x)
}
