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

#' Normalize a waypoints argument into a list of c(x, y) vectors
#'
#' Accepts:
#'   - NULL
#'   - sf / sfc POINT collection (length >= 0)
#'   - list of numeric c(x, y) vectors
#'   - list of sf/sfc single POINT features
#' @noRd
normalize_waypoints <- function(waypoints, raster_crs) {
  if (is.null(waypoints)) return(list())

  # sf or sfc of possibly-multiple points
  if (inherits(waypoints, "sf") || inherits(waypoints, "sfc")) {
    n <- if (inherits(waypoints, "sf")) nrow(waypoints) else length(waypoints)
    if (n == 0L) return(list())
    geom_types <- as.character(sf::st_geometry_type(waypoints, by_geometry = TRUE))
    if (any(geom_types != "POINT")) {
      stop("`waypoints` must contain only POINT geometries", call. = FALSE)
    }
    if (is.na(sf::st_crs(waypoints))) {
      stop(
        "`waypoints` has no CRS. Supply coordinates in the CRS of cost_surface ",
        "or set the CRS on the sf/sfc object.",
        call. = FALSE
      )
    }
    if (sf::st_crs(waypoints) != raster_crs) {
      waypoints <- sf::st_transform(waypoints, raster_crs)
    }
    coords <- sf::st_coordinates(waypoints)
    return(lapply(seq_len(n), function(i) c(coords[i, "X"], coords[i, "Y"])))
  }

  # list input
  if (is.list(waypoints)) {
    if (length(waypoints) == 0L) return(list())
    return(lapply(seq_along(waypoints), function(i) {
      normalize_corridor_point(waypoints[[i]],
                               sprintf("waypoint %d", i),
                               raster_crs)
    }))
  }

  stop(
    "`waypoints` must be NULL, an sf/sfc POINT collection, or a list of ",
    "points (numeric c(x, y) or single sf/sfc POINTs).",
    call. = FALSE
  )
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
# Multi-point (waypoints) helpers
# ---------------------------------------------------------------------------

#' Build a transient corridor graph directly from an already-prepared raster
#'
#' Used internally by route_corridor() when waypoints are supplied and the
#' input is a SpatRaster. Skips prepare_cost_surface() to avoid double
#' aggregation (the raster has already been prepared once by the caller).
#' @noRd
.build_graph_from_prepared <- function(prepared_surface, neighbours) {
  values <- terra::values(prepared_surface, mat = FALSE)
  values[is.na(values)] <- NaN
  n_rows <- terra::nrow(prepared_surface)
  n_cols <- terra::ncol(prepared_surface)
  res    <- terra::res(prepared_surface)

  ptr <- rust_corridor_build_graph(
    values, as.integer(n_rows), as.integer(n_cols),
    res[1], res[2], as.integer(neighbours)
  )
  info <- rust_corridor_graph_info(ptr)
  list(
    ptr              = ptr,
    cost_surface     = prepared_surface,
    n_edges          = info$n_edges,
    graph_build_time = info$build_time_ms / 1000
  )
}

#' Validate and snap a single point to its raster cell
#'
#' Returns 1-based cell index. Errors with `label` context.
#' @noRd
.validate_point_cell <- function(cost_surface, pt_xy, label) {
  cell <- terra::cellFromXY(cost_surface, matrix(pt_xy, ncol = 2))
  if (is.na(cell)) {
    stop(sprintf("%s falls outside the cost_surface extent", label), call. = FALSE)
  }
  val <- terra::extract(cost_surface, as.integer(cell))[1, 1]
  if (is.na(val)) {
    stop(
      sprintf("%s falls on an NA cell (impassable). ", label),
      "Ensure the point falls on a valid cell, or adjust the cost surface.",
      call. = FALSE
    )
  }
  as.integer(cell)
}

#' Solve a single segment via the Rust cached-graph entrypoint
#'
#' Returns the raw rust result augmented with cell centers and 1-based indices.
#' @noRd
.solve_segment_cached <- function(graph_ptr, cost_surface, origin_cell_1, dest_cell_1, method) {
  origin_cell_0 <- as.integer(origin_cell_1 - 1L)
  dest_cell_0   <- as.integer(dest_cell_1 - 1L)
  result <- rust_corridor_solve_cached(graph_ptr, origin_cell_0, dest_cell_0, method)

  cell_indices_1 <- result$path_cells + 1L
  origin_cell_center <- as.numeric(terra::xyFromCell(cost_surface, origin_cell_1))
  dest_cell_center   <- as.numeric(terra::xyFromCell(cost_surface, dest_cell_1))

  list(
    path_cells         = result$path_cells,
    cell_indices_1     = cell_indices_1,
    total_cost         = result$total_cost,
    n_edges_graph      = result$n_edges,
    solve_time         = result$solve_time_ms / 1000,
    graph_build_time   = result$graph_build_time_ms / 1000,
    origin_cell        = origin_cell_1,
    dest_cell          = dest_cell_1,
    origin_cell_center = origin_cell_center,
    dest_cell_center   = dest_cell_center
  )
}

#' Build sf output for multi-point (or segments-mode) corridor routing
#'
#' Handles both `output = "combined"` (single LINESTRING with per-segment
#' metadata in the `spopt` attr) and `output = "segments"` (N+1-row sf with
#' per-leg columns).
#'
#' For combined mode: interior waypoints are represented by their snapped
#' cell-center only; origin and destination use exact user-supplied coords.
#' This eliminates the spike artifact that would arise from chaining
#' single-segment outputs through off-cell-center waypoints.
#' @noRd
.build_corridor_sf <- function(segments, pts, labels, mode,
                               cost_surface, method, neighbours,
                               transient_build_time, n_waypoints_input,
                               original_waypoints_xy = list()) {
  raster_crs <- sf::st_crs(terra::crs(cost_surface))
  n_segments <- length(segments)
  n_pts      <- length(pts)

  # Per-segment geometry helpers
  seg_row_coords <- function(i) {
    seg  <- segments[[i]]
    o_xy <- pts[[i]]
    d_xy <- pts[[i + 1L]]
    sl_d <- sqrt(sum((o_xy - d_xy)^2))
    is_degen <- sl_d == 0 && seg$total_cost == 0
    if (is_degen) {
      m <- rbind(o_xy, o_xy)
      colnames(m) <- c("x", "y")
      list(coords = m, straight = 0, degen = TRUE)
    } else {
      seg_path_xy <- terra::xyFromCell(cost_surface, seg$cell_indices_1)
      coords <- rbind(
        matrix(o_xy, ncol = 2, dimnames = list(NULL, c("x", "y"))),
        seg_path_xy,
        matrix(d_xy, ncol = 2, dimnames = list(NULL, c("x", "y")))
      )
      list(coords = coords, straight = sl_d, degen = FALSE)
    }
  }

  seg_path_dist <- function(coords, degen) {
    if (degen) return(0)
    dx <- diff(coords[, 1])
    dy <- diff(coords[, 2])
    sum(sqrt(dx^2 + dy^2))
  }

  # Concatenate cells across segments, deduping at shared join cells.
  # Used by both combined-mode geometry and segments-mode aggregate n_cells.
  concat_dedup_cells <- function() {
    cells <- segments[[1]]$cell_indices_1
    if (n_segments > 1L) {
      for (i in 2:n_segments) {
        seg_cells <- segments[[i]]$cell_indices_1
        if (length(seg_cells) > 0L && length(cells) > 0L &&
            seg_cells[1] == cells[length(cells)]) {
          seg_cells <- seg_cells[-1]
        }
        cells <- c(cells, seg_cells)
      }
    }
    cells
  }

  if (mode == "segments") {
    rows <- lapply(seq_len(n_segments), function(i) {
      g <- seg_row_coords(i)
      seg_len <- seg_path_dist(g$coords, g$degen)
      sinu <- if (g$straight > 0) seg_len / g$straight else NA_real_
      df <- data.frame(
        segment            = as.integer(i),
        from_label         = labels[i],
        to_label           = labels[i + 1L],
        from_x             = pts[[i]][1],
        from_y             = pts[[i]][2],
        to_x               = pts[[i + 1L]][1],
        to_y               = pts[[i + 1L]][2],
        total_cost         = segments[[i]]$total_cost,
        n_cells            = length(segments[[i]]$cell_indices_1),
        path_dist          = seg_len,
        straight_line_dist = g$straight,
        sinuosity          = sinu,
        stringsAsFactors   = FALSE
      )
      df$geometry <- sf::st_sfc(sf::st_linestring(g$coords), crs = raster_crs)
      sf::st_as_sf(df)
    })
    combined <- do.call(rbind, rows)

    total_cost  <- sum(vapply(segments, function(s) s$total_cost, numeric(1)))
    # Deduped cell indices shared with combined mode: same semantic, same values
    cell_indices_dedup <- concat_dedup_cells()
    n_cells_tot <- length(cell_indices_dedup)
    solve_total <- sum(vapply(segments, function(s) s$solve_time, numeric(1)))
    n_rows_r    <- terra::nrow(cost_surface)
    n_cols_r    <- terra::ncol(cost_surface)
    res_r       <- terra::res(cost_surface)

    # Leg-baseline metrics match combined-mode definition
    leg_straight_line_dist_seg <- sum(vapply(seq_len(n_segments), function(i) {
      sqrt(sum((pts[[i + 1L]] - pts[[i]])^2))
    }, numeric(1)))
    path_dist_total <- sum(vapply(rows, function(r) r$path_dist, numeric(1)))
    leg_sinuosity_seg <- if (leg_straight_line_dist_seg > 0) {
      path_dist_total / leg_straight_line_dist_seg
    } else {
      NA_real_
    }

    # Waypoint bookkeeping: waypoints_input_xy preserves the full original
    # user-supplied list (length == n_waypoints_input). waypoint_cell_xy and
    # waypoint_cells are the effective routed joins (length == n_waypoints_effective).
    n_waypoints_effective <- max(n_segments - 1L, 0L)
    waypoints_input_xy <- if (length(original_waypoints_xy) > 0L) {
      do.call(rbind, original_waypoints_xy)
    } else {
      matrix(numeric(0), ncol = 2)
    }
    waypoint_cells <- if (n_waypoints_effective > 0L) {
      vapply(seq_len(n_waypoints_effective),
             function(k) segments[[k]]$dest_cell, integer(1))
    } else {
      integer(0)
    }
    waypoint_cell_xy <- if (length(waypoint_cells) > 0L) {
      unname(terra::xyFromCell(cost_surface, waypoint_cells))
    } else {
      matrix(numeric(0), ncol = 2)
    }

    attr(combined, "spopt") <- list(
      total_cost            = total_cost,
      n_cells               = as.integer(n_cells_tot),
      method                = method,
      neighbours            = as.integer(neighbours),
      n_cells_surface       = as.integer(n_rows_r) * as.integer(n_cols_r),
      n_edges_graph         = segments[[1]]$n_edges_graph,
      solve_time            = solve_total,
      graph_build_time      = transient_build_time,
      cell_indices          = cell_indices_dedup,
      origin_cell_center    = segments[[1]]$origin_cell_center,
      dest_cell_center      = segments[[n_segments]]$dest_cell_center,
      raster_dims           = c(n_rows_r, n_cols_r),
      cell_size             = c(res_r[1], res_r[2]),
      n_waypoints_input     = as.integer(n_waypoints_input),
      n_waypoints_effective = as.integer(n_waypoints_effective),
      n_segments_effective  = as.integer(n_segments),
      waypoints_input_xy    = waypoints_input_xy,
      waypoint_cell_xy      = waypoint_cell_xy,
      waypoint_cells        = waypoint_cells,
      segment_costs         = vapply(segments, function(s) s$total_cost, numeric(1)),
      segment_path_dists    = vapply(rows, function(r) r$path_dist, numeric(1)),
      segment_cells         = lapply(segments, function(s) s$cell_indices_1),
      segment_solve_times   = vapply(segments, function(s) s$solve_time, numeric(1)),
      leg_straight_line_dist = leg_straight_line_dist_seg,
      leg_sinuosity          = leg_sinuosity_seg
    )
    class(combined) <- c("spopt_corridor_segments", class(combined))
    return(combined)
  }

  # mode == "combined"
  all_cells <- concat_dedup_cells()

  origin_xy <- pts[[1]]
  dest_xy   <- pts[[n_pts]]
  straight_line_dist <- sqrt(sum((origin_xy - dest_xy)^2))
  is_degenerate <- straight_line_dist == 0 &&
    n_segments == 1L && segments[[1]]$total_cost == 0

  if (is_degenerate) {
    coords <- rbind(origin_xy, origin_xy)
    colnames(coords) <- c("x", "y")
  } else {
    path_xy <- terra::xyFromCell(cost_surface, all_cells)
    coords <- rbind(
      matrix(origin_xy, ncol = 2, dimnames = list(NULL, c("x", "y"))),
      path_xy,
      matrix(dest_xy, ncol = 2, dimnames = list(NULL, c("x", "y")))
    )
  }

  line <- sf::st_linestring(coords)
  sfc  <- sf::st_sfc(line, crs = raster_crs)
  output_sf <- sf::st_sf(geometry = sfc)

  total_cost <- sum(vapply(segments, function(s) s$total_cost, numeric(1)))
  n_cells    <- length(all_cells)
  path_dist  <- if (is_degenerate) 0 else as.numeric(sf::st_length(output_sf))
  sinuosity  <- if (straight_line_dist > 0) path_dist / straight_line_dist else NA_real_

  # Leg-baseline (sum of leg straight lines through ordered points)
  leg_straight_line_dist <- sum(vapply(seq_len(n_segments), function(i) {
    sqrt(sum((pts[[i + 1L]] - pts[[i]])^2))
  }, numeric(1)))
  leg_sinuosity <- if (leg_straight_line_dist > 0) {
    path_dist / leg_straight_line_dist
  } else {
    NA_real_
  }

  output_sf$total_cost         <- total_cost
  output_sf$n_cells            <- n_cells
  output_sf$straight_line_dist <- straight_line_dist
  output_sf$path_dist          <- path_dist
  output_sf$sinuosity          <- sinuosity

  # Expose leg_* columns only when waypoints were supplied (preserves legacy schema)
  if (n_waypoints_input > 0L) {
    output_sf$leg_straight_line_dist <- leg_straight_line_dist
    output_sf$leg_sinuosity          <- leg_sinuosity
  }

  n_rows_r <- terra::nrow(cost_surface)
  n_cols_r <- terra::ncol(cost_surface)
  res_r    <- terra::res(cost_surface)
  solve_total <- sum(vapply(segments, function(s) s$solve_time, numeric(1)))

  spopt_attr <- list(
    total_cost         = total_cost,
    n_cells            = n_cells,
    method             = method,
    neighbours         = as.integer(neighbours),
    n_cells_surface    = as.integer(n_rows_r) * as.integer(n_cols_r),
    n_edges_graph      = segments[[1]]$n_edges_graph,
    solve_time         = solve_total,
    graph_build_time   = transient_build_time,
    cell_indices       = all_cells,
    origin_cell_center = segments[[1]]$origin_cell_center,
    dest_cell_center   = segments[[n_segments]]$dest_cell_center,
    raster_dims        = c(n_rows_r, n_cols_r),
    cell_size          = c(res_r[1], res_r[2])
  )

  if (n_waypoints_input > 0L) {
    # Same semantic as segments mode: waypoints_input_xy is the full original
    # supplied list; waypoint_cell_xy / waypoint_cells are the effective joins.
    n_waypoints_effective <- max(n_segments - 1L, 0L)
    waypoints_input_xy <- if (length(original_waypoints_xy) > 0L) {
      do.call(rbind, original_waypoints_xy)
    } else {
      matrix(numeric(0), ncol = 2)
    }
    waypoint_cells <- if (n_waypoints_effective > 0L) {
      vapply(seq_len(n_waypoints_effective),
             function(k) segments[[k]]$dest_cell, integer(1))
    } else {
      integer(0)
    }
    waypoint_cell_xy <- if (length(waypoint_cells) > 0L) {
      unname(terra::xyFromCell(cost_surface, waypoint_cells))
    } else {
      matrix(numeric(0), ncol = 2)
    }

    segment_path_dists <- vapply(seq_len(n_segments), function(i) {
      g <- seg_row_coords(i)
      seg_path_dist(g$coords, g$degen)
    }, numeric(1))

    spopt_attr$n_waypoints_input     <- as.integer(n_waypoints_input)
    spopt_attr$n_waypoints_effective <- as.integer(n_waypoints_effective)
    spopt_attr$n_segments_effective  <- as.integer(n_segments)
    spopt_attr$waypoints_input_xy    <- waypoints_input_xy
    spopt_attr$waypoint_cell_xy      <- waypoint_cell_xy
    spopt_attr$waypoint_cells        <- waypoint_cells
    spopt_attr$segment_costs         <- vapply(segments, function(s) s$total_cost, numeric(1))
    spopt_attr$segment_path_dists    <- segment_path_dists
    spopt_attr$segment_cells         <- lapply(segments, function(s) s$cell_indices_1)
    spopt_attr$segment_solve_times   <- vapply(segments, function(s) s$solve_time, numeric(1))
    spopt_attr$leg_straight_line_dist <- leg_straight_line_dist
    spopt_attr$leg_sinuosity          <- leg_sinuosity
  }

  attr(output_sf, "spopt") <- spopt_attr
  class(output_sf) <- c("spopt_corridor", class(output_sf))
  output_sf
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Least-Cost Corridor Routing
#'
#' Find the minimum-cost path between an origin and destination across a
#' raster friction surface, optionally through an ordered sequence of
#' intermediate waypoints. The cost surface must be in a projected CRS
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
#' @param waypoints Optional intermediate points in the required visit order.
#'   \code{NULL} (default) for point-to-point routing. Otherwise an sf or sfc
#'   POINT collection, or a list of numeric c(x, y) vectors / single-point
#'   sf/sfc features. Each waypoint must fall on a non-NA cell. Interior
#'   waypoints are represented in the output geometry by their snapped cell
#'   center (origin and destination retain their exact user-supplied
#'   coordinates) to avoid join artifacts at off-cell-center waypoints.
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
#' @param output Character. Output shape when \code{waypoints} is non-NULL or
#'   \code{"segments"} is explicitly requested:
#'   - "combined" (default): a single-row sf LINESTRING for the full route,
#'     with per-segment breakdown in the \code{"spopt"} attribute.
#'   - "segments": an N+1-row sf, one LINESTRING per leg (class
#'     \code{spopt_corridor_segments}), with per-leg columns.
#'
#' @return An sf LINESTRING object with columns:
#'   \itemize{
#'     \item \code{total_cost}: accumulated traversal cost (friction * distance units)
#'     \item \code{n_cells}: number of cells in the path
#'     \item \code{straight_line_dist}: Euclidean distance origin to destination (CRS units)
#'     \item \code{path_dist}: actual path length (CRS units)
#'     \item \code{sinuosity}: path_dist / straight_line_dist (NA if origin == destination)
#'   }
#'
#'   When waypoints are supplied and \code{output = "combined"}, two
#'   additional columns report the leg-sum baseline:
#'   \code{leg_straight_line_dist} (sum of straight-line distances between
#'   consecutive ordered points) and \code{leg_sinuosity}
#'   (\code{path_dist / leg_straight_line_dist}). \code{sinuosity} retains
#'   its origin-to-destination semantics and will grow with the number of
#'   waypoints, as expected.
#'
#'   The returned linestring starts at the user-supplied origin coordinates,
#'   passes through the cell centers along the optimal path (including the
#'   snapped cell center at each interior waypoint), and ends at the
#'   user-supplied destination coordinates.
#'
#'   The "spopt" attribute contains metadata including \code{total_cost},
#'   \code{n_cells}, \code{method}, \code{neighbours}, \code{solve_time},
#'   \code{graph_build_time}, \code{cell_indices}, and grid dimensions.
#'   When any waypoints are supplied, the attribute also contains
#'   \code{n_waypoints_input} (count as originally supplied),
#'   \code{n_waypoints_effective} (count after eliding exact duplicates and
#'   same-cell collapses; may be less than \code{n_waypoints_input}),
#'   \code{n_segments_effective},
#'   \code{waypoints_input_xy} (\code{n_waypoints_input}-by-2 matrix of
#'   exact user-supplied coords -- always preserved, independent of
#'   elision), \code{waypoint_cell_xy} and \code{waypoint_cells}
#'   (length \code{n_waypoints_effective}, the snapped cells actually
#'   routed through), \code{segment_costs}, \code{segment_path_dists},
#'   \code{segment_cells}, and \code{segment_solve_times}.
#'   \code{cell_indices} in combined mode is the concatenated routed path
#'   after shared-join dedup -- the cells physically traversed, not a
#'   record of supplied waypoints.
#'
#'   With \code{output = "segments"}, the return is an sf with N+1 rows
#'   (one LINESTRING per leg) with columns \code{segment}, \code{from_label},
#'   \code{to_label}, per-leg coordinates, \code{total_cost}, \code{n_cells},
#'   \code{path_dist}, \code{straight_line_dist}, \code{sinuosity}. The
#'   object-level \code{"spopt"} attribute carries the same aggregate
#'   metadata as combined mode.
#'
#' @details
#' Ordered-waypoint routing returns the exact minimum-cost route through the
#' required waypoint sequence under the additive raster-cost model used here
#' -- shortest paths compose, so chaining segments between consecutive
#' required points is optimal for the ordered-via formulation. Two related
#' problems are \emph{not} solved by this API and belong to future
#' extensions: choosing the best visit order for unordered waypoints (a
#' TSP-style problem) and connecting multiple required cities with a
#' branched minimum-cost network on the cost surface (a Steiner tree
#' problem).
#'
#' @examples
#' \donttest{
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
#' # Ordered multi-city corridor
#' wp <- list(c(200000, 250000), c(350000, 350000))
#' via_path <- route_corridor(r, origin, dest, waypoints = wp)
#' via_segs <- route_corridor(r, origin, dest, waypoints = wp,
#'                            output = "segments")
#'
#' # Graph caching for multiple OD pairs or repeated waypoint routes
#' g <- corridor_graph(r, neighbours = 8L)
#' path1 <- route_corridor(g, origin, dest, method = "astar")
#' path2 <- route_corridor(g, origin, dest, waypoints = wp)
#' }
#'
#' @export
route_corridor <- function(cost_surface,
                           origin,
                           destination,
                           waypoints = NULL,
                           neighbours = 8L,
                           method = c("dijkstra", "bidirectional", "astar"),
                           resolution_factor = 1L,
                           output = c("combined", "segments")) {

  method <- match.arg(method)
  output <- match.arg(output)
  call_args <- match.call()

  # --- dispatch: cached graph vs raster ---
  if (inherits(cost_surface, "spopt_corridor_graph")) {
    return(.route_corridor_cached(
      cost_surface, origin, destination, method,
      waypoints = waypoints, output = output,
      call_args = call_args
    ))
  }

  # --- raster path ---
  prepared <- prepare_cost_surface(cost_surface, neighbours, resolution_factor)

  neighbours <- as.integer(neighbours)
  raster_crs <- sf::st_crs(terra::crs(prepared))

  origin_xy <- normalize_corridor_point(origin, "origin", raster_crs)
  dest_xy   <- normalize_corridor_point(destination, "destination", raster_crs)
  wp_list   <- normalize_waypoints(waypoints, raster_crs)
  n_wp_in   <- length(wp_list)

  # --- Legacy fast path: no waypoints + combined output ---
  # Preserves the exact pre-change single-segment code path / output.
  if (n_wp_in == 0L && output == "combined") {
    cells <- validate_corridor_cells(prepared, origin_xy, dest_xy)

    values <- terra::values(prepared, mat = FALSE)
    values[is.na(values)] <- NaN
    n_rows <- terra::nrow(prepared)
    n_cols <- terra::ncol(prepared)
    res    <- terra::res(prepared)

    origin_cell_0 <- as.integer(cells$origin_cell - 1L)
    dest_cell_0   <- as.integer(cells$dest_cell - 1L)

    result <- rust_corridor(
      values, as.integer(n_rows), as.integer(n_cols),
      res[1], res[2],
      origin_cell_0, dest_cell_0,
      neighbours, method
    )

    return(corridor_result_to_sf(result, prepared, origin_xy, dest_xy, method, neighbours))
  }

  # --- Multi-point / segments-mode path ---
  .route_corridor_multi(
    prepared_surface     = prepared,
    origin_xy            = origin_xy,
    dest_xy              = dest_xy,
    waypoint_xy_list     = wp_list,
    neighbours           = neighbours,
    method               = method,
    output               = output,
    graph_ptr            = NULL,
    external_graph_edges = NULL
  )
}

# ---------------------------------------------------------------------------
# Shared multi-point orchestration (used by raster and cached-graph paths)
# ---------------------------------------------------------------------------

#' @noRd
.route_corridor_multi <- function(prepared_surface,
                                  origin_xy,
                                  dest_xy,
                                  waypoint_xy_list,
                                  neighbours,
                                  method,
                                  output,
                                  graph_ptr,
                                  external_graph_edges) {
  n_wp_in <- length(waypoint_xy_list)

  # Build ordered points and labels
  pts <- c(list(origin_xy),
           waypoint_xy_list,
           list(dest_xy))
  labels <- c("origin",
              if (n_wp_in > 0L) sprintf("waypoint %d", seq_len(n_wp_in)),
              "destination")

  # Per-point validation + snap
  cell_ids <- integer(length(pts))
  for (i in seq_along(pts)) {
    cell_ids[i] <- .validate_point_cell(prepared_surface, pts[[i]], labels[i])
  }

  # Dedup exact-duplicate consecutive points (drop the later one)
  keep <- rep(TRUE, length(pts))
  for (i in seq.int(2L, length(pts))) {
    if (!keep[i - 1L]) next
    if (isTRUE(all.equal(pts[[i]], pts[[i - 1L]], tolerance = 0))) {
      warning(sprintf(
        "%s is identical to %s; dropping duplicate.",
        labels[i], labels[i - 1L]
      ), call. = FALSE)
      keep[i] <- FALSE
    }
  }
  pts      <- pts[keep]
  cell_ids <- cell_ids[keep]
  labels   <- labels[keep]

  if (length(pts) < 2L) {
    stop(
      "After removing duplicate points, fewer than 2 distinct points remain.",
      call. = FALSE
    )
  }

  # Collapse consecutive points that snap to the same cell. Origin (first)
  # and destination (last) are always preserved: the returned route must start
  # at the exact origin and end at the exact destination. Only intermediate
  # waypoints are elided. When a waypoint snaps to the same cell as the
  # destination, we drop the waypoint (not the destination).
  n_total <- length(pts)
  drop    <- rep(FALSE, n_total)
  last_kept <- 1L
  if (n_total >= 2L) {
    for (i in seq.int(2L, n_total)) {
      if (cell_ids[i] == cell_ids[last_kept]) {
        if (i == n_total && last_kept > 1L) {
          # Destination shares a cell with a waypoint: drop the waypoint
          warning(sprintf(
            "%s and %s snap to the same raster cell; eliding zero-length segment.",
            labels[last_kept], labels[i]
          ), call. = FALSE)
          drop[last_kept] <- TRUE
          last_kept <- i
        } else if (i < n_total) {
          # Interior waypoint shares a cell with prior kept point: drop it
          warning(sprintf(
            "%s and %s snap to the same raster cell; eliding zero-length segment.",
            labels[last_kept], labels[i]
          ), call. = FALSE)
          drop[i] <- TRUE
          # last_kept unchanged
        } else {
          # i == n_total and last_kept == 1 (origin): degenerate O == D,
          # keep both endpoints and let the solver produce a 0-cost trivial route
          last_kept <- i
        }
      } else {
        last_kept <- i
      }
    }
  }
  pts      <- pts[!drop]
  cell_ids <- cell_ids[!drop]
  labels   <- labels[!drop]

  n_segments <- length(pts) - 1L

  # --- Build solver context (transient graph when >= 2 segments and no graph supplied) ---
  transient_build_time <- 0
  n_edges_used <- external_graph_edges
  owns_graph <- FALSE

  if (is.null(graph_ptr)) {
    if (n_segments >= 2L) {
      g <- .build_graph_from_prepared(prepared_surface, neighbours)
      graph_ptr            <- g$ptr
      transient_build_time <- g$graph_build_time
      n_edges_used         <- g$n_edges
      owns_graph           <- TRUE
    }
  }

  # --- Solve each segment ---
  segments <- vector("list", n_segments)
  for (i in seq_len(n_segments)) {
    o_cell <- cell_ids[i]
    d_cell <- cell_ids[i + 1L]

    seg <- tryCatch(
      {
        if (!is.null(graph_ptr)) {
          .solve_segment_cached(graph_ptr, prepared_surface,
                                o_cell, d_cell, method)
        } else {
          # Single-segment raster solve (when waypoints empty, only segments mode reaches here)
          values <- terra::values(prepared_surface, mat = FALSE)
          values[is.na(values)] <- NaN
          n_rows <- terra::nrow(prepared_surface)
          n_cols <- terra::ncol(prepared_surface)
          res    <- terra::res(prepared_surface)
          result <- rust_corridor(
            values, as.integer(n_rows), as.integer(n_cols),
            res[1], res[2],
            as.integer(o_cell - 1L), as.integer(d_cell - 1L),
            neighbours, method
          )
          cell_indices_1 <- result$path_cells + 1L
          list(
            path_cells         = result$path_cells,
            cell_indices_1     = cell_indices_1,
            total_cost         = result$total_cost,
            n_edges_graph      = result$n_edges,
            solve_time         = result$solve_time_ms / 1000,
            graph_build_time   = result$graph_build_time_ms / 1000,
            origin_cell        = o_cell,
            dest_cell          = d_cell,
            origin_cell_center = as.numeric(terra::xyFromCell(prepared_surface, o_cell)),
            dest_cell_center   = as.numeric(terra::xyFromCell(prepared_surface, d_cell))
          )
        }
      },
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("No path exists", msg, fixed = TRUE)) {
          stop(sprintf("No path exists from %s to %s",
                       labels[i], labels[i + 1L]),
               call. = FALSE)
        }
        stop(e)
      }
    )

    # Single-segment degenerate (origin == destination cell): fabricate trivial result
    # The Rust solver should return a valid result here too; leave as-is.
    segments[[i]] <- seg
  }

  # --- Book-keeping: when only one segment & no transient graph, fill n_edges ---
  if (is.null(n_edges_used) && length(segments) >= 1L) {
    n_edges_used <- segments[[1]]$n_edges_graph
  }
  # Ensure segment n_edges_graph metadata reflects the shared graph (if transient)
  if (owns_graph || !is.null(external_graph_edges)) {
    for (i in seq_along(segments)) {
      segments[[i]]$n_edges_graph <- n_edges_used
    }
  }

  # When the raster single-segment branch ran (no transient graph built),
  # the Rust solver still constructed a graph internally. Report its time so
  # object-level graph_build_time reflects actual construction work.
  if (is.null(graph_ptr) && length(segments) == 1L) {
    transient_build_time <- segments[[1]]$graph_build_time
  }

  # --- Assemble sf output ---
  .build_corridor_sf(
    segments              = segments,
    pts                   = pts,
    labels                = labels,
    mode                  = output,
    cost_surface          = prepared_surface,
    method                = method,
    neighbours            = neighbours,
    transient_build_time  = transient_build_time,
    n_waypoints_input     = n_wp_in,
    original_waypoints_xy = waypoint_xy_list
  )
}

# ---------------------------------------------------------------------------
# Cached graph path
# ---------------------------------------------------------------------------

#' @noRd
.route_corridor_cached <- function(graph, origin, destination, method,
                                   waypoints = NULL,
                                   output = "combined",
                                   call_args) {
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
  wp_list   <- normalize_waypoints(waypoints, raster_crs)
  n_wp_in   <- length(wp_list)

  # Legacy fast path: no waypoints + combined output
  if (n_wp_in == 0L && output == "combined") {
    cells <- validate_corridor_cells(graph$cost_surface, origin_xy, dest_xy)

    origin_cell_0 <- as.integer(cells$origin_cell - 1L)
    dest_cell_0   <- as.integer(cells$dest_cell - 1L)

    result <- rust_corridor_solve_cached(
      graph$ptr, origin_cell_0, dest_cell_0, method
    )

    return(corridor_result_to_sf(
      result, graph$cost_surface, origin_xy, dest_xy, method, neighbours
    ))
  }

  # Multi-point or segments-mode path reusing the existing graph pointer
  .route_corridor_multi(
    prepared_surface     = graph$cost_surface,
    origin_xy            = origin_xy,
    dest_xy              = dest_xy,
    waypoint_xy_list     = wp_list,
    neighbours           = neighbours,
    method               = method,
    output               = output,
    graph_ptr            = graph$ptr,
    external_graph_edges = meta$n_edges
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
#' \donttest{
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
