# ---------------------------------------------------------------------------
# Internal helpers for k-diverse corridor routing
# ---------------------------------------------------------------------------

#' Try to route a corridor, returning NULL on no-path failure
#' @noRd
.try_route_corridor <- function(cost_surface, origin_xy, dest_xy,
                                neighbours, method) {
  tryCatch(
    route_corridor(cost_surface, origin_xy, dest_xy,
                   neighbours = neighbours, method = method,
                   resolution_factor = 1L),
    error = function(e) {
      if (grepl("No path exists", conditionMessage(e), fixed = TRUE)) {
        return(NULL)
      }
      stop(e)
    }
  )
}

#' Recompute traversal cost on the original (unpenalized) surface
#'
#' Matches the Rust edge-weight formula: (friction_from + friction_to) / 2 * distance
#' @noRd
.recompute_cost_on_original <- function(cell_indices, original_surface) {
  vals <- terra::extract(original_surface, cell_indices)[, 1]
  xy   <- terra::xyFromCell(original_surface, cell_indices)
  dx   <- diff(xy[, 1])
  dy   <- diff(xy[, 2])
  dists <- sqrt(dx^2 + dy^2)
  mean_friction <- (vals[-length(vals)] + vals[-1]) / 2
  sum(mean_friction * dists)
}

#' Get 1-based cell indices within penalty_radius of a path
#' @noRd
.cells_within_radius <- function(path_sf, cost_surface, penalty_radius) {
  path_vect <- terra::vect(sf::st_geometry(path_sf))
  buf <- terra::buffer(path_vect, width = penalty_radius)
  terra::cells(cost_surface, buf)[, "cell"]
}

#' Mean distance from candidate path vertices to a reference linestring
#' @noRd
.mean_spacing <- function(candidate_sf, reference_sf) {
  candidate_coords <- sf::st_coordinates(candidate_sf)
  candidate_pts <- sf::st_as_sf(
    data.frame(X = candidate_coords[, "X"], Y = candidate_coords[, "Y"]),
    coords = c("X", "Y"), crs = sf::st_crs(candidate_sf)
  )
  dists <- sf::st_distance(candidate_pts, sf::st_geometry(reference_sf))
  mean(as.numeric(dists))
}

#' Fraction of candidate's cells within rank-1's buffer zone
#' @noRd
.pct_overlap <- function(candidate_cell_indices, rank1_buffer_cells) {
  sum(candidate_cell_indices %in% rank1_buffer_cells) / length(candidate_cell_indices)
}

# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' k-Diverse Corridor Routing
#'
#' Find k spatially distinct corridors between an origin and destination
#' using iterative penalty rerouting. Each iteration penalizes cells near
#' previously found paths, forcing the solver to find geographically
#' different alternatives.
#'
#' Corridors are scored on the original (unpenalized) cost surface for
#' fair comparison. Diversity metrics (\code{mean_spacing}, \code{pct_overlap})
#' are measured against the rank-1 (optimal) path.
#'
#' @param cost_surface A terra SpatRaster (single band). Same requirements as
#'   \code{\link{route_corridor}}. Must NOT be a \code{corridor_graph} object.
#' @param origin sf/sfc POINT or numeric c(x, y). Same as \code{\link{route_corridor}}.
#' @param destination sf/sfc POINT or numeric c(x, y).
#' @param k Integer >= 1. Number of diverse corridors to find. May return fewer
#'   if no feasible alternative exists.
#' @param penalty_radius Numeric. Buffer distance (CRS units) around each found
#'   path within which cells are penalized. If NULL (default), uses 5\% of the
#'   straight-line distance. Must be supplied explicitly when origin == destination.
#' @param penalty_factor Numeric > 1.0. Multiplier applied to cells within the
#'   penalty radius. Cumulative: a cell near paths 1 and 2 is multiplied by
#'   \code{penalty_factor^2}. Default 2.0.
#' @param min_spacing Numeric or NULL. Minimum average distance (CRS units)
#'   between a candidate and its nearest prior accepted path. Candidates below
#'   this threshold are discarded and the iteration retries with increased penalty.
#'   NULL (default) disables the check.
#' @param neighbours Integer. Cell connectivity: 4, 8 (default), or 16.
#' @param method Character. Routing algorithm: "dijkstra", "bidirectional", or
#'   "astar" (default, recommended for iterative use).
#' @param resolution_factor Integer >= 1. Aggregation factor applied once before
#'   routing begins.
#'
#' @return An sf object with up to k rows, each a LINESTRING, with columns:
#'   \code{corridor_rank}, \code{total_cost} (on original surface),
#'   \code{n_cells}, \code{path_dist}, \code{straight_line_dist},
#'   \code{sinuosity}, \code{mean_spacing} (vs rank 1, NA for rank 1),
#'   \code{pct_overlap} (vs rank 1, NA for rank 1).
#'
#'   Class is \code{c("spopt_k_corridors", "sf")}. The \code{"spopt"} attribute
#'   contains: \code{k_requested}, \code{k_found}, \code{penalty_radius},
#'   \code{penalty_factor}, \code{min_spacing}, \code{total_solve_time},
#'   \code{total_graph_build_time}, \code{total_retries}, \code{method},
#'   \code{neighbours}, \code{resolution_factor}.
#'
#' @examples
#' \dontrun{
#' library(terra); library(sf)
#' r <- rast(nrows = 200, ncols = 200, xmin = 0, xmax = 200000,
#'           ymin = 0, ymax = 200000, crs = "EPSG:32614")
#' values(r) <- runif(ncell(r), 0.5, 2.0)
#' result <- route_k_corridors(r, c(10000, 10000), c(190000, 190000), k = 3)
#' print(result)
#' plot(result)
#' }
#'
#' @export
route_k_corridors <- function(cost_surface,
                              origin,
                              destination,
                              k = 5L,
                              penalty_radius = NULL,
                              penalty_factor = 2.0,
                              min_spacing = NULL,
                              neighbours = 8L,
                              method = c("dijkstra", "bidirectional", "astar"),
                              resolution_factor = 1L) {

  # --- Parameter validation ---

  if (inherits(cost_surface, "spopt_corridor_graph")) {
    stop(
      "route_k_corridors() requires a SpatRaster, not a cached graph. ",
      "The cost surface is modified each iteration.",
      call. = FALSE
    )
  }

  cost_surface <- prepare_cost_surface(cost_surface, neighbours, resolution_factor)
  method <- match.arg(method)
  neighbours <- as.integer(neighbours)

  k <- as.integer(k)
  if (is.na(k) || k < 1L) {
    stop("`k` must be an integer >= 1", call. = FALSE)
  }

  penalty_factor <- as.numeric(penalty_factor)
  if (!is.finite(penalty_factor) || penalty_factor <= 1.0) {
    stop("`penalty_factor` must be finite and > 1.0", call. = FALSE)
  }

  if (!is.null(min_spacing)) {
    min_spacing <- as.numeric(min_spacing)
    if (!is.finite(min_spacing) || min_spacing < 0) {
      stop("`min_spacing` must be finite and >= 0", call. = FALSE)
    }
  }

  raster_crs <- sf::st_crs(terra::crs(cost_surface))
  origin_xy  <- normalize_corridor_point(origin, "origin", raster_crs)
  dest_xy    <- normalize_corridor_point(destination, "destination", raster_crs)
  validate_corridor_cells(cost_surface, origin_xy, dest_xy)

  straight_line_dist <- sqrt(sum((origin_xy - dest_xy)^2))

  if (straight_line_dist == 0) {
    stop(
      "origin and destination are identical. ",
      "k-diverse routing requires distinct endpoints.",
      call. = FALSE
    )
  }

  if (is.null(penalty_radius)) {
    penalty_radius <- 0.05 * straight_line_dist
  } else {
    penalty_radius <- as.numeric(penalty_radius)
    if (!is.finite(penalty_radius) || penalty_radius <= 0) {
      stop("`penalty_radius` must be finite and > 0", call. = FALSE)
    }
  }

  # --- Main loop ---

  original_surface <- terra::deepcopy(cost_surface)
  penalty_multiplier <- terra::rast(cost_surface)
  terra::values(penalty_multiplier) <- 1.0

  results <- list()
  buffer_cache <- list()
  rank_counter <- 0L
  total_retries <- 0L
  total_solve_time <- 0
  total_graph_build_time <- 0

  for (attempt in seq_len(k)) {
    if (rank_counter >= k) break

    working_surface <- original_surface * penalty_multiplier
    saved_multiplier <- terra::deepcopy(penalty_multiplier)
    accepted <- FALSE
    retries <- 0L

    while (retries <= 3L) {
      if (retries > 0L) {
        working_surface <- original_surface * penalty_multiplier
      }

      candidate <- .try_route_corridor(
        working_surface, origin_xy, dest_xy, neighbours, method
      )

      if (is.null(candidate)) break

      meta <- attr(candidate, "spopt")
      total_solve_time <- total_solve_time + meta$solve_time
      total_graph_build_time <- total_graph_build_time + meta$graph_build_time
      cell_indices <- meta$cell_indices

      # Recompute cost on original surface
      original_cost <- .recompute_cost_on_original(cell_indices, original_surface)

      # Diversity metrics
      mean_spacing_val <- NA_real_
      pct_overlap_val  <- NA_real_

      if (rank_counter > 0L) {
        mean_spacing_val <- .mean_spacing(candidate, results[[1]])
        pct_overlap_val  <- .pct_overlap(cell_indices, buffer_cache[[1]])

        # min_spacing check: vs nearest prior
        if (!is.null(min_spacing)) {
          nearest <- min(vapply(results, function(prior) {
            .mean_spacing(candidate, prior)
          }, numeric(1)))

          if (nearest < min_spacing) {
            penalty_multiplier <- terra::deepcopy(saved_multiplier)
            escalation <- 1.5 ^ (retries + 1L)
            for (j in seq_along(buffer_cache)) {
              prior_cells <- buffer_cache[[j]]
              penalty_multiplier[prior_cells] <-
                saved_multiplier[prior_cells] * escalation
            }
            retries <- retries + 1L
            total_retries <- total_retries + 1L
            next
          }
        }
      }

      # Accept
      rank_counter <- rank_counter + 1L

      path_buffer_cells <- .cells_within_radius(candidate, cost_surface, penalty_radius)
      buffer_cache[[rank_counter]] <- path_buffer_cells

      candidate$total_cost    <- original_cost
      candidate$mean_spacing  <- mean_spacing_val
      candidate$pct_overlap   <- pct_overlap_val
      candidate$corridor_rank <- rank_counter

      # Strip single-corridor class and metadata for clean assembly
      class(candidate) <- setdiff(class(candidate), "spopt_corridor")
      attr(candidate, "spopt") <- NULL
      results[[rank_counter]] <- candidate
      accepted <- TRUE
      break
    }

    if (!accepted) break

    # Penalize cells near this accepted path
    penalty_multiplier[buffer_cache[[rank_counter]]] <-
      penalty_multiplier[buffer_cache[[rank_counter]]] * penalty_factor
  }

  # --- Assemble output ---

  if (length(results) == 0L) {
    stop("No feasible corridor could be found.", call. = FALSE)
  }

  output <- do.call(rbind, results)
  col_order <- c("corridor_rank", "total_cost", "n_cells", "path_dist",
                 "straight_line_dist", "sinuosity", "mean_spacing", "pct_overlap")
  output <- output[, c(col_order, "geometry")]

  attr(output, "spopt") <- list(
    k_requested          = k,
    k_found              = rank_counter,
    penalty_radius       = penalty_radius,
    penalty_factor       = penalty_factor,
    min_spacing          = min_spacing,
    total_solve_time     = total_solve_time,
    total_graph_build_time = total_graph_build_time,
    total_retries        = total_retries,
    method               = method,
    neighbours           = neighbours,
    resolution_factor    = resolution_factor
  )

  class(output) <- c("spopt_k_corridors", class(output))
  output
}
