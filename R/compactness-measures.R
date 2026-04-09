#' @title Second areal moment (i.e., second moment of inertia)
#' @description
#' Computes the second moment of area (also known as the second moment of inertia)
#' for polygon geometries. This is a measure of how the area of a shape is
#' distributed relative to its centroid.
#' @details
#' The second moment of area is the sum of the inertia across the x and y axes:
#'
#' The inertia for the x axis is:
#' \deqn{I_x = \frac{1}{12}\sum_{i=1}^{N} (x_i y_{i+1} - x_{i+1}y_i) (x_i^2 + x_ix_{i+1} + x_{i+1}^2)}
#'
#' While the y axis is in a similar form:
#' \deqn{I_y = \frac{1}{12}\sum_{i=1}^{N} (x_i y_{i+1} - x_{i+1}y_i) (y_i^2 + y_iy_{i+1} + y_{i+1}^2)}
#'
#' where \eqn{x_i, y_i} is the current point and \eqn{x_{i+1}, y_{i+1}} is the
#' next point, and where \eqn{x_{n+1} = x_1, y_{n+1} = y_1}.
#'
#' For multipart polygons with holes, all parts are treated as separate
#' contributions to the overall centroid, which provides the same result as if
#' all parts with holes are separately computed, and then merged together using
#' the parallel axis theorem.
#' 
#' The code and documentation are adapted from the PySAL Python package
#' (Ray and Anselin, 2007). See Hally (1987) and Li et al. (2013) for additional details.
#' 
#' @param x An sf object, sfc geometry column, or sfg geometry.
#' @param project Logical. If the geometries have geodetic coordinates,
#' then they will be projected using an Albers Equal Area Conic projection centered on the data.
#'
#' @return Numeric vector of second areal moments
#'
#' @references
#' Hally, D. 1987. "The calculations of the moments of polygons." Canadian National
#' Defense Research and Development Technical Memorandum 87/209.
#' \url{https://apps.dtic.mil/sti/tr/pdf/ADA183444.pdf}
#' 
#' Li, W., Goodchild, M.F., and Church, R.L. 2013. 
#' "An Efficient Measure of Compactness for Two-Dimensional Shapes and Its Application in Regionalization Problems." 
#' International Journal of Geographical Information Science 27 (6): 1227–50. 
#' \doi{10.1080/13658816.2012.752093}.
#'
#' Rey, Sergio J., and Luc Anselin. 2007. "PySAL: A Python Library of Spatial Analytical Methods."
#' Review of Regional Studies 37 (1): 5–27.
#' \doi{10.52324/001c.8285}.
#'
#' @export
#' @seealso [nmi()], which computes the normalized moment of inertia.
#' @examples
#' library(sf)
#' poly <- st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE)))
#' second_areal_moment(poly)
second_areal_moment <- function(x, project = TRUE) {
  # Cast to sfc geometry column
  if (inherits(x, "sf")) {
    geoms <- sf::st_geometry(x)
  } else if (inherits(x, "sfc")) {
    geoms <- x
  } else if (inherits(x, "sfg")) {
    geoms <- sf::st_sfc(x)
  } else {
    stop("Input must be an sf object, sfc, or sfg geometry")
  }

  # Auto-project if geodetic
  if (project && !is.na(sf::st_crs(geoms)) && sf::st_is_longlat(geoms)) {
    warning("Coordinates are geodetic and will be projected.")
    geoms <- .auto_project(geoms)
  }
  
  # Process each geometry
  results <- vapply(geoms, .second_moment_single_geom, numeric(1))
  results
}

#' @title Second areal moment for a single geometry
#' @description Computes the second areal moment for a single geometry (POLYGON or MULTIPOLYGON)
#' @param geom A simple feature geometry
#' @return A single number, giving the second areal moment for the geometry
#' @keywords internal
.second_moment_single_geom <- function(geom) {
  # Get the overall centroid
  overall_centroid <- sf::st_centroid(sf::st_sfc(geom))
  
  # Cast to POLYGON if needed and extract all parts
  if (sf::st_geometry_type(geom) == "MULTIPOLYGON") {
    parts <- sf::st_cast(sf::st_sfc(geom), "POLYGON")
  } else {
    parts <- list(geom)
  }
  
  # Process each polygon part
  total_moa <- 0
  
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    coords_list <- sf::st_coordinates(part)
    
    # Get unique ring IDs (L1 identifies rings within a polygon)
    ring_ids <- unique(coords_list[, "L1"])
    
    for (ring_id in ring_ids) {
      ring_coords <- coords_list[coords_list[, "L1"] == ring_id, c("X", "Y"), drop = FALSE]
      
      # Remove duplicate last point if present
      if (nrow(ring_coords) > 1 && 
          all(ring_coords[1, ] == ring_coords[nrow(ring_coords), ])) {
        ring_coords <- ring_coords[-nrow(ring_coords), , drop = FALSE]
      }
      
      # Create a polygon from this ring for centroid calculation
      ring_poly <- sf::st_polygon(list(rbind(ring_coords, ring_coords[1, ])))
      ring_centroid <- sf::st_centroid(sf::st_sfc(ring_poly))
      
      # Compute moment of area for this ring (centered at ring centroid)
      ring_moa <- .second_moa_ring(ring_coords, ring_centroid)
      
      # Compute radius (distance from ring centroid to overall centroid)
      radius <- as.numeric(sf::st_distance(ring_centroid, overall_centroid))
      
      # Compute area of ring
      ring_area <- abs(sf::st_area(ring_poly))
      
      # Determine sign: exterior rings are positive, holes are negative
      # First ring (ring_id == 1) is always exterior
      sign <- ifelse(ring_id == 1, 1, -1)
      
      # Apply parallel axis theorem: I_total = I_centroidal + A*r^2
      total_moa <- total_moa + sign * (ring_moa + ring_area * radius^2)
    }
  }
  
  total_moa
}

#' @title Second areal moment for a single ring
#' @description Computes the second areal moment for a single ring
#' @param coords A matrix of coordinates for a ring
#' @param centroid Centroid coordinates, supplied as an object of class sf, sfc, or sfg
#' @return A single number, giving the second areal moment for a ring
#' @keywords internal
.second_moa_ring <- function(coords, centroid) {
  # Center coordinates at the centroid
  centroid_coords <- sf::st_coordinates(centroid)
  centered <- coords - matrix(rep(centroid_coords, each = nrow(coords)), 
                               ncol = 2, byrow = FALSE)
  
  # Compute moment of area using shoelace-based formula
  moi <- 0
  n <- nrow(centered)
  
  for (i in seq_len(n)) {
    # Current point
    x_tail <- centered[i, 1]
    y_tail <- centered[i, 2]
    
    # Next point (wrap around)
    next_i <- ifelse(i == n, 1, i + 1)
    x_head <- centered[next_i, 1]
    y_head <- centered[next_i, 2]
    
    # Shoelace component
    cross <- x_tail * y_head - x_head * y_tail
    
    # Second moment terms
    x_term <- x_head^2 + x_head * x_tail + x_tail^2
    y_term <- y_head^2 + y_head * y_tail + y_tail^2
    
    moi <- moi + cross * (x_term + y_term)
  }
  
  abs(moi / 12)
}

#' Automatically select an appropriate equal-area projection
#' @noRd
#' @keywords internal
.auto_project <- function(geom) {
  # Get bounding box
  bbox <- sf::st_bbox(geom)
  
  # Calculate center point
  center_lon <- mean(c(bbox["xmin"], bbox["xmax"]))
  center_lat <- mean(c(bbox["ymin"], bbox["ymax"]))
  
  # Use Albers Equal Area Conic centered on the data
  # This works well for most regional-scale analyses
  proj_string <- sprintf(
    "+proj=aea +lat_1=%.6f +lat_2=%.6f +lat_0=%.6f +lon_0=%.6f +datum=WGS84 +units=m",
    center_lat - 5,  # Standard parallel 1
    center_lat + 5,  # Standard parallel 2
    center_lat,      # Latitude of origin
    center_lon       # Central meridian
  )
  message(sprintf("Projecting coordinates: %s", proj_string))
  
  sf::st_transform(geom, proj_string)
}

#' @title Normalized moment of inertia (NMI)
#' @description
#' Computes the normalized moment of inertia (NMI), a compactness measure for polygon geometries.
#' The NMI ranges between 0 and 1, where 1 is the most compact shape (a circle) and 0 is an infinitely extending shape (Feng et al. 2022).
#' @details
#' The NMI is defined as follows, where \eqn{A} is the area of a geometry,
#' and \eqn{I} is the second moment of inertia (i.e., the second areal moment):
#' \deqn{\frac{A^2}{2 \pi I}}
#' See Li et al. (2013, 2014) for additional details.
#' @param x An sf object, sfc geometry column, or sfg geometry
#' @return Numeric vector of normalized moments of inertia.
#' @export
#' @references 
#' 
#' Feng, X., Rey, S., and Wei, R. (2022). "The max-p-compact-regions problem."
#' Transactions in GIS, 26, 717–734. \doi{10.1111/tgis.12874}.
#'
#' Li, W., Goodchild, M.F., and Church, R.L. 2013.
#' "An Efficient Measure of Compactness for Two-Dimensional Shapes and Its Application in Regionalization Problems."
#' International Journal of Geographical Information Science 27 (6): 1227–50.
#' \doi{10.1080/13658816.2012.752093}.
#'
#' Li, W., Church, R.L. and Goodchild, M.F. 2014.
#' "The p-Compact-regions Problem." Geogr Anal, 46: 250-273.
#' \doi{10.1111/gean.12038}.
#' @seealso [second_areal_moment()]
nmi <- function(x) {
  if (inherits(x, "sf")) {
    geoms <- sf::st_geometry(x)
  } else if (inherits(x, "sfc")) {
    geoms <- x
  } else if (inherits(x, "sfg")) {
    geoms <- sf::st_sfc(x)
  } else {
    stop("Input must be an sf object, sfc, or sfg geometry")
  }
  
  areas <- sf::st_area(geoms)
  sam <- second_areal_moment(geoms)
  
  return(as.numeric(areas)^2 / (2 * sam * pi))
}
