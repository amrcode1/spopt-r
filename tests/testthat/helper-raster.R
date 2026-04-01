# Shared test helper for corridor tests
# Loaded automatically by testthat before any test file runs

make_test_raster <- function(nrow = 50, ncol = 50, vals = 1,
                             crs = "EPSG:32618") {
  r <- terra::rast(
    nrows = nrow, ncols = ncol,
    xmin = 0, xmax = ncol * 100,
    ymin = 0, ymax = nrow * 100,
    crs = crs
  )
  terra::values(r) <- vals
  r
}
