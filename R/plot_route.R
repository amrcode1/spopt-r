#' Plot k-Diverse Corridors
#'
#' Renders all k corridors on a single plot. Rank 1 is drawn with full opacity
#' and maximum line width; higher ranks fade and thin proportionally.
#'
#' @param x An \code{spopt_k_corridors} object from \code{\link{route_k_corridors}}.
#' @param ... Additional arguments passed to the initial \code{plot()} call.
#'
#' @export
plot.spopt_k_corridors <- function(x, ...) {
  k <- nrow(x)
  base_col <- "royalblue"
  alphas <- seq(1.0, 0.3, length.out = max(k, 2))
  widths <- seq(3, 1, length.out = max(k, 2))

  plot(sf::st_geometry(x[1, ]),
       col = grDevices::adjustcolor(base_col, alpha.f = alphas[1]),
       lwd = widths[1], ...)
  if (k > 1) {
    for (i in 2:k) {
      plot(sf::st_geometry(x[i, ]),
           col = grDevices::adjustcolor(base_col, alpha.f = alphas[i]),
           lwd = widths[i], add = TRUE)
    }
  }
}
