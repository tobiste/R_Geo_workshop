#' Extract azimuths of line segments
#'
#' @param x sf object of type `"LINESTRING"` or `"MULTILINESTRING"`
#'
#' @return sf object of type `"POINT"` with the columns and entries of the first row of `x`
#'
#' @details
#' It is recommended to perform `line_azimuth()` on single line objects, i.e.
#' type `"LINESTRING"`, instead of `"MULTILINESTRING"`. This is because the azimuth
#' of the last point of a line will be calculated to the first point of the
#' next line otherwise. This will cause a warning message. For `MULTILINESTRING`
#' objects, use `lines_azimuths()`.
#'
#' @importFrom sf st_cast st_coordinates st_as_sf st_crs st_drop_geometry
#' @export
#'
#' @name line_azimuth
#' @examples
#' data("plates", package = "tectonicr")
#' subset(plates, pair == "af-eu") |>
#'   smoothr::densify() |>
#'   line_azimuth()
#'
#' line_azimuth(plates)
NULL

#' @rdname line_azimuth
#' @export
line_azimuth <- function(x) {
  if (nrow(x) > 1) warning("It is recommended to only use single line objects")
  if (any(sf::st_geometry_type(x) == "MULTILINESTRING")) warning("It is recommended to only use single line objects")
  mat <- x |>
    sf::st_cast("POINT") |>
    sf::st_coordinates()

  n <- nrow(mat)

  a <- numeric()
  for (i in 1:(n - 1)) {
    a[i] <- tectonicr::get_azimuth(mat[i, 2], mat[i, 1], mat[i + 1, 2], mat[i + 1, 1])
  }
  data.frame(
    x = mat[1:(n - 1), 1],
    y = mat[1:(n - 1), 2],
    azi = a
  ) |>
    cbind(sf::st_drop_geometry(x[1, ])) |>
    sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(x))
}

#' @rdname line_azimuth
#' @export
lines_azimuths <- function(x){
  for(i in 1:nrow(x)){
    ai <- line_azimuth(x[i, ])
    if(i == 1) {
      a = ai
    } else {
      a = rbind(a, ai)
    }
  }
  return(a)
}
