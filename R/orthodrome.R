#' @title Orthodromic distance
#' @description The great-circle distance or orthodromic distance is the
#' shortest distance between two points on the surface of a sphere,
#' measured along the surface of the sphere (as opposed to a straight line
#' through the sphere's interior)
#' @param a lon, lat coordinate of point 1
#' @param b lon, lat coordinate of point 2
#' @return distance in km
#' @importFrom pracma acosd cosd sind deg2rad
#' @export
#' @examples
#' berlin <- c(52.517, 13.4)
#' tokyo <- c(35.7, 139.767)
#' orthodrome(berlin, tokyo)
orthodrome <- function(a, b) {
  require(pracma)
  delta <- pracma::acosd(
    pracma::sind(a[2]) * pracma::sind(b[2]) +
      pracma::cosd(a[2]) * pracma::cosd(b[2]) * pracma::cosd(b[1] - a[1])
  )

  d <- 6371.00887714 * pracma::deg2rad(delta)

  return(d)
}
