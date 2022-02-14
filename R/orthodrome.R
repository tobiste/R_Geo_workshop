#' @title Distance
#' @description This uses the "**haversine**" formula to calculate the great-circle
#' distance between two points – that is, the shortest distance over the earth’s
#' surface – giving an ‘as-the-crow-flies’ distance between the points
#' (ignoring any hills they fly over, of course!).
#' @param a lon, lat coordinate of point 1
#' @param b lon, lat coordinate of point 2
#' @return distance in km
#' @export
#' @examples
#' berlin <- c(52.517, 13.4)
#' tokyo <- c(35.7, 139.767)
#' haversine(berlin, tokyo)
haversine <- function(a, b) {
  r <- 6371.00887714 # km

  # convert deg into rad
  phi1 <- pi/180 * a[2]
  phi2 <- pi/180 * b[2]
  phi.delta <- (b[2] - a[2]) * (pi/180)
  lambda.delta <- (b[1] - a[1]) * (pi/180)

  a <- sin(phi.delta/2) * sin(phi.delta/2) +
    cos(phi1) * cos(phi2) *
    sin(lambda.delta/2) * sin(lambda.delta/2)

  c <- 2 * atan2(sqrt(a), sqrt(1-a))

  r * c
}

