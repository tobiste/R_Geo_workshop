#' @title swathR v1.0.1
#' @author V. Haburaj
#' @description Calculate swath-profile values perpendicular to a straight baseline. The
#' baseline is generated between two user-defined points (X|Y), see argument
#' \code{coords}. The distance between samples and the number of samples can be
#' specified, see arguments \code{k} and \code{dist}. Values of the swath-profile are
#' extracted from a given raster file, see argument \code{raster}. CRS of raster
#' and points have to be the same.
#' @param coords matrix(ncol=2, nrow=2) with x and y coordinates of beginning and
#' end point of the baseline; each point in one row
#' \describe{
#'   \item{column 1}{xcoordinates}
#'   \item{column 2}{ycoordinates}
#' }
#' @param raster raster file (loaded with [raster::raster()])
#' @param k integer; number of lines on each side of the baseline
#' @param dist numeric; distance between lines
#' @param crs string; CRS
#' @param method string; method for extraction of raw data, see
#' [raster::extract()]: default value: "bilinear"
#' @importFrom sp SpatialPoints CRS
#' @importFrom raster extract spLines
#' @export
swathR <- function(coords, raster, k = 1, dist = .2, crs = "+proj=longlat", method = "bilinear") {
  message("Initializing ...")
  # set default method:

  # create SpatialPoints from coords:
  spt <- sp::SpatialPoints(coords, proj4string = sp::CRS(crs))
  # get slope of baseline:
  m <- (raster::ymin(spt[1]) - raster::ymin(spt[2])) / (raster::xmin(spt[1]) - raster::xmin(spt[2]))
  # get slope of normal function:
  m1 <- -(1 / m)
  # get slope-angle from slope:
  alpha <- atan(m)
  alpha1 <- atan(m1)
  # get deltax and deltay from pythagoras:
  if ((alpha * 180) / pi < 90 & (alpha * 180) / pi > 270) {
    deltax <- cos(alpha1) * dist
  } else {
    deltax <- cos(alpha1) * dist * -1
  }
  if ((alpha * 180) / pi > 0 & (alpha * 180) / pi < 180) {
    deltay <- sqrt(dist**2 - deltax**2)
  } else {
    deltay <- sqrt(dist**2 - deltax**2) * -1
  }
  # create empty matrix:
  swath <- matrix(nrow = k * 2 + 1, ncol = 8)
  colnames(swath) <- c(
    "distance", "mean", "median",
    "std.dev.", "min", "max", "quantile(25)", "quantile(75)"
  )
  # list for spatial lines:
  allLines <- list()
  # add baseline:
  allLines[[k + 1]] <- raster::spLines(coords, crs = crs)
  # set distance for baseline:
  swath[k + 1, 1] <- 0
  # generate k lines parallel to baseline:
  for (n in 1:k) {
    # BELOW BASELINE:
    # new points
    cn <- matrix(nrow = 2, ncol = 2)
    cn[1, ] <- cbind(coords[1, 1] - (deltax * n), coords[1, 2] - (deltay * n))
    cn[2, ] <- cbind(coords[2, 1] - (deltax * n), coords[2, 2] - (deltay * n))
    # line between points:
    allLines[[k + 1 - n]] <- raster::spLines(cn, crs = crs)
    # distance value:
    swath[k + 1 - n, 1] <- -1 * n * dist
    # ABOVE BASELINE:
    # new points
    cn <- matrix(nrow = 2, ncol = 2)
    cn[1, ] <- cbind(coords[1, 1] + (deltax * n), coords[1, 2] + (deltay * n))
    cn[2, ] <- cbind(coords[2, 1] + (deltax * n), coords[2, 2] + (deltay * n))
    # line between points:
    allLines[[k + n + 1]] <- raster::spLines(cn, crs = crs)
    # distance value:
    swath[k + n + 1, 1] <- n * dist
  }
  # get raw data:
  message("Extracting raw data (this may take some time) ...")
  raw.data <- sapply(allLines, FUN = function(x) {
    raster::extract(raster, x, method = method)
  })
  # generalise data:
  message("Generalising data ...")
  swath[, 2] <- sapply(raw.data, function(x) {
    mean(x, na.rm = T)
  })
  swath[, 3] <- sapply(raw.data, function(x) {
    median(x, na.rm = T)
  })
  swath[, 4] <- sapply(raw.data, function(x) {
    sd(x, na.rm = T)
  })
  swath[, 5] <- sapply(raw.data, function(x) {
    min(x, na.rm = T)
  })
  swath[, 6] <- sapply(raw.data, function(x) {
    max(x, na.rm = T)
  })
  swath[, 7] <- sapply(raw.data, function(x) {
    quantile(x, na.rm = T)[2]
  })
  swath[, 8] <- sapply(raw.data, function(x) {
    quantile(x, na.rm = T)[4]
  })
  # return results:
  results <- list(swath = swath, data = raw.data, lines = allLines)
  message("Operation finished successfully!")
  message('Structure of results (list): "swath": swath profile data (matrix, numeric), "data": raw data (list, numeric), "lines": generated lines (list, spLines)')
  return(results)
}



#' Elevation profile
#'
#' Extracts the minimum and maximum elevation data along a swathR profile.
#'
#' @param x list. an return object from \code{swathR}
#' @return data.frame
#' @importFrom dplyr c_across rowwise ungroup mutate
#' @export
swath_profile <- function(x) {
  elevs <- c()
  center <- as.character(median(seq_along(x$data)))
  for (i in seq_along(x$data)) {
    elevs <- cbind(elevs, x$data[[i]])
  }
  elevs.df <- as.data.frame(elevs)
  names(elevs.df) <- seq_along(x$data)

  elevs.df <- elevs.df |>
    dplyr::rowwise()  |>
    dplyr::mutate(
      min = min(dplyr::c_across()),
      max = max(dplyr::c_across())
    ) %>%
    dplyr::ungroup()
  elevs.df <-
    data.frame(
      seq_along(elevs.df$min),
      elevs.df[center],
      elevs.df$min,
      elevs.df$max
    )
  names(elevs.df) <- c("distance", "elevation", "min", "max")
  return(elevs.df)
}

#' @title Distance
#' @description This uses the **haversine** formula to calculate the great-circle
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
greatcircle_distance <- function(a, b) {
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

#' Distances in degree to kilometer
#'
#' Converts distances along a great circle path from degree into kilometer
#'
#' @param x numeric vector of distances in degree
#' @param start,end start and end point as vectors with lon, lat
#' @return numeric vector
#' @importFrom scales rescale
#' @export
deg_2_km <- function(x, start, end) {
  distance.km <- greatcircle_distance(start, end)
  scales::rescale(x, to = c(0, distance.km))
}

