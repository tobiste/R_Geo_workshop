#' @title swathR v1.0.1
#' @author V. Haburaj
#' @description Calculate swath-profile values perpendicular to a straight baseline. The
#' baseline is generated between two user-defined points (X|Y), see argument
#' \code{'coords'}. The distance between samples and the number of samples can be
#' specified, see arguments \code{'k'} and \code{'dist'}. Values of the swath-profile are
#' extracted from a given raster file, see argument \code{'raster'}. CRS of raster
#' and points have to be the same.
#' @param coords matrix(ncol=2, nrow=2) with x and y coordinates of beginning and
#' end point of the baseline; each point in one row
#' \describe{
#'   \item{column 1}{xcoordinates}
#'   \item{column 2}{ycoordinates}
#' }
#' @param raster raster file (loaded with \code{raster()} from package 'raster')
#' @param k integer; number of lines on each side of the baseline
#' @param dist numeric; distance between lines
#' @param crs string; CRS
#' @param method string; method for extraction of raw data, see
#' \code{extract()} from package 'raster': default value: 'bilinear'
#' @import sp
#' @import rgeos
#' @import raster
#' @export
swathR <- function(coords, raster, k, dist, crs, method) {
  require(sp)
  require(rgeos)
  require(raster)
  message("Initializing ...")
  # set default method:
  if (missing(method)) {
    method <- "bilinear"
  }
  # create SpatialPoints from coords:
  spt <- SpatialPoints(coords, proj4string = CRS(crs))
  # get slope of baseline:
  m <- (ymin(spt[1]) - ymin(spt[2])) / (xmin(spt[1]) - xmin(spt[2]))
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
  allLines[[k + 1]] <- spLines(coords, crs = crs)
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
    allLines[[k + 1 - n]] <- spLines(cn, crs = crs)
    # distance value:
    swath[k + 1 - n, 1] <- -1 * n * dist
    # ABOVE BASELINE:
    # new points
    cn <- matrix(nrow = 2, ncol = 2)
    cn[1, ] <- cbind(coords[1, 1] + (deltax * n), coords[1, 2] + (deltay * n))
    cn[2, ] <- cbind(coords[2, 1] + (deltax * n), coords[2, 2] + (deltay * n))
    # line between points:
    allLines[[k + n + 1]] <- spLines(cn, crs = crs)
    # distance value:
    swath[k + n + 1, 1] <- n * dist
  }
  gc(verbose = FALSE)
  # get raw data:
  message("Extracting raw data (this may take some time) ...")
  raw.data <- sapply(allLines, FUN = function(x) {
    extract(raster, x, method = method)
  })
  gc(verbose = FALSE)
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
  gc(verbose = FALSE)
  return(results)
}



#' Elevation profile
#'
#' Extracts the minimum and maximum elevation data along a swathR profile.
#'
#' @param x list. an return object from \code{swathR}
#' @return data.frame
#' @import dplyr
#' @export
swath_profile <- function(x) {
  elevs <- c()
  center <- as.character(median(seq_along(x$data)))
  for (i in seq_along(x$data)) {
    elevs <- cbind(elevs, x$data[[i]])
  }
  elevs.df <- as.data.frame(elevs)
  names(elevs.df) <- seq_along(x$data)

  elevs.df <- elevs.df %>%
    rowwise() %>%
    mutate(
      min = min(c_across()),
      max = max(c_across())
    ) %>%
    ungroup()
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
