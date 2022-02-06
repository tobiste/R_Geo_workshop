clean_header <- function(x){
  x <- stringr::str_replace(x, '[ ]+', "_")
  x <- stringr::str_replace(x, pattern = '[(]', replacement = "")
  x <- stringr::str_replace(x, pattern = '[)]', replacement = "")
  x <- stringr::str_replace(x, ' ', "_")
  x <- stringr::str_replace(x, ' ', "_")
  x <- stringr::str_squish(x)
  return(x)
}


#' @title Read data from Geochron database
#'
#' @description Converts Excel spreadsheets from a Geochon download into a meta
#' data and isotope data table. \code{meta} contains all the sample information.
#' \code{isotopes} contains all the information of the measurements (i.e. isotopes, ages)
#'
#' @param path Path to the xls file.
#' @param method Character string defining the input type of the dataset. One of
#' "Ar-Ar", "Fission Track", "(U-Th)/He", and "U-Pb"
#' @return list
#' @export
#' @importFrom dplyr "%>%" filter select
#' @importFrom stringr str_replace
#' @importFrom readxl read_excel
read_geochron <- function(path, method = 'U-Pb') {
  require(dplyr)
  require(stringr)
  require(readxl)

  if (method == "Ar-Ar") {
  }
  if (method == "(U-Th)/He") {
  }
  if (method == "Fission Track") {
    details0 <- readxl::read_excel(path = path, sheet = 'Details', col_names = FALSE)

    details <- data.frame(t(details0$`...2`))
    names(details) <-  clean_header(details0$`...1`)

    meta0 <- readxl::read_excel(path = path, sheet = 'Metadata', col_names = FALSE)
    meta <- data.frame(t(meta0$`...2`))
    names(meta) <- clean_header(meta0$`...1`)

    meta <- cbind(meta, details)
    ages <- readxl::read_excel(path = path, sheet = 'Ages')

    return(list(meta=meta, tracks=ages))

  }
  if (method == "U-Pb") {
    # meta
    meta0 <- readxl::read_excel(path = path,
                                sheet = 1,
                                skip = 7)
    meta.header <- colnames(meta0)
    colnames(meta0) <- clean_header(colnames(meta0))

    meta <- meta0 %>%
      filter(!is.na(Country)) %>%
      filter(Unique_ID != 'Unique ID') %>%
      mutate(
        Longitude = as.numeric(Longitude),
        Latitude = as.numeric(Latitude),
        Max_Age_Ma = as.numeric(Max_Age_Ma),
        Min_Age_Ma = as.numeric(Min_Age_Ma),
        Oldest_Frac._Date_Ma = as.numeric(Oldest_Frac._Date_Ma),
        Youngest_Frac._Date_Ma = as.numeric(Youngest_Frac._Date_Ma)
      )

    # isotopes
    isotopes.header <- t(meta0[2,])
    rownames(isotopes.header) <- NULL
    isotopes.header[1] <- "Sample ID"
    isotopes.header[2] <- "Unique ID"
    #isotopes.header <- na.omit(isotopes.header)
    isotopes.header <-
      stringr::str_replace(isotopes.header, '[ ]', "_")
    isotopes.header <-
      stringr::str_replace(isotopes.header, "/", ".")
    isotopes.header <-
      stringr::str_replace(isotopes.header, '[ ]', "_")

    isotopes0 <-
      readxl::read_excel(path = path,
                         sheet = 1,
                         skip = 7)
    colnames(isotopes0) <- isotopes.header


    isotopes <- isotopes0[FALSE, ]
    for (i in seq_along(isotopes0$Sample_ID)) {
      if (!is.na(isotopes0$Sample_ID[i]) &
          !is.na(isotopes0$Unique_ID[i])) {
        sample_id <- isotopes0$Sample_ID[i]
        unique_id <- isotopes0$Unique_ID[i]
      }
      if (is.na(isotopes0$Sample_ID[i]) &
          is.na(isotopes0$Unique_ID[i])) {
        isotopes.i <- isotopes0[i,]
        isotopes.i$Sample_ID <- sample_id
        isotopes.i$Unique_ID <- unique_id

        isotopes <- rbind(isotopes, isotopes.i)
      }
    }
    isotopes <- isotopes %>%
      dplyr::select(na.omit(isotopes.header)) %>%
      mutate(
        `206.238_Age` = as.numeric(`206.238_Age`),
        `206.238_Age_Error` = as.numeric(`206.238_Age_Error`),
        `207.235_Age` = as.numeric(`207.235_Age`),
        `207.235_Age_Error` = as.numeric(`207.235_Age_Error`),
        `207.206_Age` = as.numeric(`207.206_Age`),
        `207.206_Age_Error` = as.numeric(`207.206_Age_Error`),
        `Pb*.Pbc` = as.numeric(`Pb*.Pbc`),
        Rho  = as.numeric(Rho),
        Rho_Error = as.numeric(Rho_Error),
        `206.238` = as.numeric(`206.238`),
        `206.238_Error` = as.numeric(`206.238_Error`),
        `206.204` = as.numeric(`206.204`),
        `208.206` = as.numeric(`208.206`),
        conc_U = as.numeric(conc_U),
        `Th.U_samp` = as.numeric(`Th.U_samp`),
        `207.235_Age` = as.numeric(`207.235_Age`),
        `207.235_Age_Error` = as.numeric(`207.235_Age_Error`),
        `207.206_Age` = as.numeric(`207.206_Age`),
        `207.206_Age_Error` = as.numeric(`207.206_Age_Error`),
        `206.238xTh_Age` = as.numeric(`206.238xTh_Age`),
        `206.238xTh_Age_Error` = as.numeric(`206.238xTh_Age_Error`),
        `207.235xPa_Age` = as.numeric(`207.235xPa_Age`),
        `207.235xPa_Age_Error` = as.numeric(`207.235xPa_Age_Error`),
        `207.206xTh_Age` = as.numeric(`207.206xTh_Age`),
        `207.206xTh_Age_Error` = as.numeric(`207.206xTh_Age_Error`),
        `207.206xPa_Age` = as.numeric(`207.206xPa_Age`),
        `207.206xPa_Age_Error` = as.numeric(`207.206xPa_Age_Error`)
      )
    colnames(isotopes) <- isotopes.header <- c(
      "Sample_ID", "Unique_ID", "Fraction_ID",
      "t.Pb206U238", "st.Pb206U238",
      "t.Pb207U235", "st.Pb207U235",
      "t.Pb207Pb206",  "st.Pb207Pb206",
      "PbPb.cor", "rho", "s.rho",
      "Pb206U238", "errPb206U238",
      "Pb206Pb204", "Pb208Pb206", "U", "ThU",
      "Age_206.238xTh", "Age_Error_206.238xTh",
      "Age_207.235xPa", "Age_Error_207.235xPar",
      "Age_207.206xTh", "Age_Error_207.206xTh",
      "Age_207.206xPa_Age", "Age_Error_207.206xPa"
    )


    return(list(meta = meta, isotopes = isotopes))
  }
}

# color palette ####
colorblind <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", '#7555C8', '#7EF791', '#85112A', '#2BE6BE', '#EE599F', '#5F9143', '#759301') # Bang Wong


geochron_to_dz <- function(x){
  lst <- list()
  for(i in x$Sample_ID){
    x.i <- subset(x, x$Sample_ID==i)
    ages <- c(x.i$t)
    #names(ages) <- i
    lst[[i]] = ages
  }

  class(lst) <- "detritals"
  return(lst)
}




#' @title swathR v1.0.1
#' @author V. Haburaj
#' @description Calculate swath-profile values perpendicular to a straight baseline. The
#' baseline is generated between two user-defined points (X|Y), see argument
#' \code{'coords'}. The distance between samples and the number of samples can be
#' specified, see arguments \code{'k'} and \code{'dist'}. Values of the swath-profile are
#' extracted from a given raster file, see argument \code{'raster'}. CRS of raster
#' and points have to be the same.
#'@param coords matrix(ncol=2, nrow=2) with x and y coordinates of beginning and
#' end point of the baseline; each point in one row
#' {\describe
#' \item{column 1}{xcoordinates}
#' \item{column 2}{ycoordinates}
#' }
#'@param raster raster file (loaded with \code{raster()} from package 'raster')
#'@param k integer; number of lines on each side of the baseline
#'@param dist numeric; distance between lines
#'@param crs string; CRS
#'@param method = string; method for extraction of raw data, see
#'\code{extract()} from package 'raster': default value: 'bilinear'
#' @import sp
#' @import rgeos
#' @import raster
#' @export
swathR <- function(coords, raster, k, dist, crs, method) {
  require(sp)
  require(rgeos)
  require(raster)
  message('Initializing ...')
  # set default method:
  if (missing(method)) {
    method <- 'bilinear'
  }
  # create SpatialPoints from coords:
  spt <- SpatialPoints(coords, proj4string = CRS(crs))
  # get slope of baseline:
  m <- (ymin(spt[1]) - ymin(spt[2])) / (xmin(spt[1]) - xmin(spt[2]))
  # get slope of normal function:
  m1 <- - (1/m)
  # get slope-angle from slope:
  alpha <- atan(m)
  alpha1 <- atan(m1)
  # get deltax and deltay from pythagoras:
  if ((alpha*180)/pi < 90 & (alpha*180)/pi > 270) {
    deltax <- cos(alpha1) * dist
  } else deltax <- cos(alpha1) * dist * -1
  if ((alpha*180)/pi > 0  & (alpha*180)/pi < 180) {
    deltay <- sqrt(dist**2 - deltax**2)
  } else deltay <- sqrt(dist**2 - deltax**2) * -1
  # create empty matrix:
  swath <- matrix(nrow = k*2+1, ncol = 8)
  colnames(swath) <- c('distance', 'mean', 'median',
                       'std.dev.', 'min', 'max', 'quantile(25)', 'quantile(75)')
  # list for spatial lines:
  allLines <- list()
  # add baseline:
  allLines[[k+1]] <- spLines(coords, crs=crs)
  # set distance for baseline:
  swath[k+1, 1] <- 0
  # generate k lines parallel to baseline:
  for (n in 1:k) {
    # BELOW BASELINE:
    # new points
    cn <- matrix(nrow=2, ncol=2)
    cn[1,] <- cbind(coords[1,1] - (deltax * n), coords[1,2] - (deltay * n))
    cn[2,] <- cbind(coords[2,1] - (deltax * n), coords[2,2] - (deltay * n))
    # line between points:
    allLines[[k+1-n]] <- spLines(cn, crs=crs)
    # distance value:
    swath[k+1-n,1] <- -1*n*dist
    # ABOVE BASELINE:
    # new points
    cn <- matrix(nrow=2, ncol=2)
    cn[1,] <- cbind(coords[1,1] + (deltax * n), coords[1,2] + (deltay * n))
    cn[2,] <- cbind(coords[2,1] + (deltax * n), coords[2,2] + (deltay * n))
    # line between points:
    allLines[[k+n+1]] <- spLines(cn, crs=crs)
    # distance value:
    swath[k+n+1, 1] <- n*dist
  }
  gc(verbose = FALSE)
  # get raw data:
  message('Extracting raw data (this may take some time) ...')
  raw.data <- sapply(allLines, FUN=function(x) {extract(raster, x, method=method)})
  gc(verbose = FALSE)
  # generalise data:
  message('Generalising data ...')
  swath[,2] <- sapply(raw.data, function(x) {mean(x, na.rm = T)})
  swath[,3] <- sapply(raw.data, function(x) {median(x, na.rm = T)})
  swath[,4] <- sapply(raw.data, function(x) {sd(x, na.rm = T)})
  swath[,5] <- sapply(raw.data, function(x) {min(x, na.rm = T)})
  swath[,6] <- sapply(raw.data, function(x) {max(x, na.rm = T)})
  swath[,7] <- sapply(raw.data, function(x) {quantile(x, na.rm = T)[2]})
  swath[,8] <- sapply(raw.data, function(x) {quantile(x, na.rm = T)[4]})
  # return results:
  results <- list(swath=swath, data=raw.data, lines=allLines)
  message('Operation finished successfully!')
  message('Structure of results (list): "swath": swath profile data (matrix, numeric), "data": raw data (list, numeric), "lines": generated lines (list, spLines)')
  gc(verbose = FALSE)
  return(results)
}
