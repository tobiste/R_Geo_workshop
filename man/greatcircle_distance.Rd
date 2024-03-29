% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swath.R
\name{greatcircle_distance}
\alias{greatcircle_distance}
\title{Distance}
\usage{
greatcircle_distance(a, b)
}
\arguments{
\item{a}{lon, lat coordinate of point 1}

\item{b}{lon, lat coordinate of point 2}
}
\value{
distance in km
}
\description{
This uses the \strong{haversine} formula to calculate the great-circle
distance between two points – that is, the shortest distance over the earth’s
surface – giving an ‘as-the-crow-flies’ distance between the points
(ignoring any hills they fly over, of course!).
}
\examples{
berlin <- c(52.517, 13.4)
tokyo <- c(35.7, 139.767)
haversine(berlin, tokyo)
}
