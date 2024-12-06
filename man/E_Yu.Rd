% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LST_terra.R
\name{E_Yu}
\alias{E_Yu}
\title{Land Surface Emissivity according to Yu et al. 2014}
\usage{
E_Yu(red = red, NDVI = NDVI, band = band)
}
\arguments{
\item{red}{SpatRaster object, red band of remote sensing imagery}

\item{NDVI}{SpatRaster object, NDVI calculated from remote sensing imagery}

\item{band}{A string specifying which Landsat 8 thermal band to use. It can be "band 10" or
"band 11"}
}
\value{
SpatRaster
}
\description{
This function calculates Land Surface Emissivity according to Yu et al. 2014
}
\examples{
red <- terra::rast(ncol=100, nrow=100)
set.seed(2)
terra::values(red) = runif(10000, min=0.1, max=0.4)
NDVI <- terra::rast(ncol=100, nrow=100)
set.seed(2)
terra::values(NDVI) = runif(10000, min=0.02, max=0.8)
E_Yu(red = red, NDVI = NDVI, band = "band 11")
}
\references{
Yu, X., Guo, X. and Wu, Z., 2014. Land surface temperature retrieval from Landsat 8 TIRS—Comparison between radiative transfer equation-based method, split window algorithm and single channel method. Remote sensing, 6(10), pp.9829-9852.
}
