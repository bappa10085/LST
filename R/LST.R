#' @title At-Sensor Temperature or brightness temperature
#'
#' @description This function calculates at-Sensor Temperature or brightness temperature
#' @importFrom terra rast
#' @param Landsat_10 SpatRaster object, Landsat band 10
#' @param Landsat_11 SpatRaster object, Landsat band 11
#' @return A list containing brightness temperature corresponding to Landsat band 10 and Landsat band 11
#' @export
#' @examples
#' a <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(a) = runif(10000, min=27791, max=30878)
#'
#' b <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(b) = runif(10000, min=25686, max=28069)
#'
#' BT(Landsat_10 = a, Landsat_11 = b)

BT <- function(Landsat_10, Landsat_11){

  #Top of Atmospheric Spectral Radiance
  l_lambda_10 <- 3.3420*10^-4*Landsat_10 + 0.1
  l_lambda_11 <- 3.3420*10^-4*Landsat_11 + 0.1

  #l_lambda = ML * Qcal + AL
  #l_lambda = TOA Spectral Reflectance (watts/m2*Srad µm),
  #ML= The band-specific multiplicative rescaling factor,
  #Qcal is the Band 10 image,
  #AL is the band-specific additive rescaling factor
  #The value for ML and AL can be found in MTL.txt file

  #Conversion of Radiance to At-Sensor Temperature
  #K1_CONSTANT_BAND_10 = 774.8853
  #K2_CONSTANT_BAND_10 = 1321.0789
  BT_10 <- (1321.0789/(log(774.8853/l_lambda_10 + 1)))

  #K1_CONSTANT_BAND_11 <- 480.8883
  #K2_CONSTANT_BAND_11 <- 1201.1442
  BT_11 <- (1201.1442/(log(480.8883/l_lambda_11 + 1)))
  rast(list(BT_10=BT_10, BT_11=BT_11))
}


#' @title NDVI
#'
#' @description Function for NDVI calculation
#' @importFrom terra rast
#' @param Red SpatRaster object, red band of remote sensing imagery
#' @param NIR SpatRaster object, NIR band of remote sensing imagery
#' @return SpatRaster
#' @export
#' @examples
#' red <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(red) = runif(10000, min=0.1, max=0.4)
#'
#' NIR <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(NIR) = runif(10000, min=0.1, max=0.6)
#'
#' NDVI(Red = red, NIR = NIR)

NDVI <- function(Red, NIR){

  #NDVI calculation
  ndvi <- (NIR - Red)/(NIR + Red)
  return(ndvi)
}

#' @title Proportion of vegetation or fractional vegetation cover
#'
#' @description Calculation of the proportion of vegetation or fractional vegetation cover from NDVI
#' @param NDVI SpatRaster object, NDVI calculated from remote sensing imagery
#' @param minNDVI = 0.2 (Ref. Sobrino et al. 2004)
#' @param maxNDVI = 0.5 (Ref. Sobrino et al. 2004)
#' @return SpatRaster
#' @references Sobrino, J.A., Jiménez-Muñoz, J.C. and Paolini, L., 2004. Land surface temperature retrieval from LANDSAT TM 5. Remote Sensing of environment, 90(4), pp.434-440.
#' @export
#' @examples
#' NDVI <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' Pv(NDVI = NDVI, minNDVI = 0.2, maxNDVI = 0.5)

Pv <- function(NDVI, minNDVI, maxNDVI){

  Pv <- ((NDVI - minNDVI)/(maxNDVI - minNDVI))^2
  return(Pv)
}

#' @title Land Surface Emissivity according to Van de Griend and Owe 1993
#'
#' @description This function calculates Land Surface Emissivity according to Van de Griend and Owe 1993
#' @param NDVI SpatRaster object, NDVI calculated from remote sensing imagery
#' @return SpatRaster
#' @references Van de Griend, A.A. and Owe, M., 1993. On the relationship between thermal emissivity and the normalized difference vegetation index for natural surfaces. International Journal of remote sensing, 14(6), pp.1119-1131.
#' @export
#' @examples
#' NDVI <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_VandeGriend(NDVI)

E_VandeGriend <- function(NDVI){

  E_VandeGriend <- 1.094 + 0.047*log(NDVI)
  #E_VandeGriend[E_VandeGriend>1] <- 1
  return(E_VandeGriend)
}

#' @title Land Surface Emissivity according to Valor and Caselles 1996
#'
#' @description This function calculates Land Surface Emissivity according to Valor and Caselles 1996
#' @param NDVI SpatRaster object, NDVI calculated from remote sensing imagery
#' @return SpatRaster
#' @references Valor, E. and Caselles, V., 1996. Mapping land surface emissivity from NDVI: Application to European, African, and South American areas. Remote sensing of Environment, 57(3), pp.167-184.
#' @export
#' @examples
#' NDVI <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_Valor(NDVI)

E_Valor <- function(NDVI){

  #Pv calculation
  Pv <- ((NDVI - 0.2)/(0.5 - 0.2))^2
  #Emissivity calculation
  E_Valor <- 0.985*Pv + 0.960*(1 - Pv) + 0.06*Pv*(1 - Pv)
  return(E_Valor)
}

#' @title Atmospheric transmittance calculation
#'
#' @description This function calculates Atmospheric transmittance from near-surface air temperature (To, °C) and relative humidity (RH, %) of the date when Landsat passed over the study area
#' @param To Near-surface air temperature (°C) of the date when Landsat passed over the study area
#' @param RH relative humidity (%) of the date when Landsat passed over the study area
#' @param band A string specifying which Landsat 8 thermal band to use. It can be "band 10" or
#' "band 11"
#' @return Atmospheric transmittance
#' @export
#' @examples
#' tau(To = 26, RH = 42, band = "band 11")

tau <- function(To = To, RH = To, band = band){
  To_K <- To + 273.15
  RH_fraction <- RH/100
  W <- 0.0981*10*0.6108*exp((17.27*(To_K-273.15))/
                              (237.3+(To_K-273.15)))*RH_fraction+0.1697
  if (band == "band 10") {
    tau <- -0.0164*W^2-0.04203*W+0.9715
  }else {
    tau <- -0.01218*W^2-0.07735*W+0.9603
  }
  return(tau)
}

#' @title Land Surface Emissivity according to Sobrino et al. 2008
#'
#' @description This function calculates Land Surface Emissivity according to Sobrino et al. 2008
#' @importFrom terra ifel
#' @param red SpatRaster object, red band of remote sensing imagery
#' @param NDVI SpatRaster object, NDVI calculated from remote sensing imagery
#' @return SpatRaster
#' @references Sobrino, J.A., Jiménez-Muñoz, J.C., Sòria, G., Romaguera, M., Guanter, L., Moreno, J., Plaza, A. and Martínez, P., 2008. Land surface emissivity retrieval from different VNIR and TIR sensors. IEEE transactions on geoscience and remote sensing, 46(2), pp.316-327.
#' @export
#' @examples
#' red <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(red) = runif(10000, min=0.1, max=0.4)
#' NDVI <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_Sobrino(red = red, NDVI = NDVI)

E_Sobrino <- function(red = red, NDVI = NDVI){

  #Pv calculation
  Pv <- ((NDVI - 0.2)/(0.5 - 0.2))^2
  E1_Sobrino <- 0.004*Pv + 0.986
  red <- (0.979 - 0.035*red)

  #where NDVI < 0.2, take the value from red band otherwise use NA:
  E_Sobrino <- terra::ifel(NDVI < 0.2, red, NA)

  #where NDVI is over 0.5 set E to 0.99, otherwise use whatever was in E:
  E_Sobrino = terra::ifel(NDVI > 0.5, 0.99, E_Sobrino)

  #if NDVI is between 0.2 and 0.5 take the value from b otherwise use whatever was in E:
  E_Sobrino = terra::ifel(NDVI>=0.2 & NDVI<=0.5, E1_Sobrino, E_Sobrino)

  return(E_Sobrino)
}

#' @title Mono window algorithm
#'
#' @description This function calculates Land Surface Temperature using mono window algorithm
#' @param BT SpatRaster object, brightness temperature
#' @param tau Atmospheric transmittance
#' @param E SpatRaster object, Land Surface Emissivity calculated according to Van de Griend and Owe 1993 or Valor and Caselles 1996 or Sobrino et al. 2008
#' @param Ta Mean atmospheric temperature (K) of the date when Landsat passed over the study area
#' @return SpatRaster
#' @references Qin, Z., Karnieli, A. and Berliner, P., 2001. A mono-window algorithm for retrieving land surface temperature from Landsat TM data and its application to the Israel-Egypt border region. International journal of remote sensing, 22(18), pp.3719-3746.
#' @export
#' @examples
#' BTemp <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(BTemp) = runif(10000, min=298, max=305)
#' E <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(E) = runif(10000, min=0.96, max=0.99)
#' MWA(BT = BTemp, tau = 0.86, E = E, Ta = 26)

MWA <- function(BT = BT, tau = tau, E = E, Ta = Ta){
  C <- E*tau
  D <- (1-tau)*(1+(1-E)*tau)
  T_MWA <- (-67.355351*(1-C-D)+(0.458606*(1-C-D)+C+D)*BT-D*Ta)/C
  return(T_MWA)
}

#' @title Mean atmospheric temperature
#'
#' @description This function calculates mean atmospheric temperature (Ta) using near-surface air
#' temperature (To)
#' @param To Near-surface air temperature (°C) of the date when Landsat passed over the study area
#' @param mod A string specifying which model to use. It can be anyone of "USA 1976 Standard" or
#' "Tropical Region" or "Mid-latitude Summer Region" or "Mid-latitude Winter Region"
#' @return Mean atmospheric temperature (K)
#' @references Sekertekin, A. and Bonafoni, S., 2020. Land surface temperature retrieval from Landsat 5, 7, and 8 over rural areas: Assessment of different retrieval algorithms and emissivity models and toolbox implementation. Remote sensing, 12(2), p.294.
#' @export
#' @examples
#' Ta(To = 26, mod = "Mid-latitude Winter Region")

Ta <- function(To = To, mod = mod){
  To_K <- To + 273.15
  if (mod == "USA 1976 Standard") {
    Ta <- 25.940 + 0.8805*To_K
  } else if (mod == "Tropical Region") {
    Ta <- 17.977 + 0.9172*To_K
  } else if (mod == "Mid-latitude Summer Region") {
    Ta <- 16.011 + 0.9262*To_K
  } else {
    Ta <- 19.270 + 0.9112*To_K
  }
  return(Ta)
}

#' @title Single channel algorithm
#'
#' @description This function calculates Land Surface Temperature using single channel algorithm
#' @param TIR SpatRaster object, Landsat band 10 or 11
#' @param tau Atmospheric transmittance
#' @param E SpatRaster object, Land Surface Emissivity calculated according to Van de Griend and Owe 1993 or Valor and Caselles 1996 or Sobrino et al. 2008
#' @param dlrad Downwelling radiance calculated from https://atmcorr.gsfc.nasa.gov/
#' @param ulrad upwelling radiance calculated from https://atmcorr.gsfc.nasa.gov/
#' @param band A string specifying which Landsat 8 thermal band to use. It can be "band 10" or
#' "band 11"
#' @return SpatRaster
#' @references Jimenez-Munoz, J.C., Cristobal, J., Sobrino, J.A., Sòria, G., Ninyerola, M. and Pons, X., 2008. Revision of the single-channel algorithm for land surface temperature retrieval from Landsat thermal-infrared data. IEEE Transactions on geoscience and remote sensing, 47(1), pp.339-349.
#' @export
#' @examples
#' TIR <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(TIR) = runif(10000, min=27791, max=30878)
#' E <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(E) = runif(10000, min=0.96, max=0.99)
#' Ts_SCA <- SCA(TIR = TIR, tau = 0.86, E = E,
#' 		dlrad = 2.17, ulrad = 1.30, band = "band 11")

SCA <- function(TIR = TIR, tau = tau, E = E, dlrad = dlrad, ulrad = ulrad, band = band){

  if (band == "band 10") {
    #Top of Atmospheric Spectral Radiance
    l_lambda <- 3.3420*10^-4*TIR + 0.1
    BT <- (1321.0789/(log(774.8853/l_lambda + 1)))
  } else{
    l_lambda <- 3.3420*10^-4*TIR + 0.1
    BT <- (1201.1442/(log(480.8883/l_lambda + 1)))
  }

  c2 <- (6.62607015*10^-34)*(2.99792458*10^8)*1000000/(1.380649*10^-23)

  #C2 = h*c/s in meter.Kelvin
  #Where, h = Planck's constant, 6.62607015*10^-34 JS,
  #S = Boltzmann constant, 1.380649*10^-23 J/K
  #C = Speed of light, 2.99792458* 10^8 m/s

  if (band == "band 10") {
    by <- c2/mean(c(10.6,11.19))
  } else{
    by <- c2/mean(c(11.5, 12.51))
  }

  #by = c2/Effective band wavelength
  #wavelength of band 10 of Landsat 8: 10.6-11.19 micrometers
  #wavelength of band 11 of Landsat 8: 11.5-12.51 micrometers

  gamma <- BT^2/(by * l_lambda)
  delta <- BT - BT^2/by

  psi1 <- 1/tau
  psi2 <- -dlrad-(ulrad/tau)
  psi3 <- dlrad

  Ts_SCA <- gamma*((1/E)*(psi1*l_lambda+psi2)+psi3)+delta
  return(Ts_SCA)
}

#' @title Radiative transfer equation method
#'
#' @description This function calculates Land Surface Temperature using radiative transfer equation method
#' @param TIR SpatRaster object, Landsat band 10 or 11
#' @param tau Atmospheric transmittance
#' @param E SpatRaster object, Land Surface Emissivity calculated according to Van de Griend and Owe 1993 or Valor and Caselles 1996 or Sobrino et al. 2008
#' @param dlrad Downwelling radiance calculated from https://atmcorr.gsfc.nasa.gov/
#' @param ulrad upwelling radiance calculated from https://atmcorr.gsfc.nasa.gov/
#' @param band A string specifying which Landsat 8 thermal band to use. It can be "band 10" or
#' "band 11"
#' @return SpatRaster
#' @references Srivastava, P.K., Majumdar, T.J. and Bhattacharya, A.K., 2009. Surface temperature estimation in Singhbhum Shear Zone of India using Landsat-7 ETM+ thermal infrared data. Advances in space research, 43(10), pp.1563-1574.
#' @export
#' @examples
#' TIR <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(TIR) = runif(10000, min=27791, max=30878)
#' BT <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(BT) = runif(10000, min=298, max=305)
#' E <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(E) = runif(10000, min=0.96, max=0.99)
#' Ts_RTE <- RTE(TIR = TIR, tau = 0.86, E = E,
#' 		dlrad = 2.17, ulrad = 1.30, band = "band 11")

RTE <- function(TIR = TIR, tau = tau, E = E, dlrad = dlrad, ulrad = ulrad, band = band){

  if (band == "band 10") {
    #Top of Atmospheric Spectral Radiance
    l_lambda <- 3.3420*10^-4*TIR + 0.1
  } else{
    l_lambda <- 3.3420*10^-4*TIR + 0.1
  }
  l_lambda_Ts <- ((l_lambda - ulrad)/(tau*E))- ((1-E)*dlrad/E)
  if (band == "band 10") {
    Ts_RTE <- 1321.08/log(774.89/l_lambda_Ts+1)
    #Ts_RTE <- K2/log(K1/l_lambda_Ts+1)
  } else{
    Ts_RTE <- 1201.14/log(480.89/l_lambda_Ts+1)
  }
  return(Ts_RTE)
}

#' @title Land Surface Emissivity according to Skokovic et al. 2014
#'
#' @description This function calculates Land Surface Emissivity according to Skokovic et al. 2014
#' @param red SpatRaster object, red band of remote sensing imagery
#' @param NDVI SpatRaster object, NDVI calculated from remote sensing imagery
#' @param band A string specifying which Landsat 8 thermal band to use. It can be "band 10" or
#' "band 11"
#' @return SpatRaster
#' @references Skoković, D., Sobrino, J.A., Jimenez-Munoz, J.C., Soria, G., Julien, Y., Mattar, C. and Cristóbal, J., 2014. Calibration and Validation of land surface temperature for Landsat8-TIRS sensor. Land product validation and evolution.
#' @export
#' @examples
#' red <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(red) = runif(10000, min=0.1, max=0.4)
#' NDVI <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_Skokovic(red = red, NDVI = NDVI, band = "band 11")

E_Skokovic <- function(red = red, NDVI = NDVI, band = band){

  #Pv calculation
  Pv <- ((NDVI - 0.2)/(0.5 - 0.2))^2
  if (band == "band 10") {
    dE_10 <- (1 - 0.9668)*(1 - Pv)*0.55*0.9863
    #Es = Emissivity of soil (0.9668) for Landsat 8 TIRS band 10
    #Ev = Emissivity of vegetation (0.9863) for Landsat 8 TIRS band 10 (Yu et al., 2014)
    E1_Skokovic_10 <- 0.987*Pv + 0.971*(1 - Pv) + dE_10

    red_10 <- (0.979 - 0.046*red)
    #where NDVI < 0.2, take the value from red band otherwise use NA:
    E_Skokovic_10 <- terra::ifel(NDVI < 0.2, red_10, NA)

    #where NDVI is over 0.5 set E to 0.987, otherwise use whatever was in E:
    E_Skokovic_10 = terra::ifel(NDVI > 0.5, 0.987+dE_10, E_Skokovic_10)

    #if NDVI is between 0.2 and 0.5 take the value from b otherwise use whatever was in E:
    E_Skokovic_10 = terra::ifel(NDVI>=0.2 & NDVI<=0.5, E1_Skokovic_10, E_Skokovic_10)

    return(E_Skokovic_10)
  }else {
    #Band 11
    dE_11 <- (1 - 0.9747)*(1 - Pv)*0.55*0.9896
    E1_Skokovic_11 <- 0.989*Pv + 0.977*(1 - Pv) + dE_11
    #TM6 soil and vegetation emissivities of 0.97 and 0.99, respectively

    red_11 <- (0.982 - 0.027*red)
    #where NDVI < 0.2, take the value from red band otherwise use NA:
    E_Skokovic_11 = terra::ifel(NDVI < 0.2, red_11, NA)

    #where NDVI is over 0.5 set E to 0.987, otherwise use whatever was in E:
    E_Skokovic_11 = terra::ifel(NDVI>0.5, 0.989+dE_11, E_Skokovic_11)

    #if NDVI is between 0.2 and 0.5 take the value from b otherwise use whatever was in E:
    E_Skokovic_11 = terra::ifel(NDVI>=0.2 & NDVI<=0.5, E1_Skokovic_11, E_Skokovic_11)
    return(E_Skokovic_11)
  }
}

#' @title Land Surface Emissivity according to Yu et al. 2014
#'
#' @description This function calculates Land Surface Emissivity according to Yu et al. 2014
#' @param red SpatRaster object, red band of remote sensing imagery
#' @param NDVI SpatRaster object, NDVI calculated from remote sensing imagery
#' @param band A string specifying which Landsat 8 thermal band to use. It can be "band 10" or
#' "band 11"
#' @return SpatRaster
#' @references Yu, X., Guo, X. and Wu, Z., 2014. Land surface temperature retrieval from Landsat 8 TIRS—Comparison between radiative transfer equation-based method, split window algorithm and single channel method. Remote sensing, 6(10), pp.9829-9852.
#' @export
#' @examples
#' red <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(red) = runif(10000, min=0.1, max=0.4)
#' NDVI <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(NDVI) = runif(10000, min=0.02, max=0.8)
#' E_Yu(red = red, NDVI = NDVI, band = "band 11")

E_Yu <- function(red = red, NDVI = NDVI, band = band){

  #Pv calculation
  Pv <- ((NDVI - 0.2)/(0.5 - 0.2))^2
  if (band == "band 10") {
    dE_10 <- (1 - 0.9668)*(1 - Pv)*0.55*0.9863
    #Es = Emissivity of soil (0.9668) for Landsat 8 TIRS band 10
    #Ev = Emissivity of vegetation (0.9863) for Landsat 8 TIRS band 10 (Yu et al., 2014)
    E1_Yu_10 <- 0.9863*Pv + 0.9668*(1 - Pv) + dE_10

    red_10 <- (0.973 - 0.047*red)
    #where NDVI < 0.2, take the value from red band otherwise use NA:
    E_Yu_10 = terra::ifel(NDVI<0.2, red_10[], NA)

    #where NDVI is over 0.5 set E to 0.9863, otherwise use whatever was in E:
    E_Yu_10 = terra::ifel(NDVI>0.5, 0.9863+dE_10, E_Yu_10)

    #if NDVI is between 0.2 and 0.5 take the value from b otherwise use whatever was in E:
    E_Yu_10 = terra::ifel(NDVI>=0.2 & NDVI<=0.5, E1_Yu_10, E_Yu_10)
    return(E_Yu_10)
  }else {
    #Band 11
    dE_11 <- (1 - 0.9747)*(1 - Pv)*0.55*0.9896
    #Es = Emissivity of soil (0.9747) for Landsat 8 TIRS band 11
    #Ev = Emissivity of vegetation (0.9896) for Landsat 8 TIRS band 11 (Yu et al., 2014)
    E1_Yu_11 <- 0.9896*Pv + 0.9747*(1 - Pv) + dE_11

    red_11 <- (0.984 - 0.0026*red)
    #where NDVI < 0.2, take the value from red band otherwise use NA:
    E_Yu_11 = terra::ifel(NDVI<0.2, red_11, NA)

    #where NDVI is over 0.5 set E to 0.9896, otherwise use whatever was in E:
    E_Yu_11 = terra::ifel(NDVI>0.5, 0.9896+dE_11, E_Yu_11)

    #if NDVI is between 0.2 and 0.5 take the value from b otherwise use whatever was in E:
    E_Yu_11 = terra::ifel(NDVI>=0.2 & NDVI<=0.5, E1_Yu_11, E_Yu_11)
    return(E_Yu_11)
  }
}

#' @title Split-window algorithm
#'
#' @description This function calculates Land Surface Temperature using split-window algorithm
#' @param TIR_10 SpatRaster object, Landsat band 10
#' @param TIR_11 SpatRaster object, Landsat band 11
#' @param tau_10 Atmospheric transmittance for Landsat band 10
#' @param tau_11 Atmospheric transmittance for Landsat band 11
#' @param E_10 SpatRaster object, Land Surface Emissivity for Landsat band 10 calculated according to Skokovic et al. 2014 or Yu et al. 2014
#' @param E_11 SpatRaster object, Land Surface Emissivity for Landsat band 11 calculated according to Skokovic et al. 2014 or Yu et al. 2014
#' @return SpatRaster
#' @references Yu, X., Guo, X. and Wu, Z., 2014. Land surface temperature retrieval from Landsat 8 TIRS—Comparison between radiative transfer equation-based method, split window algorithm and single channel method. Remote sensing, 6(10), pp.9829-9852.
#' @export
#' @examples
#' TIR_10 <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(TIR_10) = runif(10000, min=27791, max=30878)
#' TIR_11 <- terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(TIR_11) = runif(10000, min=25686, max=28069)
#' E_10 <- terra::rast(ncol=100, nrow=100)
#' set.seed(1)
#' terra::values(E_10) = runif(10000, min=0.96, max=0.99)
#' E_11 <-terra::rast(ncol=100, nrow=100)
#' set.seed(2)
#' terra::values(E_11) = runif(10000, min=0.96, max=0.99)
#' Ts_SWA <- SWA(TIR_10=TIR_10, TIR_11=TIR_11, tau_10=0.86,
#' 		tau_11=0.87, E_10=E_10, E_11=E_11)

SWA <- function(TIR_10 = TIR_10, TIR_11 = TIR_11, tau_10=tau_10, tau_11=tau_11, E_10 = E_10, E_11 = E_11){

  #Top of Atmospheric Spectral Radiance
  l_lambda_10 <- 3.3420*10^-4*TIR_10 + 0.1
  l_lambda_11 <- 3.3420*10^-4*TIR_11 + 0.1

  #l_lambda = ML * Qcal + AL
  #l_lambda = TOA Spectral Reflectance (watts/m2*Srad µm),
  #ML= The band-specific multiplicative rescaling factor,
  #Qcal is the Band 10 image,
  #AL is the band-specific additive rescaling factor
  #The value for ML and AL can be found in MTL.txt file

  #Conversion of Radiance to At-Sensor Temperature
  #K1_CONSTANT_BAND_10 = 774.8853
  #K2_CONSTANT_BAND_10 = 1321.0789
  BT_10 <- (1321.0789/(log(774.8853/l_lambda_10 + 1)))

  #K1_CONSTANT_BAND_11 <- 480.8883
  #K2_CONSTANT_BAND_11 <- 1201.1442
  BT_11 <- (1201.1442/(log(480.8883/l_lambda_11 + 1)))
  a10 <- E_10*tau_10
  c10 <- (1-tau_10)*(1+(1-E_10)*tau_10)
  L10 <- 0.4464*BT_10 - 66.61
  a11 <- E_11*tau_11
  c11 <- (1-tau_11)*(1+(1-E_11)*tau_11)
  L11 <- 0.4831*BT_11 - 71.23
  b0 <- (c11*(1-a10-c10)*L10-c10*(1-a11-c11)*L11)/(c11*a10-c10*a11)
  b1 <- c10/(c11*a10-c10*a11)
  T_SWA <- BT_10 + b1*(BT_10 - BT_11) + b0
  return(T_SWA)
}
