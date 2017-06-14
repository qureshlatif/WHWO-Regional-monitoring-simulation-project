## Applies published habitat model (Latif et al. 2012, JWM) to generate habitat layers by forest ##

library(raster)
library(rgdal)
library(rgeos)
library(shapefiles)
library(maptools)
library(maps)
library(spatstat)
library(sp)
library(spdep)

source("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-WHWO/Qs/WHWO_nest_models/Maxent_HSI_functions.r")

NF.codes <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")
wrkspc <- "T:/FS/RD/RMRS/Science/WTE/Research/RMRS-WHWO/Qs/ArcGIS/Occ_sims/R6/"

for(i in 1:length(NF.codes)) {
  setwd(paste(wrkspc,"NF_",NF.codes[i],sep=""))
  slope <- raster("slope.tif",format="GTiff")
  cosasp <- raster("cosasp.tif",format="GTiff")
  loccc <- raster("ccov_1ha.tif",format="GTiff")
  lndcc <- raster("ccov_1km.tif",format="GTiff")
  pipo <- raster("pipo_1km.tif",format="GTiff")
  
  HSI <- hab <- slope
  lcc <- getValues(loccc)
  ldc <- getValues(lndcc)
  slp <- getValues(slope)
  ppo <- (getValues(pipo)/3409)*100
  csp <- getValues(cosasp)
  mxnt.scrs <- MxntSimplified_scr(lcc,ldc,csp,ppo,slp)
  hsi <- mxnt.scrs$hsi
  hsi[which(ppo<10)] <- 0
  HSI[1:length(HSI)] <- hsi
  writeRaster(HSI,filename="HSI.tif",format="GTiff",dataType="FLT4S",overwrite=T,update=T)
  hsi[which(hsi>=0.36)] <- 1
  hsi[which(hsi<0.36)] <- 0
  hsi[is.na(hsi)] <- 0
  hab[1:length(hab)] <- hsi
  writeRaster(hab,filename="habitat.tif",format="GTiff",dataType="FLT4S",overwrite=T,update=T)
}
