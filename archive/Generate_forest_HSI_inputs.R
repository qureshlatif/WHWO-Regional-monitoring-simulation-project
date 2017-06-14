## Looks like initially intended to generate input layers for applying HSI model, but doesn't look finished. Must not have been used. ##
## Archiving for now ##

setwd("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-WHWO/Qs/ArcGIS/Occ_sims/R6/")

library(raster)
library(rgdal)
library(rgeos)
library(shapefiles)
library(maptools)
library(maps)
library(spatstat)
library(sp)
library(spdep)

# Import NF polygons.
forests <- readShapeSpatial("Eastern_R6_forests_n83z10.shp")
proj4string(forests) <-CRS("+proj=utm +zone=10 +datum=NAD83")
for.names <- as.character(forests@data$FORESTNAME)
for.abbrv <- c("FROR","GPWA","OKWA","MAOR","OCOR","DEOR","WAOR","COWA","UMOR","MHOR")

# Get sampling points and raw layers for calculating habitat
samp_pnts <- readShapeSpatial("R6_sampling_grid.shp")
slope <- raster("slope.tif",format="GTiff")
cosasp <- raster("cosasp.tif",format="GTiff")
loccc <- raster("ccov_1ha.tif",format="GTiff")
lndcc <- raster("ccov_1km.tif",format="GTiff")
pipo <- raster("pipo_1km.tif",format="GTiff")

for(i in 1:nrow(forests@data)) {
  frst <- forests[i,]
  # 2. Apply 1km buffer around each NF polygon.
  fbuff <- gBuffer(frst,width=1000)
  #___ Not sure how to write buffer to shapefile. Prob don't need to though___#
  #fbuff <- SpatialPolygonsDataFrame(fbuff,data=frst@data)
  #dir.create(paste(getwd(),"/",for.abbrv[i],sep=""))
  #writeOGR(fbuff,dsn=paste(getwd(),"/",for.abbrv[i],sep=""),layer="frst_bound_2km",driver="ESRI Shapefile")
  #___________________________________________________________________________#
  slp <- extract(slope,fbuff)
  casp <- extract(cosasp,fbuff)
}
