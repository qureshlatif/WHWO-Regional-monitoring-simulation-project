## Looks like script for generating sampling grid when simulations were centered around established transects. Obsolete.
## Archived Nov 16, 2016.

setwd("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/")

library(raster)
library(rgdal)
library(shapefiles)
library(maptools)

## Create point shapefile for individual transects ##
#r6_pnts <- readShapeSpatial("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/R6/R6_points_2011.shp")
#TR.names <- unique(as.character(r6_pnts$Transect_I))

#for(i in 1:length(TR.names)) {
#  path <- paste("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/",TR.names[i],"/",sep="")
#  dir.create(path)
#  pull <- r6_pnts[r6_pnts$Transect_I==TR.names[i],]
#  writeSpatialShape(pull,paste(path,"points",sep=""))
#}

r6_pnts <- readShapeSpatial("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/R6/R6_points_2011.shp")
TR.names <- unique(as.character(r6_pnts$Transect_I))
TR.or <- TR.names[-which(is.element(TR.names,c("TR601","TR602","TR603","TR604")))]
rm(r6_pnts)

library(spatstat)
library(sp)
library(raster)

## Extract habitat rasters for Oregon transects from Maxent layer ##
#HSI <- raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/Mxnt_habitat_or_utm10.tif",band=1)
#for(tr in 1:length(TR.or)) {
#  pnts <- readShapeSpatial(paste(TR.or[tr],"/points.shp",sep=""))
#  ## Generate 5km buffer on points and apply as mask to Maxent raster to generate habitat layer. ##
#  polys <- list()
#  for (i in 1:nrow(pnts)) {
#    discbuff <- disc(radius=5000, centre=c(pnts$coords_x1[i],pnts$coords_x2[i]))
#    discpoly <- Polygon(rbind(cbind(discbuff$bdry[[1]]$x,y=discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1],
#                                                                                    y=discbuff$bdry[[1]]$y[1])))
#    polys<-c(polys, discpoly)
#  }
#  spolys<-list()
#  for(i in 1:length(polys)) {
#    spolybuff<-Polygons(list(polys[[i]]), ID=row.names(pnts)[i])
#    spolys<-c(spolys, spolybuff)
#  }
#  polybuff5km<-SpatialPolygons(spolys)
#
#  ID <- rep(1,length(polybuff5km))
#  poly5km <- unionSpatialPolygons(polybuff5km,IDs=ID)
#  
#  hab <- crop(HSI,poly5km)
#  hab <- mask(hab,poly5km)
#  hab[is.na(hab)] <- 0.99
#  hab[which(values(hab)==0)] <- 0.99
#  
#  writeRaster(hab,filename=paste("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/",TR.or[tr],"/hab.tif",sep=""),format="GTiff",overwrite=T)  
#}

### Generate Maxent HSI layer for WA transects ###
## Convert forest type to pipo indicator layer for WA ##
#fortyp <- raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/fortyp.tif",band=1)
#pipo <- fortyp
#vals <- pipo@data@attributes[[1]]$ID[grep("PIPO",pipo@data@attributes[[1]]$FORTYPBA)]
#pipo <- calc(pipo,fun=function(x) ifelse(is.element(x,vals),1,0))
#writeRaster(pipo,filename="T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/pipo.tif",format="GTiff",overwrite=T)

#source("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/WHWO_nest_models/Maxent_HSI_functions.r")
#LocCC <- raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/loccc.tif",band=1)
#LandCC<-raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/landcc.tif",band=1)
#PIPO<-raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/pipo_1km.tif",band=1)
#Slp<-raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/slope.tif",band=1)
#cosasp<-raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/cosasp.tif",band=1)

#Ldc <- values(LandCC)
#lcc <- values(LocCC)
#ppo <- (values(PIPO)/3409)*100
#slp <- values(Slp)
#casp <- values(cosasp)
#mxnt.scrs <- MxntSimplified_scr(lcc,Ldc,casp,ppo,slp)
#hsi <- as.numeric(mxnt.scrs$hsi>=0.36)
#HSI <- LocCC
#values(HSI) <- hsi
#writeRaster(HSI,filename="T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/Mxnt36.tif",format="GTiff",overwrite=T)

##################################################

## Extract habitat rasters for Washington transects from Maxent layer ##
#HSI <- raster("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/WA/Mxnt36.tif",band=1)
#TR.WA <- c("TR601","TR602","TR603","TR604")
#for(tr in 1:length(TR.WA)) {
#  pnts <- readShapeSpatial(paste(TR.WA[tr],"/points.shp",sep=""))
#  ## Generate 5km buffer on points and apply as mask to Maxent raster to generate habitat layer. ##
#  polys <- list()
#  for (i in 1:nrow(pnts)) {
#    discbuff <- disc(radius=5000, centre=c(pnts$coords_x1[i],pnts$coords_x2[i]))
#    discpoly <- Polygon(rbind(cbind(discbuff$bdry[[1]]$x,y=discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1],
#                                                                                    y=discbuff$bdry[[1]]$y[1])))
#    polys<-c(polys, discpoly)
#  }
#  spolys<-list()
#  for(i in 1:length(polys)) {
#    spolybuff<-Polygons(list(polys[[i]]), ID=row.names(pnts)[i])
#    spolys<-c(spolys, spolybuff)
#  }
#  polybuff5km<-SpatialPolygons(spolys)
#
#  ID <- rep(1,length(polybuff5km))
#  poly5km <- unionSpatialPolygons(polybuff5km,IDs=ID)
#
#  hab <- crop(HSI,poly5km)
#  hab <- mask(hab,poly5km)
#  hab[is.na(hab)] <- 0.99
#  hab[which(values(hab)==0)] <- 0.99
#
#  writeRaster(hab,filename=paste("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/",TR.WA[tr],"/hab.tif",sep=""),format="GTiff",overwrite=T)  
#}


## Generate 150m buffer on points, convert to raster, and export for survey raster layer ##
for(tr in 1:length(TR.names)) {
  pnts <- readShapeSpatial(paste(TR.names[tr],"/points.shp",sep=""))
  hab <- raster(paste("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/",TR.names[tr],"/hab.tif",sep=""),band=1)
  PID <- paste(pnts$Transect_I,"_0",pnts$Point_ID,sep="")
  PID[which(pnts$Point_ID==10)] <-
    paste(pnts$Transect_I[which(pnts$Point_ID==10)],"_",pnts$Point_ID[which(pnts$Point_ID==10)],sep="")
  row.names(pnts) <- PID
  polys <- list()
  for (i in 1:nrow(pnts)) {
    discbuff <- disc(radius=150, centre=c(pnts$coords_x1[i],pnts$coords_x2[i]))
    discpoly <- Polygon(rbind(cbind(discbuff$bdry[[1]]$x,y=discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1],
                                                                                    y=discbuff$bdry[[1]]$y[1])))
    polys<-c(polys, discpoly)
  }
  spolys<-list()
  for(i in 1:length(polys)) {
    spolybuff<-Polygons(list(polys[[i]]), ID=PID[i])
    spolys<-c(spolys, spolybuff)
  }
  polybuff150m <- SpatialPolygons(spolys)
  polybuff150m <- SpatialPolygonsDataFrame(polybuff150m,data=pnts@data)
  bg <- hab
  values(bg) <- 0
  survey <- rasterize(polybuff150m,bg,background=0)
  writeRaster(survey,filename=paste("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-R6-WHWO/Qs/ArcGIS/Occ_sims/R6/",TR.names[tr],"/surveys.tif",sep=""),format="GTiff",overwrite=T)
}

