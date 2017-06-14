## Generates sampling grid input layer for rSPACE by forest ##

library(raster)
library(rgdal)
library(rgeos)
library(shapefiles)
library(maptools)
library(maps)
library(spatstat)
library(sp)
library(spdep)

wrksp <- "T:/FS/RD/RMRS/Science/WTE/Research/RMRS-WHWO/Qs/ArcGIS/Occ_sims/R6/"
NF_codes <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")

## Generate grid rasters for simulations ##
for(i in 1:length(NF_codes)) {
  wrk <- paste(wrksp,"NF_",NF_codes[i],sep="")
  setwd(wrk)
  pnts <- readShapeSpatial("Sampling_pnts.shp") # Reads in shapefile with 300m point grid 
  mark.dlt.x <- mark.dlt.y <- numeric(length=length(pnts))
  for(j in 1:length(mark.dlt.x)) {  # Find all points that can't form a transect (continuous line of 10 points) and mark for deletion #
    xy <- pnts@data[j,c("X","Y")]
    
    # First mark all points incapable for forming a transect oriented east-west
    chck.x <- pnts@data[-j,c("X","Y")]
    chck.x <- chck.x[which(chck.x$Y==xy$Y),]
    if(nrow(chck.x)>0) chck.x <- chck.x[which(chck.x$X<=(xy$X+2700)&chck.x$X>=(xy$X-2700)),]
    ifelse(nrow(chck.x)<9,keep.x <- FALSE,keep.x <- TRUE)
    if(keep.x==TRUE) {  
      all.x <- seq(xy$X-2700,xy$X+2700,by=300)
      x.present <- is.element(all.x,c(chck.x$X,xy$X))
      keep.x <- any(sum(x.present[1:10])==10,sum(x.present[2:11])==10,sum(x.present[3:12])==10,sum(x.present[4:13])==10,
                    sum(x.present[5:14])==10,sum(x.present[6:15])==10,sum(x.present[7:16])==10,sum(x.present[8:17])==10,
                    sum(x.present[9:18])==10,sum(x.present[10:19])==10)
    }
    if(keep.x==F) mark.dlt.x[j] <- 1
    
    # Then mark all points incapable for forming a transect oriented north-south
    chck.y <- pnts@data[-j,c("X","Y")]
    chck.y <- chck.y[which(chck.y$X==xy$X),]
    if(nrow(chck.y)>0) chck.y <- chck.y[which(chck.y$Y<=(xy$Y+2700)&chck.y$Y>=(xy$Y-2700)),]
    ifelse(nrow(chck.y)<9,keep.y <- FALSE,keep.y <- TRUE)
    if(keep.y==TRUE) {
      all.y <- seq(xy$Y-2700,xy$Y+2700,by=300)
      y.present <- is.element(all.y,c(chck.y$Y,xy$Y))
      keep.y <- any(sum(y.present[1:10])==10,sum(y.present[2:11])==10,sum(y.present[3:12])==10,sum(y.present[4:13])==10,
                    sum(y.present[5:14])==10,sum(y.present[6:15])==10,sum(y.present[7:16])==10,sum(y.present[8:17])==10,
                    sum(y.present[9:18])==10,sum(y.present[10:19])==10)
    }
    if(keep.y==F) mark.dlt.y[j] <- 1
  }
  pnts <- pnts[-which(mark.dlt.x==1&mark.dlt.y==1),]  # Remove points that can't form either a north-south or east-west transect.
  write.table(pnts@data,paste0(getwd(),"/Survey_pnt_coords.txt"),row.names=F,sep=",")  # Write coordinates to table for producing point shapefile
  survey.area <- gBuffer(pnts,width=150,quadsegs=20,byid=T,id=pnts@data$POINTID)  # Generate 150m buffer around remaining points to define sampling grid for rSPACE
  hab <- raster(paste(getwd(),"/habitat.tif",sep=""))  # Get habitat layer to serve as rasterization template
  surv.rast <- rasterize(survey.area, hab, field=survey.area@data$POINTID)  # Convert sampling grid to raster format
  surv.rast[is.na(surv.rast)] <- 0
  writeRaster(surv.rast,paste0(getwd(),"/surveys.tif"),format="GTiff",overwrite=T) # Save sampling grid raster to file.
}
