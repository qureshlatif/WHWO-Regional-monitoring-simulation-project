##*** Redundant with 'Generate_forest_sampling_grid' script.***##
## Consolidated and archived Nov 16, 2016 ##

##########################################################################################
# Identifies potential survey points (i.e., those capable of forming a transect) and     #
# eliminates all others (points not belonging to any set of 10 continuous line of points #
# oriented either north-south or east-west                                               #
##########################################################################################

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

## Identify and store grid point coordinates ##
for(i in 1:length(NF_codes)) {
  wrk <- paste(wrksp,"NF_",NF_codes[i],sep="")
  setwd(wrk)
  pnts <- readShapeSpatial("Sampling_pnts.shp")
  mark.dlt.x <- mark.dlt.y <- numeric(length=length(pnts))
  for(j in 1:length(mark.dlt.x)) {
    xy <- pnts@data[j,c("X","Y")]
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
  pnts <- pnts[-which(mark.dlt.x==1&mark.dlt.y==1),]
  write.table(pnts@data,"Survey_pnt_coords.txt",row.names=F,sep=",")
}
  