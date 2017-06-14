## Extracts and compiles layers identifying large-cone pine and ponderosa pine dominated forests ##
## Large-cone pine dominated forest defines "potential habitat" within which simulations were constrained. ##
## Ponderosa pine dominated forest informs application of nest habitat suitability model, which constrains home range centers ##

library(raster)
library(rgdal)
library(rgeos)
library(shapefiles)
library(maptools)
library(maps)
library(spatstat)
library(sp)
library(spdep)

setwd("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-WHWO/Qs/ArcGIS/Occ_sims/R6/interm")

# Get GNN forest type layer and attribute table #
#fortyp <- raster("fortyp_NF.tif",format="GTiff")
fortyp <- raster("fortyp_all.tif",format="GTiff")
key <- read.table("T:/FS/RD/RMRS/Science/WTE/Research/RMRS-WHWO/Qs/ArcGIS/gnn_sppsz_2014_04_21/fortyp_table.txt",header=T,sep=",",stringsAsFactors=F)

# Designate output raster #
lcpine <- ppine <- fortyp

# Identify raster values that correspond with large-cone pine dominated forest #
lcpine.codes <- c("PIPO","PILA","PIMO3","PIJE")
vals.lcpine <- c()
for(i in 1:length(lcpine.codes)) vals.lcpine <- c(vals.lcpine,grep(lcpine.codes[i],key$FORTYPBA))
vals.lcpine <- unique(vals.lcpine)
vals.lcpine <- key$VALUE[vals.lcpine]

# Generate large-cone pine binary values in output layer and write to file #
ind1 <- which(is.element(getValues(lcpine),vals.lcpine))
ind0 <- which(!is.element(getValues(lcpine),vals.lcpine))
lcpine[ind1] <- 1
lcpine[ind0] <- 0
#writeRaster(lcpine,filename="lcpine.tif",format="GTiff",overwrite=T)
writeRaster(lcpine,filename="lcpine_all.tif",format="GTiff",overwrite=T)
rm(lcpine)

# Identify raster values that correspond with ponderosa pine dominated forest #
vals.PIPO <- grep("PIPO",key$FORTYPBA)
vals.PIPO <- key$VALUE[vals.PIPO]

# Generate PIPO binary values in output layer and write to file #
ind1 <- which(is.element(getValues(ppine),vals.PIPO))
ind0 <- which(!is.element(getValues(ppine),vals.PIPO))
ppine[ind1] <- 1
ppine[ind0] <- 0
writeRaster(ppine,filename="PIPO.tif",format="GTiff",overwrite=T)
