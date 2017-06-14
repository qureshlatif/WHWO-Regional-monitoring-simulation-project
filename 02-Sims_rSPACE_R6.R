###################################################################################################################
#               Code for simulating detection-nondetection data for 30 regional monitoring transects              #
###################################################################################################################
library(raster)
library(rSPACE)
library(rgdal)
library(spatstat)
library(maptools)

forest <- "OKWA"
frst_area <- c(COWA=583861.9,DEOR=784698.1,FROR=1303163.5,GPWA=440452.1,MAOR=837544,MHOR=258748.1,OCOR=394800.5,
               OKWA=1319493.5,UMOR=733938.4,WAOR=1020271.1) #In hectares
trnd <- 0.9 #Can be 0.9, 0.95, 0.98, or 1.0
HR <- 1 #Home range radius (km; can be 1 or 0.6)
D <- 0.2522 #Density = number of home ranges per 70 ha
HR_path <- "HR1000mD2522/" #subfolder corresponding with selected home range size and initial density
path <- paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/",HR_path,sep="")
#path <- paste("C:/qs_files/occ_sims/",HR_path,sep="")
path.grid <- "E:/GISData/WHWO/Occ_sims/R6"

library(R.utils)

pth <- paste(path,forest,sep="")
if(!dir.exists(pth)) dir.create(pth)
setwd(pth)
N <- as.numeric(round(D*(frst_area[forest]/70)))
HabitatMap <- raster(paste(path.grid,"/NF_",forest,"/habitat.tif",sep=""),band=1)
HabitatMap[which(getValues(HabitatMap)==0)] <-
  0.99 # Replace 0s with 0.99 so that WHWO movement is not restricted by habitat (only HR center should be restricted)
GridMap<-raster(paste(path.grid,"/NF_",forest,"/surveys.tif",sep=""), band=1)
pth <- paste(pth,"/Lmbd",trnd*100,sep="")
if(!dir.exists(pth)) dir.create(pth)
setwd(pth)
#-- Manual input version
BaseParameters<-list(
  N               = N,               # Initial population size
  trendtype       ="abundance-exponential", #Not a necessary addition, but we're moving towards requiring this to be set explicitly.  You can also use "abundance-linear".
  lmda            = trnd,              # Population growth rate
  n_yrs           = 20,                # Maximum number of years in simulation
  n_visits        = 2,                  # Maximum number of visits per year
  grid_size       = 25,        ##IGNORED!         # Cell size in grid
  MFratio         = c(1),               # Ratio of types of individuals
  buffer          = c(0.3),               # Distance between individual center locations (km)
  moveDist        = c(HR),               # Movement radius
  moveDistQ       = c(0.95),             # Proportion of time in radius (defines meaning of "how far")
  maxDistQ        = c(0),               # Truncate movements above 1 SD
  habitat.cutoff  = 1,                  # Minimum habitat value required for individual center locations
  sample.cutoff   = 0.5,      ##IGNORED!          # % pixels in cell above habitat value to include cell in sample frame
  repeat.groups   = T,                 #Doubles availability to represent pair of individuals
  wghts         = T)         #Indicates whether the habitat suitability values should be used to weight the probability of choosing  individual center locations

#from <- paste(forest,"/rSPACEy",seq(1,8),".txt",sep="")
#to <- paste(forest,"/rSPACEx",seq(1,8)+92,".txt",sep="")
#file.rename(from,to)

if(!file.exists("Ppres.txt")) file.create("Ppres.txt",showWarnings=F)

#as.list(body(encounter.history)) # To see where to put the trace (at=?), run this and place it at the [[#]] after 'P.pres[,tt] <- probPRES(useLayer, grid_layer)'

#untrace(encounter.history) #Run this if re-running within a session.
trace(encounter.history, tracer=quote(if(file.exists("Ppres.txt")){
  cat(round(P.pres,3),'\n', file="Ppres.txt", append=T)
} ), at=28, where=rSPACE::createReplicates)

createReplicates(n_runs=30, map=HabitatMap, Parameters=BaseParameters,
                 filter.map=GridMap,skipConfirm=T, run.label=forest, base.name="rSPACEx", add=T)
