###################################################################################################################
#               Code for simulating detection-nondetection data for 30 regional monitoring transects              #
###################################################################################################################
library(raster)
library(rSPACE)
library(rgdal)
library(spatstat)
library(maptools)
library(R.utils)

TR.names <- c("TR601","TR602","TR603","TR604","TR605","TR606","TR607","TR608","TR609","TR611","TR612",
              "TR613","TR614","TR615","TR616","TR617","TR618","TR619","TR620","TR621","TR623","TR624","TR625","TR626",
              "TR628","TR629","TR630","TR640","TR646","TR647")

#### Variables ####
trnd <- 1
N <- 36 #seq(10,90,by=4)
HR <- 1 #home range radius in kilometers
###################

for(n in N) {
  ##Name working directory 
  path <- paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR",HR*1000,"m_yr1/constN",n,sep="")
  dir.create(path)
  setwd(path)

  for(tr in 1:length(TR.names)) { #For original transects
    HabitatMap <- raster(paste("E:/GISData/WHWO/Occ_sims/R6/archive/",TR.names[tr],"/hab.tif",sep=""),band=1)
    GridMap<-raster(paste("E:/GISData/WHWO/Occ_sims/R6/archive/",TR.names[tr],"/surveys.tif",sep=""), band=1)
  
    #-- Manual input version
    BaseParameters<-list(
      N               = n,               # Initial population size
      trendtype       ="abundance-exponential", #Not a necessary addition, but we're moving towards requiring this to be set explicitly.  You can also use "abundance-linear".
      lmda            = trnd,              # Population growth rate
      n_yrs           = 1,                # Maximum number of years in simulation
      n_visits        = 2,                  # Maximum number of visits per year
      grid_size       = 25,        ##IGNORED!         # Cell size in grid
      MFratio         = c(1),               # Ratio of types of individuals
      buffer          = c(0.3),               # Distance between individual center locations (km)
      moveDist        = c(HR),               # Movement radius (km)
      moveDistQ       = c(0.95),             # Proportion of time in radius (defines meaning of "how far")
      maxDistQ        = c(0),               # Truncate movements above 1 SD
      habitat.cutoff  = 1,                  # Minimum habitat value required for individual center locations
      sample.cutoff   = 0.5,      ##IGNORED!          # % pixels in cell above habitat value to include cell in sample frame
      repeat.groups   = T,                 #Doubles availability to represent pair of individuals
      wghts         = T)         #Indicates whether the habitat suitability values should be used to weight the probability of choosing  individual center locations

    createReplicates(n_runs=30, map=HabitatMap, Parameters=BaseParameters,filter.map=GridMap,
      skipConfirm=T, run.label=TR.names[tr])
  }
}
