###################################################################################################################
#               Code for simulating detection-nondetection data for 30 regional monitoring transects              #
###################################################################################################################
library(raster)
library(rSPACE)
library(rgdal)
library(spatstat)
library(maptools)

forests <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")
frst_area <- c(COWA=583861.9,DEOR=784698.1,FROR=1303163.5,GPWA=440452.1,MAOR=837544,MHOR=258748.1,OCOR=394800.5,
               OKWA=1319493.5,UMOR=733938.4,WAOR=1020271.1) #In hectares
trnd <- c(0.9,0.95,0.98,1)
HR <- 1 #Home range radius (km; can be 1 or 0.6)
D <- 0.2522 #Density = number of home ranges per 70 ha
HR_path <- "HR1000mD2522/" #subfolder corresponding with selected home range size and initial density
path <- paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/",HR_path,sep="")

library(R.utils)

for(i in 1:length(forests)) {
  pth <- paste(path,forests[i],sep="")
  if(!dir.exists(pth)) dir.create(pth)
  setwd(pth)
  N <- as.numeric(round(D*(frst_area[i]/70)))
  HabitatMap <- raster(paste("E:/GISData/WHWO/Occ_sims/R6/NF_",forests[i],"/habitat.tif",sep=""),band=1)
  GridMap<-raster(paste("E:/GISData/WHWO/Occ_sims/R6/NF_",forests[i],"/surveys.tif",sep=""), band=1)
  for(j in 1:length(trnd)) {
    pth <- paste(pth,"/Lmbd",trnd[j]*100,sep="")
    if(!dir.exists(pth)) dir.create(pth)
    setwd(pth)
    #-- Manual input version
    BaseParameters<-list(
      N               = N,               # Initial population size
      lmda            = trnd[j],              # Population growth rate
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

    if(!file.exists("Ppres.txt")) file.create("Ppres.txt",showWarnings=F)
  
    #untrace(encounter.history) #Run this if re-running within a session.
    trace(encounter.history, tracer=quote(if(file.exists("/Ppres.txt")){
      cat(round(P.pres,3),'\n', file="Ppres.txt", append=T)
      } ), at=25, where=rSPACE::create.landscapes)

    create.landscapes(n_runs=100, map=HabitatMap, Parameters=BaseParameters,
      filter.map=GridMap,skipConfirm=T, run.label=forests[i])
  }
}
