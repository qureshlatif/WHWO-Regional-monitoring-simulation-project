###########################################################################################
#  Code for simulating detection-nondetection data for 30 regional monitoring transects   #
###########################################################################################
library(raster)
library(rSPACE)
library(rgdal)
library(spatstat)
library(maptools)
library(R.utils)

forest <- "OKWA"    # Run script for one NF at a time.
frst_area <- c(COWA=583861.9,DEOR=784698.1,FROR=1303163.5,GPWA=440452.1,MAOR=837544,MHOR=258748.1,OCOR=394800.5,
               OKWA=1319493.5,UMOR=733938.4,WAOR=1020271.1) #In hectares
trnd <- 0.9 #Can be 0.9, 0.95, 0.98, or 1.0
HR <- 1 #Home range radius (km; can be 1 or 0.6)
D <- 0.2522 #Density = number of home ranges per 70 ha
path <- paste("F:/...",HR_path,sep="")  # path where simulation results are to be stored
path.grid <- "E:/GISData/WHWO/Occ_sims/R6" # path where habitat and sampling grid layers are stored

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

#-- Set parameters for running rSPACE simulations
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

if(!file.exists("Ppres.txt")) file.create("Ppres.txt",showWarnings=F)

#as.list(body(encounter.history)) # To see where to put the trace (at=?), run this and place it at the [[#]] after 'P.pres[,tt] <- probPRES(useLayer, grid_layer)'
trace(encounter.history, tracer=quote(if(file.exists("Ppres.txt")){  # Places trace to store point-level encounter probabilities from each simulation
  cat(round(P.pres,3),'\n', file="Ppres.txt", append=T)
} ), at=28, where=rSPACE::createReplicates)

createReplicates(n_runs=30, map=HabitatMap, Parameters=BaseParameters,
                 filter.map=GridMap,skipConfirm=T, run.label=forest, base.name="rSPACEx", add=T)

### Additional scripts for compiling transect detection datasets available upon request.

########################################################################################################
# BUGS specification of constant-p occupancy model for analyzing repeat-survey transect detection data #
########################################################################################################

model { 
  # prior distributions
  p~dunif(0,1) # Detection probability (point scale)
  
  for (t in 1:n.yrs) {
    B0[t] ~ dnorm(0,0.1)T(-10,10) # Transect scale logit occupancy for points
    logit(PSI[t]) <- B0[t]
  }
  
  for (i in 1:n.trns) {
    for (t in 1:n.yrs) {
      Z[i,t] ~ dbin(PSI[t],1) # Occupancy state at transects following year 1
      prob.y[i,t] <- Z[i,t]*p     
      Y.mat[i,t] ~ dbin(prob.y[i,t],n.vsts[i,t]) # OBSERVATION MODEL
    }
  }
}

######################################################################################################
# BUGS specification of yearly-p occupancy model for analyzing repeat-survey transect detection data #
######################################################################################################

model { 
  # prior distributions
  for (t in 1:n.yrs) {
    B0[t] ~ dnorm(0,0.1)T(-10,10) # Transect scale logit occupancy for points
    logit(PSI[t]) <- B0[t]
    p[t]~dunif(0,1)
  }
  
  for (i in 1:n.trns) {
    for (t in 1:n.yrs) {
      Z[i,t] ~ dbin(PSI[t],1) # Occupancy state at transects following year 1
      prob.y[i,t] <- Z[i,t]*p[t]     
      Y.mat[i,t] ~ dbin(prob.y[i,t],n.vst) # OBSERVATION MODEL
    }
  }
}

#######################################################################################################
# BUGS specification of logistic regression model for analyzing single-survey transect detection data #
#######################################################################################################

model { 
  # prior distributions
  for (t in 1:n.yrs) {
    B0[t] ~ dnorm(0,0.1)T(-10,10) # Transect scale logit occupancy for points
    logit(PSI[t]) <- B0[t]
  }
  
  for (i in 1:n.trns) {
    for (t in 1:n.yrs) {
      Y.lgreg[i,t] ~ dbin(PSI[t],1) # OBSERVATION MODEL
    }
  }
}


### Additional scripts for implementing analyses of simulated data and storing results available upon request. ###