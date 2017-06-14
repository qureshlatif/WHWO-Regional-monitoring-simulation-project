################################################################################
#               Test of new code for having a user-specified grid              #
################################################################################
setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation")
library(raster)
library(rSPACE)
library(rgdal)

#### STEP 1: LOAD HABITAT AND GRID LAYERS

HabitatMap <- raster("E:/GISData/WHWO/Occ_sims/R6/Plot_exmpl/Habitat.tif",band=1)
HabitatMap[which(is.na(getValues(HabitatMap)))] <- 0
HabitatMap[which(getValues(HabitatMap)==0)] <- 0.99
GridMap<-raster("E:/GISData/WHWO/Occ_sims/R6/Plot_exmpl/survey.tif", band=1)
GridMap[which(is.na(getValues(GridMap)))] <- 0

# Toy version to run faster
library(spatstat)
library(maptools)

#### STEP 2: ENTER PARAMETER VALUES (I just made these up)

#-- Alternative dialogue input for parameters
# BaseParameters <- enter.parameters()

N <- round(0.2522*(7166.7/70))
trnd <- 1 #Can be 0.9, 0.95, 0.98, or 1.0
HR <- 1 #Home range radius (km; can be 1 or 0.6)

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

#  #-- Optional: specify subsetting levels for later use.  
#  BaseParameters<-c(BaseParameters,
#    list(
#      n_visit_test=c(2,4,6),                # Number of sampling occasions per year
#      detP_test   =c(0.8),                  # Per visit detection probability (if spp present)
#      grid_sample =c(0.05,0.15,0.25,0.35,   # Percent of grid to sample
#                    0.45,0.55,0.75,0.95),
#      sample_yrs  =0                        # Possible alternative models specifications
#    ))

set.seed(1)
#debug(encounter.history) # For exporting values to generate plots in ArcGIS.
encounter.history(map=HabitatMap, Parameters=BaseParameters, filter.map=GridMap, showSteps=T)
