## Calculates true abundance for each simulation after applying correction factor for GPWA.                           ##
## GPWA was the only forest for which true N was constrained below the prescribed max density by habitat availability ##

# Set working directory to folder with simulated datasets:
library(foreign)
library(R.utils)

#### Variables ############################################
trnd <- 1 # Can be 1, 0.98, 0.95, or 0.9
D <- 0.2522 #Initial territory density (number per 70 ha averaged across National Forest in year 1)
HR <- 1000 #Home range radius in meters.
###########################################################

n.yrs <- 20 # Fixed for now.

parent.folder <- paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR",
              HR,"mD",D*10000,sep="")

# Set working path specific to each scenario
path <- paste(parent.folder,"/obj1/lmb",trnd*100,sep="")
if(!dir.exists(path)) dir.create(path)

forest <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")
trueN.Y1 <- sum(c(2104,2755,4695,1340,3018,932,1422,4754,2644,3676)) # Copied from 'output/N_final.txt' files.

# Apply correction factor to get trueN.Y1 for each sim
CF.tab <- read.csv(paste0(parent.folder,"/N_deviations_from_mean_by_sim.csv"),header=T)
CF <- CF.tab[,paste("C",trnd*100,sep="")]
trueN.Y1 <- trueN.Y1+CF
trueN <- matrix(NA,nrow=length(trueN.Y1),ncol=n.yrs)
for(i in 0:(n.yrs-1)) trueN[,i+1] <- trueN.Y1*(trnd^i)
trueN <- round(trueN)

saveObject(trueN,paste(path,"/TrueN",sep=""))
