# Set working directory to folder with simulated datasets:
library(foreign)
library(R.utils)

#### Variables ############################################
D <- 0.2522 #Initial territory density (number per 70 ha averaged across National Forest in year 1)
HR <- 1000 #Home range radius in meters.
###########################################################


forest <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")

for(f in 1:length(forest)) {
  frst <- forest[f]  
  setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR",
              HR,"mD",D*10000,"/",frst,"/Lmbd90/",frst,sep=""))
  dat <- read.delim("rSPACEx1.txt",header=F,sep='*',as.is=T)
  PID <- dat$V2
  saveObject(PID,paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR",
              HR,"mD",D*10000,"/",frst,"/PIDs_in_sampled_order",sep=""))
}
