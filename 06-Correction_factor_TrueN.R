################################################################################################################
# True abundance at Gifford Pinchot NF (GPWA) was limited by suitable nesting habitat (rather than potential   #
# habitat), so it varied slightly with each simulation. This script tabulates a correction factor based on the #
# deviation from mean true abundance at GPWA to derive true abundance for each simulation under each trend     #
# scenario.                                                                                                    #
################################################################################################################

library(R.utils)
setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522")

dat100 <- read.table("GPWA/Lmbd100/GPWA/output/N_raw.txt",header=F,sep=" ")
dat98 <- read.table("GPWA/Lmbd98/GPWA/output/N_raw.txt",header=F,sep=" ")
dat95 <- read.table("GPWA/Lmbd95/GPWA/output/N_raw.txt",header=F,sep=" ")
dat90 <- read.table("GPWA/Lmbd90/GPWA/output/N_raw.txt",header=F,sep=" ")

N.mean <- round(mean(c(dat100$V1,dat98$V1,dat95$V1,dat90$V1))) # Mean abundance in year 1 at GPWA

C100 <- dat100$V1 - N.mean # Deviations from the mean for sims with lambda = 1
C98 <- dat98$V1 - N.mean # Deviations from the mean for sims with lambda = 0.98
C95 <- dat95$V1 - N.mean # Deviations from the mean for sims with lambda = 0.95
C90 <- dat90$V1 - N.mean # Deviations from the mean for sims with lambda = 0.9

out <- cbind(sim=seq(1,30),C90,C95,C98,C100) # Tabulate deviations from the mean for all sims
write.csv(out,"N_deviations_from_mean_by_sim.csv",row.names=F)
