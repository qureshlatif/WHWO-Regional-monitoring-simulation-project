#########################################################################
# Objective 2 scenario with yearly-p model, 3-pt transect, panel design #
#########################################################################

# Set working directory to folder with simulated datasets:
library(foreign)
require(data.table)
library(R.utils)

#### Variables ############################################
trnd <- 0.98 # Can be 1, 0.98, 0.95, or 0.9
trnd.psi <- 0.9902002
                      # (averaged across 30 simulations for each trend scenario)
sims <- 100 #Set to number of simulated analyses desired per scenario.
pops <- 30 # The number of populations and point detection datasets generated per forest.
n.vst <- 2 #Number of surveys conducted at each transect per season.
n.trns <- 600 #Select from 60, 90, 120, 150
n.pnts <- 3 #No. surveyed points per transect (Set to 5 or 10)
n.yrs <- 20 #Number of years (Currently fixed).
D <- 0.2522 #Initial territory density (number per 70 ha averaged across National Forest in year 1)
HR <- 1000 #Home range radius in meters.
mon.freq <- 3 #Select from 1, 2, or 3. Indicates whether all (1), half (2), or a third (3) of transects surveyed each year.
###########################################################

parent.folder <- paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR",
              HR,"mD",D*10000,sep="")
parent.GIS <- "E:/GISData/WHWO/Occ_sims/R6"

source("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/scripts/Functions.R")

# Set working path specific to each scenario
path <- paste(parent.folder,"/obj2/yp",sep="")
if(!dir.exists(path)) dir.create(path)
path <- paste(path,"/Tr",n.trns,"Pt",n.pnts,"Mf",mon.freq,sep="")
if(!dir.exists(path)) dir.create(path)

forest <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")
frst_area <- c(COWA=583861.9,DEOR=784698.1,FROR=1303163.5,GPWA=440452.1,MAOR=837544,MHOR=258748.1,OCOR=394800.5,
               OKWA=1319493.5,UMOR=733938.4,WAOR=1020271.1) #In hectares

n.yrs <- 20 #Number of years will now be fixed.

pops.draw <- sample(pops, sims, replace = T)
dat.files <- paste("rSPACEx",pops.draw,".txt",sep="")

###### Analysis using single-scale occupancy model with constant p ###################################

model_yp <- "F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/scripts/Model_jags_yearly_p.R"

nc <- 4
nb <- 5000
ni <- 20000
nt <- 1

library(R2jags)
# Number of transects established in each forest #
MaxN <- Maxtrns(parent.GIS,forest) # Get max number of transects per forest
frst_n <- round((frst_area/sum(frst_area))*n.trns) 
if(sum(frst_n)!=n.trns) { #Ensures sum of number of transects established in each forest = n.trns
  diff <- n.trns-sum(frst_n)
  ind <- sample.int(10,size=(abs(diff)),prob=frst_area/sum(frst_area))
  for(j in 1:length(ind)) frst_n[ind[j]] <- frst_n[ind[j]] + (abs(diff)/diff)
}
if(any((MaxN-frst_n)<0)) {
  ind <- which((MaxN-frst_n)<0)
  count <- 0
  for(j in ind) {
    count <- count + as.numeric(frst_n[j]-MaxN[j])
    frst_n[j] <- MaxN[j]
  }
  ind <- (1:length(frst_n))[-ind]
  ind <- sample(ind,size=count,prob=frst_area[ind]/sum(frst_area[ind]),replace=T)
  for(j in 1:length(ind)) frst_n[ind[j]] <- frst_n[ind[j]] + 1
}
#________________________________________________#

for(i in 1:sims) {
  # Compile data for analysis
  Y.arry <- array(NA,dim=c(n.trns,10,2,20)) #dim = transects, max points, max visits, max years
  for(f in 1:length(forest)) {
    frst <- forest[f]
    n <- as.numeric(frst_n[f])
    sets <- read.table(paste(parent.GIS,"/NF_",frst,"/","Transect_sets.txt",sep=""),header=T,sep=",")
    
    #Compile table that matches transect with point IDs
    tp_match <- fread(paste(parent.GIS,"/NF_",frst,"/","Survey_pnt_coords.txt",sep=""),header=T,sep=",")
    tp_match <- tp_match[,.(POINTID,tr=c(H1tr,H2tr,V1tr,V2tr),pt=c(H1pt,H2pt,V1pt,V2pt))]
    
    #Randomly select transects from NF pool and identify survey points along those transects
    ntrns.by.set <- apply(sets[,-1],2,sum)
    set.col <- sample((2:ncol(sets))[which(ntrns.by.set>=n)],1)
    trns.available <- which(sets[,set.col]==1)
    transects <- sets$transects[sample(trns.available,n)]
    pnts <- matrix(NA,nrow=length(transects),ncol=10)
    for(tr in 1:length(transects)) pnts[tr,] <- tp_match$POINTID[which(tp_match$tr==transects[tr])]

    #Retrieve data for points along selected transects and add to Y.arry
    pth <- paste(parent.folder,"/",frst,"/Lmbd",trnd*100,"/",frst,sep="")
    dat <- fread(paste(pth,"/",dat.files[i],sep=""),header=F,sep='*')
    dat <- dat[V2 %in% as.numeric(pnts),]
    ID <- dat$V2
    dat <- substr(dat[,V3],3,(n.vst*20+2))
    dat <- strsplit(dat,split='')
    y <- array(NA,dim=c(length(transects),ncol(pnts),2,20))
    for(tr in 1:dim(y)[1]) for(pt in 1:dim(y)[2]) {
      ind <- which(ID==pnts[tr,pt])
      obs <- as.numeric(dat[[ind]])
      y[tr,pt,,] <- matrix(obs,nrow=n.vst,ncol=20)
    }
    if(f==1) Y.arry[1:frst_n[1],,,] <- y
    if(f>1) Y.arry[(sum(frst_n[1:(f-1)])+1):sum(frst_n[1:f]),,,] <- y
  }

  # Condense data for analysis
  Y.mat <- apply(Y.arry[,1:n.pnts,,],c(1,3,4),max)
  Y.mat <- apply(Y.mat, c(1, 3), sum,na.rm=T)
  
  if(mon.freq==2) {
    Y.mat[(1:(nrow(Y.mat)/2)),seq(2,20,by=2)] <- NA
    Y.mat[((nrow(Y.mat)/2)+1):nrow(Y.mat),seq(1,19,by=2)] <- NA
  }
  if(mon.freq==3) {
    Y.mat[(1:(nrow(Y.mat)/3)),seq(3,18,by=3)] <- NA
    Y.mat[((nrow(Y.mat)/3)+1):((nrow(Y.mat)/3)*2),seq(2,20,by=3)] <- NA
    Y.mat[(((nrow(Y.mat)/3)*2)+1):nrow(Y.mat),seq(1,19,by=3)] <- NA
  }
  
  #Yearly p analysis
  data <- list("Y.mat","n.trns","n.yrs","n.vst")
  z.init <- (Y.mat>0)*1
  inits <- function()
    list (Z=z.init)
  parameters <- c("p","PSI")
  out<-jags(data,inits,parameters,model_yp,n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni)
  if(any(round(out$BUGSoutput$summary[,"Rhat"],digits=1)>=1.1)|any(out$BUGSoutput$summary[,"n.eff"]<100))
    ni.new <- ni
  while(any(round(out$BUGSoutput$summary[,"Rhat"],digits=1)>=1.1)|any(out$BUGSoutput$summary[,"n.eff"]<100)){
    ni.new <- ni.new*2
    out<-update(out,n.iter=ni.new)
  }
  saveObject(out,paste(path,"/model_yp_sim",i,sep="")) 
}

# Compile parameter estimates #
library(R.utils)
library(abind)
library(dplyr)
library(stringr)

PSIprd_yp.sim <- pprd_yp.sim <- array(NA,dim=c(sims,n.yrs,3))
PSIprd_yp.trend <- matrix(NA,sims,5)
for(i in 1:sims) {
  mod_yp <- loadObject(paste(path,"/model_yp_sim",i,sep="")) 

  # Store occupancy estimates from yearly-p model
  PSI <- mod_yp$BUGSoutput$sims.list$PSI
  PSIprd_yp.sim[i,,1] <- apply(PSI,2,median)
  PSIprd_yp.sim[i,,2] <- apply(PSI,2,function(x) quantile(x,prob=0.025,type=8))
  PSIprd_yp.sim[i,,3] <- apply(PSI,2,function(x) quantile(x,prob=0.975,type=8))
  PSI[which(PSI==1)] <- 0.9999
  PSI[which(PSI==0)] <- 0.0001
  assign(paste0("PSI_yp",str_pad(i, width = 3, pad = "0")), PSI)
  
  # Store detectability estimates from yearly-p model
  p <- mod_yp$BUGSoutput$sims.list$p
  pprd_yp.sim[i,,1] <- apply(p,2,median)
  pprd_yp.sim[i,,2] <- apply(p,2,function(x) quantile(x,prob=0.025,type=8))
  pprd_yp.sim[i,,3] <- apply(p,2,function(x) quantile(x,prob=0.975,type=8))
}

# Trend estimates for yearly-p model
PSI <- mget(paste0("PSI_yp", str_pad(1:sims, width = 3, pad = "0")))
rm(list = paste0("PSI_yp", str_pad(1:sims, width = 3, pad = "0")))
PSIprd.tnd <- lapply(PSI, function(x) apply(x, 1, trend))

PSIprd_yp.trend[,1] <- sapply(PSIprd.tnd, median)
PSIprd_yp.trend[,2] <- sapply(PSIprd.tnd, function(x) quantile(x,prob=0.025,type=8))
PSIprd_yp.trend[,3] <- sapply(PSIprd.tnd, function(x) quantile(x,prob=0.975,type=8))
PSIprd_yp.trend[,4] <- sapply(PSIprd.tnd, function(x) mean((exp(x) - trnd)^2))
PSIprd_yp.trend[,5] <- sapply(PSIprd.tnd, function(x) mean((exp(x) - trnd.psi)^2))

dimnames(PSIprd_yp.trend)[[2]] <- c("median","lo95","hi95","MSE","MSE.psi")

rm(PSIprd.tnd,ID,Y.arry,dat,data,f,i,ind,frst_n,frst_area,frst,mod_yp,
   n,nb,nc,nt,ni,obs,out,parameters,pt,pth,set.col,tr,transects,trns.available,y,
   Y.mat,PSI,p,pnts,sets,tp_match,z.init)
save.image(paste(path,"Results.RData",sep=""))
####################################################################################################################
