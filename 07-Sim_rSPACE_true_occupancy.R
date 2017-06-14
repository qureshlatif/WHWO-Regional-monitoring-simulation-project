# Set working directory to folder with simulated datasets:
require(foreign)
require(data.table)
require(R.utils)
require(dplyr)

#### Variables ############################################
trnd <- 0.95 # Can be 1, 0.98, 0.95, or 0.9
nsim <- 30 # Set to total number of simulated populations (or landscapes).
###########################################################

#### Fixed for now, but could change in the future ########
D <- 0.2522 #Initial territory density (number per 70 ha averaged across National Forest in year 1)
HR <- 1000 #Home range radius in meters.
###########################################################

# Directories set for running on server space. #
#parent.folder <- paste("C:/qs_files/occ_sims/HR",
#              HR,"mD",D*10000,sep="")
#parent.GIS <- "C:/qs_files/occ_sims"
parent.folder <- paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR",
              HR,"mD",D*10000,sep="")
parent.GIS <- "E:/GISData/WHWO/Occ_sims/R6"

# Set working path specific to each trend scenario
path <- paste(parent.folder,"/obj1/lmb",trnd*100,sep="")
if(!dir.exists(path)) dir.create(path)

forest <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")

#### Compile ID lists and total number of transects for each forest. ####
saturated_n <- numeric(length=length(forest))
tp_match <- PID <- list()
for(f in 1:length(forest)) {
  ## Tables that matche transect with point IDs. ##
  dat <- fread(paste(parent.GIS,"/NF_",forest[f],"/","Survey_pnt_coords.txt",sep=""),header=T,sep=",")
  dat <- dat[,.(POINTID,tr=as.numeric(c(H1tr,H2tr,V1tr,V2tr)),
                          pt=as.numeric(c(H1pt,H2pt,V1pt,V2pt)))]
  tp_match[[f]] <- dat[!is.na(dat$tr),]

  ## Total number of transects identified by forest ##
  saturated_n[f] <- length(unique(tp_match[[f]]$tr))
  
  ## Point IDs in sampled order ##
  dat <- fread(paste(parent.folder,"/",forest[f],"/Lmbd",trnd*100,"/",forest[f],
                     "/rSPACEx1.txt",sep=""),header=F,sep='*')
  PID[[f]] <- as.integer(dat$V2)
}

############ Compile true values (abundance, occupancy, detectability) #################################
#***Note: use file.append() to combine Ppres chunks for FROR as necessary before running this.***
cols <- c("occtru","p.tru.med","p.tru.lo","p.tru.hi")
truth_10pts <- truth_3pts <-
  array(NA,c(nsim,20,length(cols)),dimnames=list(NULL,NULL,cols))
for(sim in 1:nsim) {
  ep_10pts <- ep_3pts <- matrix(NA,nrow=sum(saturated_n),ncol=20)
  for(f in 1:length(forest)) {
    frst <- forest[f]
    sets <- fread(paste(parent.GIS,"/NF_",frst,"/","Transect_sets.txt",sep=""),header=T,sep=",")
    trns <- sets$transects

    #Read and organize encounter probabilities for each point (rows) x year (cols) occassion
      #Note: use file.append() to combine Ppres chunks as necessary before running this.
    path <- paste(parent.folder,"/",frst,"/Lmbd",trnd*100,sep="")
    ppres <- fread(paste(path,"/Ppres.txt",sep=""),skip=(sim-1),nrows=1,header=F,sep=' ')
    ppres <- as.numeric(ppres)
    ppres <- ppres[which(!is.na(ppres))] #Remove NAs (May be one at the end)
    ppres <- round(ppres,digits=2) #Round encounter probabilities so that occupied when prob >= 0.01
    ppres <- matrix(ppres,nrow=(length(ppres)/20),ncol=20) #Convert to matrix with rows = points, cols = years
    dimnames(ppres)[[1]] <- PID[[f]] # Attach point IDs to point encounter probability matrix
    
    #Calculate encounter probabilities for each transect (rows) in a given year (cols)
    tpres_10pts <- tpres_3pts <- matrix(NA,nrow=length(trns),ncol=20,dimnames=list(trns,NULL))
    tp <- tp_match[[f]]
    ppres <- ppres[as.character(tp$POINTID),]
    ppres <- suppressWarnings(as.tbl(data.frame(ppres)))
    ppres$tr <- tp$tr
    ppres$pt <- tp$pt
    ppres_comp <- ppres %>% select(X1:X20)
    ppres_comp <- as.tbl(1 - ppres_comp)
    ppres_comp <- select(ppres,tr,pt) %>% bind_cols(ppres_comp)
    
    tpres_10pts <- 1 - (ppres_comp %>% select(tr,X1:X20) %>% group_by(tr) %>%
      summarise_each(funs(prod)) %>% select(X1:X20))
    tpres_3pts <- 1 - (ppres_comp %>% filter(pt %in% 1:3) %>% select(tr,X1:X20) %>% group_by(tr) %>%
      summarise_each(funs(prod)) %>% select(X1:X20))
    
    ifelse(f==1, st <- 1, st <- (sum(saturated_n[1:(f-1)])+1))
    end <- sum(saturated_n[1:f])
    ep_10pts[st:end,] <- as.matrix(tpres_10pts)
    ep_3pts[st:end,] <- as.matrix(tpres_3pts)
    }
  count <- (ep_10pts>=0.05)*1
  truth_10pts[sim,,"occtru"] <- apply(count,2,sum)/nrow(count)
  truth_10pts[sim,,"p.tru.med"] <- apply(ep_10pts,2,function(x) median(x[which(x>0)]))
  truth_10pts[sim,,"p.tru.lo"] <- apply(ep_10pts,2,function(x) quantile(x[which(x>0)],prob=0.025,type=8))
  truth_10pts[sim,,"p.tru.hi"] <- apply(ep_10pts,2,function(x) quantile(x[which(x>0)],prob=0.975,type=8))
  count <- (ep_3pts>=0.05)*1
  truth_3pts[sim,,"occtru"] <- apply(count,2,sum)/nrow(count)
  truth_3pts[sim,,"p.tru.med"] <- apply(ep_3pts,2,function(x) median(x[which(x>0)]))
  truth_3pts[sim,,"p.tru.lo"] <- apply(ep_3pts,2,function(x) quantile(x[which(x>0)],prob=0.025,type=8))
  truth_3pts[sim,,"p.tru.hi"] <- apply(ep_3pts,2,function(x) quantile(x[which(x>0)],prob=0.975,type=8))
  }

#Save true occupancy and true detectability values
path <- paste(parent.folder,"/obj1/lmb",trnd*100,sep="")
saveObject(truth_10pts,paste(path,"/True_psi_10PntsPerTrns",sep=""))
saveObject(truth_3pts,paste(path,"/True_psi_3PntsPerTrns",sep=""))

#Save true occupancy and true detectability values (old; delete if above works)
#path <- paste(parent.folder,"/lmb",trnd*100,"/True_psi_10PntsPerTrns/",sep="")
#if(!dir.exists(path)) dir.create(path)
#saveObject(truth_10pts,paste(path,"sim",sim,sep=""))
#path <- paste(parent.folder,"/lmb",trnd*100,"/True_psi_3PntsPerTrns/",sep="")
#if(!dir.exists(path)) dir.create(path)
#saveObject(truth_3pts,paste(path,"sim",sim,sep=""))
#######################################################################################################