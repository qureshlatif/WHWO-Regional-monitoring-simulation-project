library(R.utils)

#### Variables ####
HR <- 1 #home range radius in kilometers
N <- c(seq(10,34,by=4),36,seq(38,90,by=4)) #Transect-specific abundances
sims <- 30 #Set to number of simulated datasets.
###################

# Set working directory to folder with simulated datasets:
path <- paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR",HR*1000,"m_yr1",sep="")
setwd(path)

TR.names <- c("TR601","TR602","TR603","TR604","TR605","TR606","TR607","TR608","TR609","TR611","TR612",
              "TR613","TR614","TR615","TR616","TR617","TR618","TR619","TR620","TR621","TR623","TR624","TR625","TR626",
              "TR628","TR629","TR630","TR640","TR646","TR647")
dat.files <- paste("rSPACEx",seq(1:sims),".txt",sep="")
n.vsts <- 2

library(foreign)
TID <- TR.names
n.trns <- length(TID)
n.yrs <- 1
n.pnts <- 10 #No. points per transect

#Descriptives
cols <- c("psi_tr.app","psi.app","det.tot","p.app")
out <- matrix("",nrow=length(N),ncol=length(cols),dimnames=list(N,cols))

for (n in N) {
  out.desc <- array(NA,dim=c(sims,length(cols)),dimnames=list(NULL,cols))
  for(i in 1:sims) {
    # Compile data
    Y.mat <- array(NA,dim=c(length(TID),10),dimnames=list(TR.names,NULL))
    path <- getwd()
    for(tr in 1:length(TID)) {
      dat <- read.delim(paste(path,"/constN",n,"/",TR.names[tr],"/",dat.files[i],sep=""),header=F,sep='*',as.is=T)
      dat <- dat[order(dat$V2),]
      dat <- substr(dat[,3],3,(n.vsts+2))
      dat <- strsplit(dat,split='')
      for(pt in 1:10) {
        obs <- as.numeric(dat[[pt]])
        obs <- matrix(obs,nrow=n.vsts,ncol=n.yrs)
        Y.mat[tr,pt] <- apply(obs,2,sum)
      }
    }
    #Descriptives
    out.desc[i,"psi_tr.app"] <- sum(apply(Y.mat,1,sum)>0)/n.trns
    out.desc[i,"psi.app"] <- length(which(Y.mat>0))/length(Y.mat)
    out.desc[i,"det.tot"] <- sum(Y.mat)
    out.desc[i,"p.app"] <- sum(Y.mat)/(length(which(Y.mat>0))*2)
    }
  out[as.character(n),"psi_tr.app"] <- paste(round(median(out.desc[,"psi_tr.app"]),digits=3),"(",
                                             round(quantile(out.desc[,"psi_tr.app"],prob=0.025,type=8),digits=3),",",
                                             round(quantile(out.desc[,"psi_tr.app"],prob=0.975,type=8),digits=3),")",sep="")
  out[as.character(n),"psi.app"] <- paste(round(median(out.desc[,"psi.app"]),digits=3),"(",
                                             round(quantile(out.desc[,"psi.app"],prob=0.025,type=8),digits=3),",",
                                             round(quantile(out.desc[,"psi.app"],prob=0.975,type=8),digits=3),")",sep="")
  out[as.character(n),"det.tot"] <- paste(round(median(out.desc[,"det.tot"])),"(",
                                             round(quantile(out.desc[,"det.tot"],prob=0.025,type=8)),",",
                                             round(quantile(out.desc[,"det.tot"],prob=0.975,type=8)),")",sep="")
  out[as.character(n),"p.app"] <- paste(round(median(out.desc[,"p.app"]),digits=3),"(",
                                             round(quantile(out.desc[,"p.app"],prob=0.025,type=8),digits=3),",",
                                             round(quantile(out.desc[,"p.app"],prob=0.975,type=8),digits=3),")",sep="")
  }

row.names(out) <- paste("N = ",row.names(out)," (D = ",round(as.numeric(row.names(out))/150.664,digits=4)," per 70 ha)",sep="")
write.csv(out,"Summary.csv",row.names=T)
