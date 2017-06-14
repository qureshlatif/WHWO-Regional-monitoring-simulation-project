####################################################################
# Identify sets of transects with >= 2km spacing between neighbors #
####################################################################

library(R.utils)

path.files <- "E:/GISData/WHWO/Occ_sims/R6/"
forests <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")

# Assign transects to sets whose members are adequately spaced #
for(i in 1:length(forests)) {
  dat <- read.table(paste(path.files,"NF_",forests[i],"/Survey_pnt_coords.txt",sep=""),header=T,sep=",",stringsAsFactors=F)
  dat$V2pt <- dat$V2tr <- dat$V1pt <- dat$V1tr <- dat$H2pt <- dat$H2tr <- dat$H1pt <- dat$H1tr <- ""
  for(j in 1:nrow(dat)) {
  # Could use sapply here instead to speed this up a bit if needed. Would need to work out how.
    if(dat$H1ID[j]!="") dat[j,c("H1tr","H1pt")] <- strsplit(dat$H1ID[j],"-")[[1]]
    if(dat$H2ID[j]!="") dat[j,c("H2tr","H2pt")] <- strsplit(dat$H2ID[j],"-")[[1]]
    if(dat$V1ID[j]!="") dat[j,c("V1tr","V1pt")] <- strsplit(dat$V1ID[j],"-")[[1]]
    if(dat$V2ID[j]!="") dat[j,c("V2tr","V2pt")] <- strsplit(dat$V2ID[j],"-")[[1]]
  }
  transects <- c(dat$H1tr[which(dat$H1tr!="")],dat$H2tr[which(dat$H2tr!="")],dat$V1tr[which(dat$V1tr!="")],
                 dat$V2tr[which(dat$V2tr!="")])
  transects <- unique(transects)
  dist.mat <- loadObject(paste(path.files,"NF_",forests[i],"/dist_mat",sep=""))
  S <- numeric(length=length(transects))
  S[1] <- 1
  for(j in 2:length(transects)) if(!any(dist.mat[j,which(S==1)]<=2000)) S[j] <- 1
  out <- data.frame(cbind(transects,S1=S),stringsAsFactors=F)
  out$S1 <- as.numeric(out$S1)
  counter <- 2
  while(ifelse(ncol(out)==2,any(out[,2]==0),any(apply(out[,2:ncol(out)],1,sum)==0))) {
    S <- numeric(length=length(transects))
    ifelse(ncol(out)==2,ind.left <- which(out[,2]==0),ind.left <- which(apply(out[,2:ncol(out)],1,sum)==0))
    S[ind.left[1]] <- 1
    if(length(ind.left)>1)
      for(j in ind.left[2:length(ind.left)]) if(!any(dist.mat[j,which(S==1)]<=2000)) S[j] <- 1
    out <- cbind(out,S)
    names(out)[counter+1] <- paste("S",counter,sep="")
    counter <- counter + 1
  }
  for(j in 2:ncol(out)) for(k in 1:nrow(out))
    if(!any(dist.mat[k,which(out[,j]==1)]<=2000)) out[k,j] <- 1
  write.table(dat,paste(path.files,"NF_",forests[i],"/Survey_pnt_coords.txt",sep=""),row.names=F,sep=",")
  write.table(out,paste(path.files,"NF_",forests[i],"/Transect_sets.txt",sep=""),row.names=F,sep=",")
}
 
  
