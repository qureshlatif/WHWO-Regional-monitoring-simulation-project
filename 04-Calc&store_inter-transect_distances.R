## Calculates and stores all inter-transect distances ##

library(R.utils)

path.files <- "E:/GISData/WHWO/Occ_sims/R6/"
forests <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")

#Distance calculation functions
dist <- function(from,to) sqrt(((as.numeric(from[1])-as.numeric(to[,1]))^2) +
                                 ((as.numeric(from[2])-as.numeric(to[,2]))^2)) #From must be a vector of length 2 (X,Y), to can be a matrix with 2 columns (X,Y) and any number of rows
dist.min <- function(tr1,tr2) { #tr1 and tr2 must be matrices with two columns (X,Y)
  d <- apply(tr1,1,function(x) min(dist(x,tr2)))
  return(min(d))
}

for(i in 1:length(forests)) {
  dat <- read.table(paste(path.files,"NF_",forests[i],"/Survey_pnt_coords.txt",sep=""),header=T,sep=",",stringsAsFactors=F)
  dat$V2pt <- dat$V2tr <- dat$V1pt <- dat$V1tr <- dat$H2pt <- dat$H2tr <- dat$H1pt <- dat$H1tr <- ""
  for(j in 1:nrow(dat)) { # Could use sapply here instead to speed this up a bit if needed? Would need to work out how.
    if(dat$H1ID[j]!="") dat[j,c("H1tr","H1pt")] <- strsplit(dat$H1ID[j],"-")[[1]]
    if(dat$H2ID[j]!="") dat[j,c("H2tr","H2pt")] <- strsplit(dat$H2ID[j],"-")[[1]]
    if(dat$V1ID[j]!="") dat[j,c("V1tr","V1pt")] <- strsplit(dat$V1ID[j],"-")[[1]]
    if(dat$V2ID[j]!="") dat[j,c("V2tr","V2pt")] <- strsplit(dat$V2ID[j],"-")[[1]]
  }
  transects <- c(dat$H1tr[which(dat$H1tr!="")],dat$H2tr[which(dat$H2tr!="")],dat$V1tr[which(dat$V1tr!="")],
                 dat$V2tr[which(dat$V2tr!="")])
  transects <- unique(transects)
  trans.coords <- array(NA,dim=c(10,2,length(transects)))
  for(j in 1:length(transects)) {
    coords <- dat[which(dat$H1tr==transects[j]|dat$H2tr==transects[j]|
                          dat$V1tr==transects[j]|dat$V2tr==transects[j]),c("X","Y")]
    trans.coords[,,j] <- as.matrix(coords)
  }
  dist.mat <- matrix(NA,nrow=length(transects),ncol=length(transects))
  for(j in 1:length(transects)) dist.mat[j,] <- apply(trans.coords,3,function(x) dist.min(x,trans.coords[,,j]))
  saveObject(dist.mat,paste(path.files,"NF_",forests[i],"/dist_mat",sep=""))
}
