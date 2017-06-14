####################################################################################################
# Identifies sets of points capable of forming a transect and assigns a transect ID to each set.   #
# A sufficient number of transects were identified to include each point in at least one transect. #
####################################################################################################

forests <- c("COWA","DEOR","FROR","GPWA","MAOR","MHOR","OCOR","OKWA","UMOR","WAOR")

# Assign all points to one or more transects #
for(i in 1:length(forests)) {
  dat <- read.table(paste(path.files,"NF_",forests[i],"/Survey_pnt_coords.txt",sep=""),header=T,sep=",",stringsAsFactors=F)
  mark.dlt.x <- mark.dlt.y <- numeric(length=nrow(dat))
  for(j in 1:length(mark.dlt.x)) {
    xy <- dat[j,c("X","Y")]
    chck.x <- dat[-j,c("X","Y")]
    chck.x <- chck.x[which(chck.x$Y==xy$Y),]
    if(nrow(chck.x)>0) chck.x <- chck.x[which(chck.x$X<=(xy$X+2700)&chck.x$X>=(xy$X-2700)),]
    ifelse(nrow(chck.x)<9,keep.x <- FALSE,keep.x <- TRUE)
    if(keep.x==TRUE) {
      all.x <- seq(xy$X-2700,xy$X+2700,by=300)
      x.present <- is.element(all.x,c(chck.x$X,xy$X))
      keep.x <- any(sum(x.present[1:10])==10,sum(x.present[2:11])==10,sum(x.present[3:12])==10,sum(x.present[4:13])==10,
                    sum(x.present[5:14])==10,sum(x.present[6:15])==10,sum(x.present[7:16])==10,sum(x.present[8:17])==10,
                    sum(x.present[9:18])==10,sum(x.present[10:19])==10)
    }
    if(keep.x==F) mark.dlt.x[j] <- 1
    chck.y <- dat[-j,c("X","Y")]
    chck.y <- chck.y[which(chck.y$X==xy$X),]
    if(nrow(chck.y)>0) chck.y <- chck.y[which(chck.y$Y<=(xy$Y+2700)&chck.y$Y>=(xy$Y-2700)),]
    ifelse(nrow(chck.y)<9,keep.y <- FALSE,keep.y <- TRUE)
    if(keep.y==TRUE) {
      all.y <- seq(xy$Y-2700,xy$Y+2700,by=300)
      y.present <- is.element(all.y,c(chck.y$Y,xy$Y))
      keep.y <- any(sum(y.present[1:10])==10,sum(y.present[2:11])==10,sum(y.present[3:12])==10,sum(y.present[4:13])==10,
                    sum(y.present[5:14])==10,sum(y.present[6:15])==10,sum(y.present[7:16])==10,sum(y.present[8:17])==10,
                    sum(y.present[9:18])==10,sum(y.present[10:19])==10)
    }
    if(keep.y==F) mark.dlt.y[j] <- 1
  }
  dat.cols <- dat[-which(mark.dlt.y==1),] # Points that can form vertical transects
  cols <- sort(unique(dat.cols$X)) # Possible positions for vertical transects
  dat.cols$V2ID <- dat.cols$V1ID <- "" #ID columns for vertical transects
  tr <- 1
  for(j in 1:length(cols)) {
    ind <- which(dat.cols[,"X"]==cols[j]) 
    position <- dat.cols[ind,"Y"]
    position <- position - min(position)
    run.length <- numeric(length=length(position))
    for(k in 1:length(run.length)) run.length[k] <-
      length(which(is.element(position,seq(position[k],position[k]+2700,by=300))))
    possible.starts <- position[which(run.length==10)]
    st <- min(possible.starts)
    while(!is.na(st)) {
      dat.cols[ind[which(is.element(position,seq(st,st+2700,by=300)))],"V1ID"] <- paste(tr,"-",1:10,sep="")
      ifelse(any(possible.starts>=st+3000),
             st <- min(possible.starts[which(possible.starts>=st+3000)]),
             st <- NA)
      tr <- tr + 1
    }
    if(any(dat.cols[ind,"V1ID"]=="")) {
      run.length <- numeric(length=length(position))
      for(k in 1:length(run.length)) run.length[k] <-
        length(which(is.element(position,seq(position[k],position[k]-2700,by=-300))))
      possible.starts <- position[which(run.length==10)]
      st <- max(possible.starts)
      while(!is.na(st)) {
        dat.cols[ind[which(is.element(position,seq(st,st-2700,by=-300)))],"V2ID"] <- paste(tr,"-",1:10,sep="")
        ifelse(any(possible.starts<=st-3000),
             st <- max(possible.starts[which(possible.starts<=st-3000)]),
             st <- NA)
        tr <- tr + 1
      }
    }
  }
  dat.rows <- dat[-which(mark.dlt.x==1),]
  rows <- sort(unique(dat.rows$Y)) # Possible positions for vertical transects
  dat.rows$H2ID <- dat.rows$H1ID <- "" #ID columns for vertical transects
  tr <- tr + 1
  for(j in 1:length(rows)) {
    ind <- which(dat.rows[,"Y"]==rows[j]) 
    position <- dat.rows[ind,"X"]
    position <- position - min(position)
    run.length <- numeric(length=length(position))
    for(k in 1:length(run.length)) run.length[k] <-
      length(which(is.element(position,seq(position[k],position[k]+2700,by=300))))
    possible.starts <- position[which(run.length==10)]
    st <- min(possible.starts)
    while(!is.na(st)) {
      dat.rows[ind[which(is.element(position,seq(st,st+2700,by=300)))],"H1ID"] <- paste(tr,"-",1:10,sep="")
      ifelse(any(possible.starts>=st+3000),
             st <- min(possible.starts[which(possible.starts>=st+3000)]),
             st <- NA)
      tr <- tr + 1
    }
    if(any(dat.rows[ind,"H1ID"]=="")) {
      run.length <- numeric(length=length(position))
      for(k in 1:length(run.length)) run.length[k] <-
        length(which(is.element(position,seq(position[k],position[k]-2700,by=-300))))
      possible.starts <- position[which(run.length==10)]
      st <- max(possible.starts)
      while(!is.na(st)) {
        dat.rows[ind[which(is.element(position,seq(st,st-2700,by=-300)))],"H2ID"] <- paste(tr,"-",1:10,sep="")
        ifelse(any(possible.starts<=st-3000),
             st <- max(possible.starts[which(possible.starts<=st-3000)]),
             st <- NA)
        tr <- tr + 1
      }
    }
  }
  dat$V2ID <- dat$V1ID <- dat$H2ID <- dat$H1ID <- ""
  for(j in 1:nrow(dat)) {
    if(is.element(dat$POINTID[j],dat.cols$POINTID)) {
      dat$V1ID[j] <- dat.cols$V1ID[which(dat.cols$POINTID==dat$POINTID[j])]
      dat$V2ID[j] <- dat.cols$V2ID[which(dat.cols$POINTID==dat$POINTID[j])]
    }
    if(is.element(dat$POINTID[j],dat.rows$POINTID)) {
      dat$H1ID[j] <- dat.rows$H1ID[which(dat.rows$POINTID==dat$POINTID[j])]
      dat$H2ID[j] <- dat.rows$H2ID[which(dat.rows$POINTID==dat$POINTID[j])]
    }
  }
  write.table(dat,paste(path.files,"NF_",forests[i],"/Survey_pnt_coords.txt",sep=""),row.names=F,sep=",")
}
