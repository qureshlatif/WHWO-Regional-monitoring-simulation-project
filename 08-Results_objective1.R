require(dplyr)
setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6")

#############################################################
## Tabulate and plot power for detecting population trends ##
#############################################################

library(R.utils)

base.folder <- "HR1000mD2522/obj1/"
trnds <- c(0.9,0.95,0.98,1)
no.trns <- c(60,90,120,150)

for(ld in trnds) for(nt in no.trns) {
  load(paste(base.folder,"lmb",ld*100,"/Tr",nt,"Pt10Mf1Results.RData",sep=""))
  pow.cp <- (length(which(PSIprd_cp.trend[,3]<0))/dim(PSIprd_cp.trend)[1])*100
  pow.yp <- (length(which(PSIprd_yp.trend[,3]<0))/dim(PSIprd_yp.trend)[1])*100
  if(ld == 0.98 & nt %in% c(60, 90) | nt == 150 | ld == 1) {
    load(paste(base.folder,"lmb",ld*100,"/Tr",nt*2,"Pt10Mf1Results_lr.RData",sep=""))
    pow.lr <- length(which(PSIprd_lr.trend[,3]<0))
    spur.pos.lr <- (length(which(PSIprd_lr.trend[,2]>0))/dim(PSIprd_lr.trend)[1])*100
  } else {
    pow.lr <- 100
    spur.pos.lr <- 0
  }
  ifelse(ld==0.9&nt==60,
         out.power <- c(nt, ld, pow.cp, pow.yp, pow.lr),
         out.power <- rbind(out.power, c(nt, ld, pow.cp, pow.yp, pow.lr))) 
  spur.pos.cp <- (length(which(PSIprd_cp.trend[,2]>0))/dim(PSIprd_cp.trend)[1])*100
  spur.pos.yp <- (length(which(PSIprd_yp.trend[,2]>0))/dim(PSIprd_yp.trend)[1])*100
  ifelse(ld==0.9&nt==60,
         out.spurious <- c(nt, ld, spur.pos.cp, spur.pos.yp, spur.pos.lr),
         out.spurious <- rbind(out.spurious, c(nt, ld, spur.pos.cp, spur.pos.yp, spur.pos.lr))) 
}
dimnames(out.power) <- dimnames(out.spurious) <-
  list(NULL, c("no.trns", "lambda", "constant-p", "yearly-p", "logistic"))

## Add two lambda = 0.98 scenarios for log reg with lower sampling effort ##
out.power <- rbind(out.power, c(45, 0.98, NA, NA, 0), c(30, 0.98, NA, NA, 0))
load("HR1000mD2522/obj1/lmb98/Tr90Pt10Mf1Results_lr.RData")
out.power[17, "logistic"] <- length(which(PSIprd_lr.trend[,3]<0))
load("HR1000mD2522/obj1/lmb98/Tr60Pt10Mf1Results_lr.RData")
out.power[18, "logistic"] <- length(which(PSIprd_lr.trend[,3]<0))

## Tabulate spurious trend detection ##
out <- out.power[which(out.power[,"lambda"]==1),c("constant-p", "yearly-p", "logistic")] + # Add spurious negative trends...
  out.spurious[which(out.spurious[,"lambda"]==1),c("constant-p", "yearly-p", "logistic")]  # with spurious positive trends.
out <- cbind(out.power[which(out.power[,"lambda"]==1),c("no.trns","lambda")],out)
write.csv(out,"manuscript/Spurious_trend_summary.csv")

## Plot power curves ##
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

  # lambda = 0.9
dat <- data.frame(rbind(out.power[which(out.power[,"lambda"]==0.9),c("no.trns","constant-p")],
              out.power[which(out.power[,"lambda"]==0.9),c("no.trns","yearly-p")],
              out.power[which(out.power[,"lambda"]==0.9),c("no.trns","logistic")]))
names(dat)[2] <- "power"
dat$Model <- c(rep("Constant-p",4),rep("Yearly-p",4),rep("Logistic reg",4))
dat$Model <- factor(dat$Model, levels = c("Constant-p", "Yearly-p", "Logistic reg"))

dat[which(dat$Model == "Constant-p"), "power"] <- 99.5    #Offset for visibility
dat[which(dat$Model == "Logistic reg"), "power"] <- 100.5 #Offset for visibility

theme_set(theme_bw())
pA <- ggplot(data = dat,aes(x=no.trns,y=power,linetype=Model,shape=Model)) + 
  geom_line(aes(linetype=Model),size=1.5) +
  geom_point(aes(shape=Model),size=7) +
  geom_hline(yintercept=80,linetype="dotted",size=1) +
  scale_x_continuous(breaks=c(60,90,120,150),
                     labels = c("60 (120)", "90 (180)", "120 (240)", "150 (300)")) +
  scale_y_continuous(limits=c(0,101),breaks=seq(0,100,by=20)) +
  ylab(NULL) + xlab(NULL) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotdash')) +
  scale_shape_manual(values = c(16, 18, 15)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20)) +
#  theme(axis.title.x=element_text(size=30)) +
#  theme(axis.title.y=element_text(size=30)) +
  guides(linetype=FALSE) +
  guides(shape=FALSE) +
  geom_text(aes(x=130,y=5,label="lambda[N]"),size=10,parse=T) +
  geom_text(aes(x=145,y=5,label="= 0.90"),size=9)

  # lambda = 0.95
dat <- data.frame(rbind(out.power[which(out.power[,"lambda"]==0.95),c("no.trns","constant-p")],
              out.power[which(out.power[,"lambda"]==0.95),c("no.trns","yearly-p")],
              out.power[which(out.power[,"lambda"]==0.95),c("no.trns","logistic")]))
names(dat)[2] <- "power"
dat$Model <- c(rep("Constant-p",4),rep("Yearly-p",4),rep("Logistic reg",4))
dat$Model <- factor(dat$Model, levels = c("Constant-p", "Yearly-p", "Logistic reg"))

dat[which(dat$Model == "Constant-p"), "power"] <- 99.5    #Offset for visibility
dat[which(dat$Model == "Logistic reg"), "power"] <- 100.5 #Offset for visibility

theme_set(theme_bw())
pB <- ggplot(data = dat,aes(x=no.trns,y=power,linetype=Model,shape=Model)) + 
  geom_line(aes(linetype=Model),size=1.5) +
  geom_point(aes(shape=Model),size=7) +
  geom_hline(yintercept=80,linetype="dotted",size=1) +
  scale_x_continuous(breaks=c(60,90,120,150),
                     labels = c("60 (120)", "90 (180)", "120 (240)", "150 (300)")) +
  scale_y_continuous(limits=c(0,101),breaks=seq(0,100,by=20)) +
  ylab(NULL) + xlab(NULL) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotdash')) +
  scale_shape_manual(values = c(16, 18, 15)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20)) +
#  theme(axis.title.x=element_text(size=30)) +
#  theme(axis.title.y=element_text(size=30)) +
  guides(linetype=FALSE) +
  guides(shape=FALSE) +
  geom_text(aes(x=130,y=5,label="lambda[N]"),size=10,parse=T) +
  geom_text(aes(x=145,y=5,label="= 0.95"),size=9)

  # lambda = 0.98
dat <- data.frame(rbind(out.power[which(out.power[,"lambda"]==0.98),c("no.trns","constant-p")],
              out.power[which(out.power[,"lambda"]==0.98),c("no.trns","yearly-p")],
              out.power[which(out.power[,"lambda"]==0.98),c("no.trns","logistic")]))
names(dat)[2] <- "power"
dat <- dat[-which(is.na(dat$power)),]
dat$Model <- c(rep("Constant-p",4),rep("Yearly-p",4),rep("Logistic regression",6))
dat$Model <- factor(dat$Model, levels = c("Constant-p", "Yearly-p", "Logistic regression"))

theme_set(theme_bw())
pC <- ggplot(data = dat,aes(x=no.trns,y=power,linetype=Model,shape=Model)) + 
  geom_line(aes(linetype=Model),size=1.5) +
  geom_point(aes(shape=Model),size=7) +
  geom_hline(yintercept=80,linetype="dotted",size=1) +
  scale_x_continuous(breaks=c(30, 45, 60, 90, 120, 150),
                     labels = c("30 (60)", "45 (90)", "60 (120)", "90 (180)", "120 (240)", "150 (300)")) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=20)) +
  ylab(NULL) + xlab(NULL) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotdash')) +
  scale_shape_manual(values = c(16, 18, 15)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.title.y=element_text(size=30)) +
  theme(legend.text=element_text(size=25)) +
  theme(legend.title=element_text(size=30)) +
  geom_text(aes(x=130,y=5,label="lambda[N]"),size=10,parse=T) +
  geom_text(aes(x=145,y=5,label="= 0.98"),size=9)


theme_set(theme_bw())
p <- ggdraw() + 
    draw_plot(pA, x = 0.05, y = 0.55, width = 0.43, height = 0.45) +
    draw_plot(pB, x = 0.5, y = 0.55, width = 0.43, height = 0.45) +
    draw_plot(pC, x = 0.05, y = 0.1, width = 0.94, height = 0.45) +
    draw_plot_label(label = c("Number of Transects (for logistic regression)",
                              "Power to observe decline"),size=30,
                    x=c(0,0), y=c(0.1,0.2), angle=c(0, 90))
#p
save_plot("manuscript/Figure_Obj1.tiff", p, ncol = 3.5, nrow = 3.5)


### Reduced figure for presentation ###

# lambda = 0.98
dat <- data.frame(rbind(out.power[which(out.power[,"lambda"]==0.98),c("no.trns","constant-p")],
                        out.power[which(out.power[,"lambda"]==0.98),c("no.trns","logistic")]))
names(dat)[2] <- "power"
dat <- dat[-which(is.na(dat$power)),]
dat$Model <- c(rep("Constant-p",4),rep("Logistic regression",6))
dat$Model <- factor(dat$Model, levels = c("Constant-p", "Logistic regression"))

theme_set(theme_bw())
p <- ggplot(data = dat,aes(x=no.trns,y=power,linetype=Model,shape=Model)) + 
  geom_line(aes(linetype=Model),size=1.5) +
  geom_point(aes(shape=Model),size=7) +
  geom_hline(yintercept=80,linetype="dotted",size=1) +
  scale_x_continuous(breaks=c(30, 45, 60, 90, 120, 150),
                     labels = c("60", "90", "120", "180", "240", "300")) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,by=20)) +
  ylab("Power to observe decline") + xlab("Number of surveys per year") +
  scale_linetype_manual(values = c('solid', 'dotdash')) +
  scale_shape_manual(values = c(16, 15)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.title.y=element_text(size=30)) +
  theme(legend.text=element_text(size=25)) +
  theme(legend.title=element_text(size=30)) +
  geom_text(aes(x=130,y=5,label="lambda[N]"),size=10,parse=T) +
  geom_text(aes(x=145,y=5,label="= 0.98"),size=9)

save_plot("manuscript/Figure_Power_MTTWS_2017.tiff", p, ncol = 3, nrow = 2)


#########################################################
## Plot occupancy estimates and true occupancy by year ##
#########################################################
theme_set(theme_bw())
n.yrs <- 20
SIMS <- 30

### 10 points per transect ###
## Lambda = 0.9 ##
  # constant-p model
trnd <- 0.9 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")
dat.true <- data.frame(cbind(year = seq(1,20), PSI = apply(trvls[,,"occtru"],2,mean)))

PSIprd.sim <- PSIprd_cp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="J",size=9)
p.90cp <- p

  # yearly-p model
PSIprd.sim <- PSIprd_yp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="K",size=9)
p.90yp <- p

  # logistic regression
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model
sims <- SIMS

PSIprd.sim <- PSIprd_lr.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="L",size=9)
p.90lr <- p

## Lambda = 0.95 ##
  # constant-p model
trnd <- 0.95 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")
dat.true <- data.frame(cbind(year = seq(1,20), PSI = apply(trvls[,,"occtru"],2,mean)))

PSIprd.sim <- PSIprd_cp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="G",size=9)
p.95cp <- p

  # yearly-p model
PSIprd.sim <- PSIprd_yp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="H",size=9)
p.95yp <- p

  # logistic regression
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model
sims <- SIMS

PSIprd.sim <- PSIprd_lr.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="I",size=9)
p.95lr <- p

## Lambda = 0.98 ##
  # constant-p model
trnd <- 0.98 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")
dat.true <- data.frame(cbind(year = seq(1,20), PSI = apply(trvls[,,"occtru"],2,mean)))

PSIprd.sim <- PSIprd_cp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="D",size=9)
p.98cp <- p

  # yearly-p model
PSIprd.sim <- PSIprd_yp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="E",size=9)
p.98yp <- p

  # logistic regression
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model
sims <- SIMS

PSIprd.sim <- PSIprd_lr.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="F",size=9)
p.98lr <- p

## Lambda = 1 ##
  # constant-p model
trnd <- 1 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")
dat.true <- data.frame(cbind(year = seq(1,20), PSI = apply(trvls[,,"occtru"],2,mean)))

PSIprd.sim <- PSIprd_cp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="A",size=9)
p.100cp <- p

  # yearly-p model
PSIprd.sim <- PSIprd_yp.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="B",size=9)
p.100yp <- p

  # logistic regression
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model
sims <- SIMS

PSIprd.sim <- PSIprd_lr.sim[1:sims,,] # Adjust to match desired model (constant-p or yearly-p)

cols <- c("sim","Year","PSI","PSI.lo","PSI.hi","PSI.true")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$PSI <- as.numeric(PSIprd.sim[,,1]) 
dat$PSI.lo <- as.numeric(PSIprd.sim[,,2]) 
dat$PSI.hi <- as.numeric(PSIprd.sim[,,3]) 
p <- ggplot(data = dat,aes(x=Year,y=PSI))
for(i in 1:sims) {
  p <- p + geom_line(data=dat[which(dat$sim==i),],aes(x=Year,y=PSI),alpha=0.3,size=1)
}
p <- p + geom_errorbar(aes(ymin=PSI.lo,ymax=PSI.hi),width=0.1,alpha=0.3,color="blue")
p <- p + geom_point(size=3,alpha=0.4,color="black")
p <- p + geom_point(dat = dat.true, aes(x=year, y=PSI), size=4, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=20,y=1.05),label="C",size=9)
p.100lr <- p

theme_set(theme_bw())
p <- ggdraw() + 
    draw_plot(p.100cp, x = 0.05, y = 0.7625, width = 0.316, height = 0.2375) +
    draw_plot(p.100yp, x = 0.366, y = 0.7625, width = 0.316, height = 0.2375) +
    draw_plot(p.100lr, x = 0.683, y = 0.7625, width = 0.316, height = 0.2375) +
    draw_plot(p.98cp, x = 0.05, y = 0.525, width = 0.316, height = 0.2375) +
    draw_plot(p.98yp, x = 0.366, y = 0.525, width = 0.316, height = 0.2375) +
    draw_plot(p.98lr, x = 0.683, y = 0.525, width = 0.316, height = 0.2375) +
    draw_plot(p.95cp, x = 0.05, y = 0.2875, width = 0.316, height = 0.2375) +
    draw_plot(p.95yp, x = 0.366, y = 0.2875, width = 0.316, height = 0.2375) +
    draw_plot(p.95lr, x = 0.683, y = 0.2875, width = 0.316, height = 0.2375) +
    draw_plot(p.90cp, x = 0.05, y = 0.05, width = 0.316, height = 0.2375) +
    draw_plot(p.90yp, x = 0.366, y = 0.05, width = 0.316, height = 0.2375) +
    draw_plot(p.90lr, x = 0.683, y = 0.05, width = 0.316, height = 0.2375) +
    draw_plot_label(label = c("Year","hat(psi)"), size=c(40,45),
                    x=c(0.45,0), y=c(0.05,0.6), parse = T)
#p

setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/")
save_plot("manuscript/Figure_Psi_estimates.tiff", p, ncol = 4, nrow = 4, dpi = 50)

###########################
# Detectability estimates #
###########################
theme_set(theme_bw())

## constant-p models ##
n.yrs <- 20
SIMS <- 30

trnd <- 1 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","p","p.lo","p.hi","lambda")
dat.plot <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat.plot$sim <- seq(1,sims)
dat.plot$lambda <- trnd
dat.plot$p <- as.numeric(pprd_cp.sim[1:sims,1])
dat.plot$p.lo <- as.numeric(pprd_cp.sim[1:sims,2])
dat.plot$p.hi <- as.numeric(pprd_cp.sim[1:sims,3])

dat.true <- data.frame(cbind(lmb = rep(trnd, dim(trvls)[2]), p = as.numeric(trvls[,,"p.tru.med"])))

trnd <- 0.98 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","p","p.lo","p.hi","lambda")
add.plot <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
add.plot$sim <- seq(1,sims)
add.plot$lambda <- trnd
add.plot$p <- as.numeric(pprd_cp.sim[1:sims,1])
add.plot$p.lo <- as.numeric(pprd_cp.sim[1:sims,2])
add.plot$p.hi <- as.numeric(pprd_cp.sim[1:sims,3])
dat.plot <- rbind(dat.plot,add.plot)
rm(add.plot)

add.true <- data.frame(cbind(lmb = rep(trnd, dim(trvls)[2]), p = as.numeric(trvls[,,"p.tru.med"])))
dat.true <- rbind(dat.true, add.true)
rm(add.true)

trnd <- 0.95 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","p","p.lo","p.hi","lambda")
add.plot <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
add.plot$sim <- seq(1,sims)
add.plot$lambda <- trnd
add.plot$p <- as.numeric(pprd_cp.sim[1:sims,1])
add.plot$p.lo <- as.numeric(pprd_cp.sim[1:sims,2])
add.plot$p.hi <- as.numeric(pprd_cp.sim[1:sims,3])
dat.plot <- rbind(dat.plot,add.plot)
rm(add.plot)

add.true <- data.frame(cbind(lmb = rep(trnd, dim(trvls)[2]), p = as.numeric(trvls[,,"p.tru.med"])))
dat.true <- rbind(dat.true, add.true)
rm(add.true)

trnd <- 0.9 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","p","p.lo","p.hi","lambda")
add.plot <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
add.plot$sim <- seq(1,sims)
add.plot$lambda <- trnd
add.plot$p <- as.numeric(pprd_cp.sim[1:sims,1])
add.plot$p.lo <- as.numeric(pprd_cp.sim[1:sims,2])
add.plot$p.hi <- as.numeric(pprd_cp.sim[1:sims,3])
dat.plot <- rbind(dat.plot,add.plot)
rm(add.plot)

add.true <- data.frame(cbind(lmb = rep(trnd, dim(trvls)[2]), p = as.numeric(trvls[,,"p.tru.med"])))
dat.true <- rbind(dat.true, add.true)
rm(add.true)

dat.plot$lambda <- dat.plot$lambda + runif(nrow(dat.plot),-0.008,0.008)
dat.true$lmb <- dat.true$lmb + runif(nrow(dat.true),-0.008,0.008)

p <- ggplot(data = dat.plot,aes(x=lambda,y=p))
p <- p + geom_errorbar(aes(ymin=p.lo,ymax=p.hi),width=0.001,alpha=0.15,color="blue")
p <- p + geom_point(size=3,alpha=0.3,color="black")
p <- p + geom_point(data = dat.true, aes(x=lmb, y=p), size=3, alpha=0.3, color="red")
p <- p + ylab(NULL) + xlab(expression(lambda[N])) +
  scale_x_continuous(breaks=c(0.9, 0.95, 0.98, 1)) +
  scale_y_continuous(lim = c(0, 1.05), breaks = seq(0, 1, by = 0.2)) +
  theme(axis.text.x=element_text(size=30)) +
  theme(axis.title.x=element_text(size=45)) +
  theme(axis.text.y=element_text(size=30)) +
  geom_text(aes(x=0.88,y=1.05),label="A",size=12)
#p
p.cp <- p

## yearly-p model ##
trnd <- 1 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","Year","p","p.lo","p.hi")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$p <- as.numeric(pprd_yp.sim[1:sims,,1])
dat$p.lo <- as.numeric(pprd_yp.sim[1:sims,,2])
dat$p.hi <- as.numeric(pprd_yp.sim[1:sims,,3])

dat.true <- data.frame(cbind(year = rep(seq(1,20), each = 30), p = as.numeric(trvls[,,"p.tru.med"])))
dat.true$year <- dat.true$year + runif(nrow(dat.true), -0.2, 0.2)

p <- ggplot(data = dat,aes(x=Year,y=p))
p <- p + geom_errorbar(aes(ymin=p.lo,ymax=p.hi),width=0.1,alpha=0.15,color="blue")
p <- p + geom_point(size=3,alpha=0.3,color="black")
p <- p + geom_point(data = dat.true, aes(x=year, y=p), size=3, alpha=0.3, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks=c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +  
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=1,y=1.05),label="B",size=12)
p.l100yp <- p

## yearly-p model ##
trnd <- 0.98 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","Year","p","p.lo","p.hi")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$p <- as.numeric(pprd_yp.sim[1:sims,,1])
dat$p.lo <- as.numeric(pprd_yp.sim[1:sims,,2])
dat$p.hi <- as.numeric(pprd_yp.sim[1:sims,,3])

dat.true <- data.frame(cbind(year = rep(seq(1,20), each = 30), p = as.numeric(trvls[,,"p.tru.med"])))
dat.true$year <- dat.true$year + runif(nrow(dat.true), -0.2, 0.2)

p <- ggplot(data = dat,aes(x=Year,y=p))
p <- p + geom_errorbar(aes(ymin=p.lo,ymax=p.hi),width=0.1,alpha=0.15,color="blue")
p <- p + geom_point(size=3,alpha=0.3,color="black")
p <- p + geom_point(data = dat.true, aes(x=year, y=p), size=3, alpha=0.3, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks=c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +  
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=1,y=1.05),label="C",size=12)
p.l98yp <- p

## yearly-p model ##
trnd <- 0.95 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","Year","p","p.lo","p.hi")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$p <- as.numeric(pprd_yp.sim[1:sims,,1])
dat$p.lo <- as.numeric(pprd_yp.sim[1:sims,,2])
dat$p.hi <- as.numeric(pprd_yp.sim[1:sims,,3])

dat.true <- data.frame(cbind(year = rep(seq(1,20), each = 30), p = as.numeric(trvls[,,"p.tru.med"])))
dat.true$year <- dat.true$year + runif(nrow(dat.true), -0.2, 0.2)

p <- ggplot(data = dat,aes(x=Year,y=p))
p <- p + geom_errorbar(aes(ymin=p.lo,ymax=p.hi),width=0.1,alpha=0.15,color="blue")
p <- p + geom_point(size=3,alpha=0.3,color="black")
p <- p + geom_point(data = dat.true, aes(x=year, y=p), size=3, alpha=0.3, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks=c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +  
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=1,y=1.05),label="D",size=12)
p.l95yp <- p

## yearly-p model ##
trnd <- 0.9 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model
sims <- SIMS

trvls <- loadObject("True_psi_10PntsPerTrns")

cols <- c("sim","Year","p","p.lo","p.hi")
dat <- data.frame(matrix(NA,nrow=n.yrs*sims,ncol=length(cols),dimnames=list(NULL,cols)))
dat$sim <- rep(seq(1,sims),n.yrs)
dat$Year <- rep(seq(1,n.yrs),each=sims)
dat$Year <- dat$Year + runif(sims,-0.2,0.2)
dat$p <- as.numeric(pprd_yp.sim[1:sims,,1])
dat$p.lo <- as.numeric(pprd_yp.sim[1:sims,,2])
dat$p.hi <- as.numeric(pprd_yp.sim[1:sims,,3])

dat.true <- data.frame(cbind(year = rep(seq(1,20), each = 30), p = as.numeric(trvls[,,"p.tru.med"])))
dat.true$year <- dat.true$year + runif(nrow(dat.true), -0.2, 0.2)

p <- ggplot(data = dat,aes(x=Year,y=p))
p <- p + geom_errorbar(aes(ymin=p.lo,ymax=p.hi),width=0.1,alpha=0.15,color="blue")
p <- p + geom_point(size=3,alpha=0.3,color="black")
p <- p + geom_point(data = dat.true, aes(x=year, y=p), size=3, alpha=0.3, color="red")
p <- p + ylab(NULL) + xlab(NULL) +
  scale_x_continuous(breaks=c(1, 5, 10, 15, 20)) +  
  scale_y_continuous(lim = c(0,1.05), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +  
  theme(axis.text.x=element_text(size=17)) +
  theme(axis.text.y=element_text(size=20)) +
  geom_text(aes(x=1,y=1.05),label="E",size=12)
p.l90yp <- p

p <- ggdraw() + 
    draw_plot(p.cp, x = 0.05, y = 0.682, width = 0.95, height = 0.316) +
    draw_plot(p.l100yp, x = 0.05, y = 0.366, width = 0.475, height = 0.316) +
    draw_plot(p.l98yp, x = 0.525, y = 0.366, width = 0.475, height = 0.316) +
    draw_plot(p.l95yp, x = 0.05, y = 0.05, width = 0.475, height = 0.316) +
    draw_plot(p.l90yp, x = 0.525, y = 0.05, width = 0.475, height = 0.316) +
    draw_plot_label(label = c("Year", "hat(p)"), size=c(40, 45),
                    x=c(0.45,0), y=c(0.05,0.6), parse = T)
#p

setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/")
save_plot("manuscript/Figure_p_estimates.tiff", p, ncol = 3, nrow = 4, dpi = 50)

####################################################
## Plot occupancy trend estimates vs true lambdas ##
####################################################
library(ggplot2)
library(cowplot)
n.yrs <- 20
SIMS <- 30

  # Lambda = 0.9
    # Occupancy models
trnd <- 0.9 #Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model

dat.cp <- data.frame(exp(PSIprd_cp.trend[,1:3]))
dat.yp <- data.frame(exp(PSIprd_yp.trend[,1:3]))
names(dat.cp) <- names(dat.yp) <- c("med","lo","hi")
dat.cp$lmbN <- dat.yp$lmbN <- rep(trnd,nrow(dat.cp))

RMSE.cp <- paste(".",
  substr(as.character(round(sqrt(mean(PSIprd_cp.trend[, "MSE"])), digits = 3)),3,5),sep="")
RMSE.yp <- paste(".",
  substr(as.character(round(sqrt(mean(PSIprd_yp.trend[, "MSE"])), digits = 3)),3,5),sep="")

RMSE.cp.psi <- paste(".",
  substr(as.character(round(sqrt(mean(PSIprd_cp.trend[, "MSE.psi"])), digits = 3)),3,5),sep="")
RMSE.yp.psi <- paste(".",
  substr(as.character(round(sqrt(mean(PSIprd_yp.trend[, "MSE.psi"])), digits = 3)),3,5),sep="")

    #Logistic regression
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model

dat.lr <- data.frame(exp(PSIprd_lr.trend[,1:3]))
names(dat.lr) <- c("med","lo","hi")
dat.lr$lmbN <- rep(trnd,nrow(dat.cp))

RMSE.lr <- paste(".",
  substr(as.character(round(sqrt(mean(PSIprd_lr.trend[, "MSE"])), digits = 3)),3,5),sep="")
RMSE.lr.psi <- paste(".",
  substr(as.character(round(sqrt(mean(PSIprd_lr.trend[, "MSE.psi"])), digits = 3)),3,5),sep="")

  # Lambda = 0.95
    # Occupancy models
trnd <- 0.95 # Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model

dat.cp <- dat.cp %>% rbind(exp(PSIprd_cp.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_cp.trend))) %>%
                             data.frame() %>% setNames(names(dat.cp)))
dat.yp <- dat.yp %>% rbind(exp(PSIprd_yp.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_yp.trend))) %>%
                             data.frame() %>% setNames(names(dat.yp)))

RMSE.cp <- RMSE.cp %>% c("." %>% paste(mean(PSIprd_cp.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))
RMSE.yp <- RMSE.yp %>% c("." %>% paste(mean(PSIprd_yp.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

RMSE.cp.psi <- RMSE.cp.psi %>% c("." %>% paste(mean(PSIprd_cp.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))
RMSE.yp.psi <- RMSE.yp.psi %>% c("." %>% paste(mean(PSIprd_yp.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

    # Logistic regression
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model

dat.lr <- dat.lr %>% rbind(exp(PSIprd_lr.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_lr.trend))) %>%
                             data.frame() %>% setNames(names(dat.lr)))

RMSE.lr <- RMSE.lr %>% c("." %>% paste(mean(PSIprd_lr.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

RMSE.lr.psi <- RMSE.lr.psi %>% c("." %>% paste(mean(PSIprd_lr.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

  # Lambda = 0.98
    #Occupancy
trnd <- 0.98 # Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model

dat.cp <- dat.cp %>% rbind(exp(PSIprd_cp.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_cp.trend))) %>%
                             data.frame() %>% setNames(names(dat.cp)))
dat.yp <- dat.yp %>% rbind(exp(PSIprd_yp.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_yp.trend))) %>%
                             data.frame() %>% setNames(names(dat.yp)))

RMSE.cp <- RMSE.cp %>% c("." %>% paste(mean(PSIprd_cp.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))
RMSE.yp <- RMSE.yp %>% c("." %>% paste(mean(PSIprd_yp.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

RMSE.cp.psi <- RMSE.cp.psi %>% c("." %>% paste(mean(PSIprd_cp.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))
RMSE.yp.psi <- RMSE.yp.psi %>% c("." %>% paste(mean(PSIprd_yp.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

    # Logistic regression
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model

dat.lr <- dat.lr %>% rbind(exp(PSIprd_lr.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_lr.trend))) %>%
                             data.frame() %>% setNames(names(dat.lr)))

RMSE.lr <- RMSE.lr %>% c("." %>% paste(mean(PSIprd_lr.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

RMSE.lr.psi <- RMSE.lr.psi %>% c("." %>% paste(mean(PSIprd_lr.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

  # Lambda = 1
trnd <- 1 # Can be 0.9, 0.95, 0.98, or 1
setwd(paste("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1/lmb",trnd*100,sep=""))
load("Tr150Pt10Mf1Results.RData") #Adjust to get results for desired model

dat.cp <- dat.cp %>% rbind(exp(PSIprd_cp.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_cp.trend))) %>%
                             data.frame() %>% setNames(names(dat.cp)))
dat.yp <- dat.yp %>% rbind(exp(PSIprd_yp.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_yp.trend))) %>%
                             data.frame() %>% setNames(names(dat.yp)))

RMSE.cp <- RMSE.cp %>% c("." %>% paste(mean(PSIprd_cp.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))
RMSE.yp <- RMSE.yp %>% c("." %>% paste(mean(PSIprd_yp.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

RMSE.cp.psi <- RMSE.cp.psi %>% c("." %>% paste(mean(PSIprd_cp.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))
RMSE.yp.psi <- RMSE.yp.psi %>% c("." %>% paste(mean(PSIprd_yp.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

    # Logistic regression
load("Tr300Pt10Mf1Results_lr.RData") #Adjust to get results for desired model

dat.lr <- dat.lr %>% rbind(exp(PSIprd_lr.trend[,1:3]) %>% cbind(rep(trnd,nrow(PSIprd_lr.trend))) %>%
                             data.frame() %>% setNames(names(dat.lr)))

RMSE.lr <- RMSE.lr %>% c("." %>% paste(mean(PSIprd_lr.trend[, "MSE"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

RMSE.lr.psi <- RMSE.lr.psi %>% c("." %>% paste(mean(PSIprd_lr.trend[, "MSE.psi"]) %>%
            sqrt() %>% round(digits = 3) %>% as.character() %>% substr(3,5),sep=""))

dat.cp$X <- dat.yp$X <- dat.lr$X <- dat.cp$lmbN + runif(nrow(dat.cp),-0.005,0.005) # jitter for plotting

## Estimated occupancy vs true abundance trends
pcp <- ggplot(data = dat.cp,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.003,alpha=0.3) +
  geom_point(x=0.9,y=0.9,size=6,color="red") + # Add accuracy targets for true abundance trends
  geom_point(x=0.95,y=0.95,size=6,color="red") +
  geom_point(x=0.98,y=0.98,size=6,color="red") +
  geom_point(x=1,y=1,size=6,color="red") +
  geom_point(x=0.9,y=0.9525030,size=6,color="blue") + # Add accuracy targets for true occupancy trends
  geom_point(x=0.95,y=0.9854334,size=6,color="blue") +
  geom_point(x=0.98,y=0.9963604,size=6,color="blue") +
  scale_x_continuous(breaks=c(0.9,0.95,0.98,1)) +
  scale_y_continuous(breaks=c(0.8,0.9,1)) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.text.x=element_text(size=26)) +
  theme(axis.title.x=element_text(size=35)) +
  theme(axis.title.y=element_text(size=35)) +
  theme(axis.text.y=element_text(size=26)) +
  geom_text(aes(x=0.9,y=1.082),label="A",size=10) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.062,1.062,1.062,1.062),label=RMSE.cp.psi,size=c(8,8,8,8)) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.042,1.042,1.042,1.042),label=RMSE.cp,size=c(8,8,8,8))
  
#  geom_text(aes(x=0.9,y=1.142),label="A",size=10) +
#  annotate("text", x=0.897, y=1.11, label = "RMSE[psi]*' = '", size=8, hjust=0, parse = T) +
#  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.09,1.09,1.09,1.09),label=RMSE.cp.psi,size=c(8,8,8,8)) +
#  annotate("text", x=0.897, y=1.062, label = "RMSE[N]*' = '", size=8, hjust=0, parse = T) +
#  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.042,1.042,1.042,1.042),label=RMSE.cp,size=c(8,8,8,8))

pyp <- ggplot(data = dat.yp,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.003,alpha=0.3) +
  geom_point(x=0.9,y=0.9,size=6,color="red") + # Add accuracy targets for true occupancy trends
  geom_point(x=0.95,y=0.95,size=6,color="red") +
  geom_point(x=0.98,y=0.98,size=6,color="red") +
  geom_point(x=1,y=1,size=6,color="red") +
  geom_point(x=0.9,y=0.9525030,size=6,color="blue") + # Add accuracy targets for true occupancy trends
  geom_point(x=0.95,y=0.9854334,size=6,color="blue") +
  geom_point(x=0.98,y=0.9963604,size=6,color="blue") +
  scale_x_continuous(breaks=c(0.9,0.95,0.98,1)) +
  scale_y_continuous(breaks=c(0.8,0.9,1,1.1)) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.text.x=element_text(size=26)) +
  theme(axis.title.x=element_text(size=35)) +
  theme(axis.title.y=element_text(size=35)) +
  theme(axis.text.y=element_text(size=26)) +
  geom_text(aes(x=0.9,y=1.17),label="B",size=10) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.15,1.15,1.15,1.15),label=RMSE.yp.psi,size=c(8,8,8,8)) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.13,1.13,1.13,1.13),label=RMSE.yp,size=c(8,8,8,8))

#  geom_text(aes(x=0.9,y=1.2),label="B",size=10) +
#  annotate("text", x=0.897, y=1.17, label = "RMSE[psi]*' = '", size=8, hjust=0, parse = T) +
#  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.15,1.15,1.15,1.15),label=RMSE.yp.psi,size=c(8,8,8,8)) +
#  annotate("text", x=0.897, y=1.12, label = "RMSE[N]*' = '", size=8, hjust=0, parse = T) +
#  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.1,1.1,1.1,1.1),label=RMSE.yp,size=c(8,8,8,8))

plr <- ggplot(data = dat.lr,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.003,alpha=0.3) +
  geom_point(x=0.9,y=0.9,size=6,color="red") + # Add accuracy targets for true occupancy trends
  geom_point(x=0.95,y=0.95,size=6,color="red") +
  geom_point(x=0.98,y=0.98,size=6,color="red") +
  geom_point(x=1,y=1,size=6,color="red") +
  geom_point(x=0.9,y=0.9525030,size=6,color="blue") + # Add accuracy targets for true occupancy trends
  geom_point(x=0.95,y=0.9854334,size=6,color="blue") +
  geom_point(x=0.98,y=0.9963604,size=6,color="blue") +
  scale_x_continuous(breaks=c(0.9,0.95,0.98,1)) +
  scale_y_continuous(breaks=c(0.9,1)) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.text.x=element_text(size=26)) +
  theme(axis.title.x=element_text(size=35)) +
  theme(axis.title.y=element_text(size=35)) +
  theme(axis.text.y=element_text(size=26)) +
  geom_text(aes(x=0.9,y=1.05),label="C",size=10) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.04,1.04,1.04,1.04),label=RMSE.lr.psi,size=c(8,8,8,8)) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.03,1.03,1.03,1.03),label=RMSE.lr,size=c(8,8,8,8))

#  geom_text(aes(x=0.9,y=1.082),label="C",size=10) +
#  annotate("text", x=0.897, y=1.066, label = "RMSE[psi]*' = '", size=8, hjust=0, parse = T) +
#  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.056,1.056,1.056,1.056),label=RMSE.lr.psi,size=c(8,8,8,8)) +
#  annotate("text", x=0.897, y=1.04, label = "RMSE[N]*' = '", size=8, hjust=0, parse = T) +
#  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.03,1.03,1.03,1.03),label=RMSE.lr,size=c(8,8,8,8))

p <- ggdraw() + 
    draw_plot(pcp, x = 0.07, y = 0.525, width = 0.465, height = 0.475) +
    draw_plot(pyp, x = 0.535, y = 0.525, width = 0.465, height = 0.475) +
    draw_plot(plr, x = 0.07, y = 0.05, width = 0.465, height = 0.475) +
    draw_plot_label(label = c("lambda[N]", "hat(lambda)[psi]", "RMSE[psi]*' = '",
                              "RMSE[N]*' = '",
                              "RMSE[psi]*' = '",
                              "RMSE[N]*' = '"),
                    size=c(40, 40, 20, 20, 20, 20),
                    x=c(0.5, 0.01, 0.02, 0.02, 0.02, 0.02),
                    y=c(0.06, 0.6, 0.975, 0.945, 0.5, 0.475),
                    hjust = c(0, 0, 0, 0, 0, 0), parse = T)
#p

setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/")
save_plot("manuscript/Figure_trend_estimates.tiff", p, ncol = 3, nrow = 3.5, dpi = 50)

#### Abridged version for presentation ####
## Estimated occupancy vs true abundance trends
pcp <- ggplot(data = dat.cp,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.003,alpha=0.3) +
  geom_point(x=0.9,y=0.9,size=6,color="red") + # Add accuracy targets for true abundance trends
  geom_point(x=0.95,y=0.95,size=6,color="red") +
  geom_point(x=0.98,y=0.98,size=6,color="red") +
  geom_point(x=1,y=1,size=6,color="red") +
  scale_x_continuous(breaks=c(0.9,0.95,0.98,1)) +
  scale_y_continuous(breaks=c(0.8,0.9,1)) +
  ylab(NULL) + xlab(expression(lambda[N])) +
  theme(axis.text.x=element_text(size=26)) +
  theme(axis.title.x=element_text(size=40)) +
  theme(axis.text.y=element_text(size=26)) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.042,1.042,1.042,1.042),label=RMSE.cp,size=c(8,8,8,8))

plr <- ggplot(data = dat.lr,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.003,alpha=0.3) +
  geom_point(x=0.9,y=0.9,size=6,color="red") + # Add accuracy targets for true occupancy trends
  geom_point(x=0.95,y=0.95,size=6,color="red") +
  geom_point(x=0.98,y=0.98,size=6,color="red") +
  geom_point(x=1,y=1,size=6,color="red") +
  scale_x_continuous(breaks=c(0.9,0.95,0.98,1)) +
  scale_y_continuous(breaks=c(0.9,1)) +
  ylab(NULL) + xlab(expression(lambda[N])) +
  theme(axis.text.x=element_text(size=26)) +
  theme(axis.title.x=element_text(size=40)) +
  theme(axis.text.y=element_text(size=26)) +
  annotate("text",x=c(0.9,0.95,0.98,1),y=c(1.03,1.03,1.03,1.03),label=RMSE.lr,size=c(8,8,8,8))

p <- ggdraw() + 
  draw_plot(pcp, x = 0.05, y = 0, width = 0.465, height = 1) +
  draw_plot(plr, x = 0.525, y = 0, width = 0.465, height = 1) +
  draw_plot_label(label = c("hat(lambda)[psi]", "RMSE*' = '"),
                  size=c(40, 20),
                  x=c(0.01, 0.02),
                  y=c(0.7, 0.98),
                  hjust = c(0, 0), parse = T)
#p

setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/")
save_plot("manuscript/Figure_trend_estimates_MTTWS_2017.tiff", p, ncol = 3.5, nrow = 1.75, dpi = 50)
