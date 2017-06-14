library(ggplot2)
library(cowplot)
library(grid)
setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6")

######### Tabulate power for detecting population trends ##########
library(R.utils)

base.folder <- "HR1000mD2522/Obj2/"
scn <- 1:28
ntrns <- c(60,120,180,75,150,225,120,240,360,200,400,600,60,67,80,90,120,240,360,150,300,450,240,480,
           720,400,800,1200)
PrpSurvPerYr <- c(rep(seq(1,3),4),rep(1,4),rep(seq(1,3),4)) #1 = all, 2 = half, 3 = third
PrpSurv2x <- c(rep(1,12),c(1,0.8,0.5,0.33),rep(1,12)) #Proportion of transects surveyed twice per year
npnts <- c(rep(10,3),rep(8,3),rep(5,3),rep(3,3),rep(10,4),
           rep(10,3),rep(8,3),rep(5,3),rep(3,3))
model <- c(rep("occupancy",16),rep("glm binom",12))

out <- data.frame(cbind(Scenario=scn,ntrns,PrpSurvPerYr,PrpSurv2x,npnts))
out$model <- model
out$RMSE <- out$power <- (NA)*1
trend.estimates <- array(NA,c(100,3,length(scn)))
rm(scn,ntrns,PrpSurvPerYr,PrpSurv2x,npnts,model) #Cleanup

for(i in (1:nrow(out))) {
  if(out$model[i]=="occupancy") {
    load(paste(base.folder,"cp/Tr",out$ntrns[i],"Pt",out$npnts[i],
               "Mf",out$PrpSurvPerYr[i],"P2V",out$PrpSurv2x[i]*100,
               "Results.RData",sep=""))
    nsim <- dim(PSIprd_cp.trend)[1]
    out$power[i] <- (length(which(PSIprd_cp.trend[,3]<0))/nsim)*100
    trend.estimates[1:nsim,,i] <- exp(PSIprd_cp.trend[,1:3])
    out$RMSE[i] <- sqrt(mean(PSIprd_cp.trend[,4]))
  }
  if(out$model[i]=="glm binom") {
    load(paste(base.folder,"np/Tr",out$ntrns[i],"Pt",out$npnts[i],
               "Mf",out$PrpSurvPerYr[i],"Results.RData",sep=""))
    nsim <- dim(PSIprd_lr.trend)[1]
    out$power[i] <- (length(which(PSIprd_lr.trend[,3]<0))/nsim)*100
    trend.estimates[1:nsim,,i] <- exp(PSIprd_lr.trend[,1:3])
    out$RMSE[i] <- sqrt(mean(PSIprd_lr.trend[,4]))
  }
}

######## Plot power for panel design and variable points per transect ##########

## Constant-p model
dat <- out[1:12,]
dat$npnts <- as.factor(dat$npnts)
dat$PrpSurvPerYr <- rep(c(100,50,33),4)
dat$PrpSurvPerYr <- as.factor(dat$PrpSurvPerYr)
dat$ntrns <- paste("(",dat$ntrns,")",sep="")

dodge=position_dodge(width=0.9)
pbar <- ggplot(data = dat,aes(x=npnts,y=power,fill=PrpSurvPerYr)) +
  geom_bar(aes(fill=PrpSurvPerYr),stat="identity",colour="black",position=dodge) +
  geom_hline(yintercept=80,linetype="dotted",size=1) +
  ylab("Power") + xlab("Number of points per transect") +
  scale_fill_manual(values = c("white","gray","black")) +
  scale_y_continuous(breaks=seq(0,100,by=20)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.title.y=element_text(size=30)) +
  geom_text(aes(label=ntrns),vjust=-0.3,size=8,position=dodge)+
  labs(fill="% transects \nsurveyed \nper year") +
  theme(legend.title=element_text(size=20)) +
  theme(legend.text=element_text(size=18))

dat.trests <- trend.estimates[,,1:12]
x <- c(0.68,1,1.32,1.68,2,2.32,2.68,3,3.32,3.68,4,4.32)
ind <- replace(as.character(dat$PrpSurvPerYr),which(is.element(dat$PrpSurvPerYr,c("33","50"))),
                paste("0",dat$PrpSurvPerYr[which(is.element(as.character(dat$PrpSurvPerYr),c("33","50")))],sep=""))
ind <- paste(dat$npnts,ind)
ind[4:12] <- paste("0",ind[4:12],sep="")
x <- x[order(ind)]
dat.plt <- data.frame(cbind(dat.trests[,,1],rep(x[1],100)))
for(i in 2:12) dat.plt <- rbind(dat.plt,setNames(data.frame(cbind(dat.trests[,,i],rep(x[i],100))),
                                                 names(dat.plt)))
rm(ind)
names(dat.plt) <- c("med","lo","hi","index")
dat.plt$X <- dat.plt$index + runif(nrow(dat.plt),-0.05,0.05) # jitter for plotting

RMSE <- substr(as.character(round(out$RMSE[1:12],digits=3)),2,5)

ptrd <- ggplot(data = dat.plt,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0,alpha=0.3) +
  geom_hline(yintercept=0.98,linetype="dotted",size=1) +
  scale_y_continuous(lim=c(0.82,1.16),breaks=c(0.9,1,1.1)) +
  scale_x_continuous(lim=c(0.59,4.4),breaks=c()) +
  ylab(expression(hat(lambda)[psi])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=35,angle=0)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=0.6,y=1.16),label="RMSE = ",size=10,hjust=0) +
  annotate("text",x=unique(dat.plt$index),y=rep(1.122,length(x)),label=RMSE,size=rep(10,length(x)))

theme_set(theme_bw())
p <- ggdraw() +
  draw_plot(pbar, x = 0.011, y = 0, width = 0.989, height = 0.65) + 
  draw_plot(ptrd, x = 0, y = 0.65, width = 0.878, height = 0.35) +
  draw_plot_label("lambda[N]",size=25,x=0.865,y=0.84,parse=T)

save_plot("manuscript/Figure_npnts_&_PanelDesign_cp.tiff", p, ncol = 3.5, nrow = 3.5, dpi = 50)

## Logistic regression (apparent occupancy)
dat <- out[17:28,]
dat$npnts <- as.factor(dat$npnts)
dat$PrpSurvPerYr <- rep(c(100,50,33),4)
dat$PrpSurvPerYr <- as.factor(dat$PrpSurvPerYr)
dat$ntrns <- paste("(",dat$ntrns,")",sep="")

dodge=position_dodge(width=0.9)
pbar <- ggplot(data = dat,aes(x=npnts,y=power,fill=PrpSurvPerYr)) +
  geom_bar(aes(fill=PrpSurvPerYr),stat="identity",colour="black",position=dodge) +
  geom_hline(yintercept=80,linetype="dotted",size=1) +
  ylab("Power") + xlab("Number of points per transect") +
  scale_fill_manual(values = c("white","gray","black")) +
  scale_y_continuous(breaks=seq(0,100,by=20)) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.title.y=element_text(size=30)) +
  geom_text(aes(label=ntrns),vjust=-0.3,size=8,position=dodge)+
  labs(fill="% transects \nsurveyed \nper year") +
  theme(legend.title=element_text(size=20)) +
  theme(legend.text=element_text(size=18))

dat.trests <- trend.estimates[,,17:28]
x <- c(0.68,1,1.32,1.68,2,2.32,2.68,3,3.32,3.68,4,4.32)
ind <- replace(as.character(dat$PrpSurvPerYr),which(is.element(dat$PrpSurvPerYr,c("33","50"))),
                paste("0",dat$PrpSurvPerYr[which(is.element(as.character(dat$PrpSurvPerYr),c("33","50")))],sep=""))
ind <- paste(dat$npnts,ind)
ind[4:12] <- paste("0",ind[4:12],sep="")
x <- x[order(ind)]
dat.plt <- data.frame(cbind(dat.trests[,,1],rep(x[1],100)))
for(i in 2:12) dat.plt <- rbind(dat.plt,setNames(data.frame(cbind(dat.trests[,,i],rep(x[i],100))),
                                                 names(dat.plt)))
rm(ind)
names(dat.plt) <- c("med","lo","hi","index")
dat.plt$X <- dat.plt$index + runif(nrow(dat.plt),-0.05,0.05) # jitter for plotting

RMSE <- substr(as.character(round(out$RMSE[17:28],digits=3)),2,5)

ptrd <- ggplot(data = dat.plt,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0,alpha=0.3) +
  geom_hline(yintercept=0.98,linetype="dotted",size=1) +
  scale_y_continuous(lim=c(0.82,1.16),breaks=c(0.9,1,1.1)) +
  scale_x_continuous(lim=c(0.59,4.4),breaks=c()) +
  ylab(expression(hat(lambda)[psi])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=35,angle=0)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=0.6,y=1.16),label="RMSE = ",size=10,hjust=0) +
  annotate("text",x=unique(dat.plt$index),y=rep(1.122,length(x)),label=RMSE,size=rep(10,length(x)))

theme_set(theme_bw())
p <- ggdraw() +
  draw_plot(pbar, x = 0.011, y = 0, width = 0.989, height = 0.65) + 
  draw_plot(ptrd, x = 0, y = 0.65, width = 0.878, height = 0.35) +
  draw_plot_label("lambda[N]",size=25,x=0.865,y=0.84,parse=T)

save_plot("manuscript/Figure_npnts_&_PanelDesign_logreg.tiff", p, ncol = 3.5, nrow = 3.5, dpi = 50)

############# Plot power for reduced repeat-visits ###############
dat <- out[13:16,]
dat$PrpSurv2x <- dat$PrpSurv2x*100
dat$PrpSurv2x <- as.factor(dat$PrpSurv2x)
dat$ntrns <- paste("(",dat$ntrns,")",sep="")
dat$no_surveys <- as.factor(rep(300,4))

pbar <- ggplot(data = dat,aes(x=PrpSurv2x,y=power)) +
  geom_bar(stat="identity",colour="black",fill="gray") +
  geom_hline(yintercept=80,linetype="dotted",size=1) +
  ylab("Power") + xlab("% transects surveyed twice (total transects surveyed)") +
  scale_x_discrete(labels=paste(as.character(dat$PrpSurv2x),dat$ntrns)[order(dat$PrpSurv2x)]) +
  scale_y_continuous(breaks=seq(0,100,by=20)) +
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.text.y=element_text(size=25)) +
  theme(axis.title.x=element_text(size=35)) +
  theme(axis.title.y=element_text(size=35))

dat.trests <- trend.estimates[,,13:16]
x <- 1:4
ind <- replace(as.character(dat$PrpSurv2x),which(is.element(dat$PrpSurv2x,c("33","50","80"))),
                paste("0",dat$PrpSurv2x[which(is.element(as.character(dat$PrpSurv2x),c("33","50","80")))],sep=""))
x <- x[order(ind)]
dat.plt <- data.frame(cbind(dat.trests[,,1],rep(x[1],100)))
for(i in 2:4) dat.plt <- rbind(dat.plt,setNames(data.frame(cbind(dat.trests[,,i],rep(x[i],100))),
                                                 names(dat.plt)))
names(dat.plt) <- c("med","lo","hi","index")
dat.plt$X <- dat.plt$index + runif(nrow(dat.plt),-0.2,0.2) # jitter for plotting

RMSE <- substr(as.character(round(out$RMSE[13:16],digits=3)),2,5)

ptrd <- ggplot(data = dat.plt,aes(x=X,y=med)) +
  geom_point(size=3,alpha=0.3) +
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0,alpha=0.3) +
  geom_hline(yintercept=0.98,linetype="dotted",size=1) +
  scale_y_continuous(lim=c(0.77,1.17),breaks=c(0.8,0.9,1,1.1)) +
  scale_x_continuous(lim=c(0.6,4.4),breaks=c()) +
  ylab(expression(hat(lambda)[psi])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=35,angle=0)) +
  theme(axis.text.y=element_text(size=26)) +
  geom_text(aes(x=0.9,y=1.17),label="RMSE = ",size=10,hjust=0) +
  annotate("text",x=unique(dat.plt$index),y=rep(1.14,length(x)),label=RMSE,size=rep(10,length(x)))

theme_set(theme_bw())
p <- ggdraw() +
  draw_plot(pbar, x = 0.02, y = 0, width = 0.94, height = 0.65) + 
  draw_plot(ptrd, x = 0, y = 0.65, width = 0.96, height = 0.35) +
  draw_plot_label("lambda[N]",size=25,x=0.95,y=0.855,parse=T)

save_plot("manuscript/Figure_vary_PrpRepVisits.tiff", p, ncol = 3.5, nrow = 3.5, dpi=50) #thumbnail
