library(R.utils)
library(abind)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

### Plot true abundance, occupancy, and detectability ###

setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/HR1000mD2522/obj1")
odds <- function(x) x/(1-x)

trueN <- abind(loadObject("lmb90/TrueN"),loadObject("lmb95/TrueN"),loadObject("lmb98/TrueN"),along=3)

true_psi <- abind(loadObject("lmb90/True_psi_10PntsPerTrns")[,,"occtru"],
                  loadObject("lmb95/True_psi_10PntsPerTrns")[,,"occtru"],
                  loadObject("lmb98/True_psi_10PntsPerTrns")[,,"occtru"],along=3)

## True abundance vs. true proportion occupancy ##
dat <- data.frame(cbind(as.numeric(trueN),as.numeric(true_psi)))
names(dat) <- c("N","Psi")
pPsiN_10pt <- ggplot(data = dat,aes(x=N,y=Psi)) +
  geom_point(size=3,alpha=0.2,color="black") +
  ylab(expression(psi["10-point transects"])) + xlab(NULL) +
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.title.y=element_text(size=35)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=3700,y=1),label="A",size=15)

true_psi <- abind(loadObject("lmb90/True_psi_3PntsPerTrns")[,,"occtru"],
                  loadObject("lmb95/True_psi_3PntsPerTrns")[,,"occtru"],
                  loadObject("lmb98/True_psi_3PntsPerTrns")[,,"occtru"],along=3)

dat <- data.frame(cbind(as.numeric(trueN),as.numeric(true_psi)))
names(dat) <- c("N","Psi")
pPsiN_3pt <- ggplot(data = dat,aes(x=N,y=Psi)) +
  geom_point(size=3,alpha=0.2,color="black") +
  ylab(expression(psi["3-point transects"])) + xlab("N (across all forests)") +
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.title.y=element_text(size=35)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=3700,y=1),label="C",size=15)

## Abundance vs. Occupancy lambdas ##
lmbN <- c(apply(trueN[1:30,2:20,1]/trueN[1:30,1:19,1],1,mean),
          apply(trueN[1:30,2:20,2]/trueN[1:30,1:19,2],1,mean),
          apply(trueN[1:30,2:20,3]/trueN[1:30,1:19,3],1,mean))
true_psi <- abind(loadObject("lmb90/True_psi_10PntsPerTrns")[,,"occtru"],
                  loadObject("lmb95/True_psi_10PntsPerTrns")[,,"occtru"],
                  loadObject("lmb98/True_psi_10PntsPerTrns")[,,"occtru"],along=3)
lmbPsi <- c(apply(true_psi[1:30,2:20,1]/true_psi[1:30,1:19,1],1,mean,na.rm=T),
          apply(true_psi[1:30,2:20,2]/true_psi[1:30,1:19,2],1,mean,na.rm=T),
          apply(true_psi[1:30,2:20,3]/true_psi[1:30,1:19,3],1,mean,na.rm=T))

dat_10pt <- data.frame(cbind(lmbN,lmbPsi))
for(i in 1:30) dat <- rbind(dat, c(1, 1)) # Add representation of no-trend scenario.

dseg <- data.frame(x = 0.9, y = 0.9, xe = 1, ye = 1)
pLmb_10pt <- ggplot(data = dat_10pt,aes(x=lmbN,y=lmbPsi)) +
  geom_segment(aes(x = x, y = y, xend = xe, yend = ye), data = dseg, color = "#D55E00", size = 2) +
  geom_point(size=3, alpha=0.2, color="black") +
  ylab(expression(lambda[psi["10-point transects"]])) + xlab(NULL) +
  scale_x_continuous(breaks=c(0.9, 0.95, 0.98, 1)) +
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.title.y=element_text(size=35)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=0.9,y=1),label="B",size=15)
#pLmb

lmbN <- c(apply(trueN[1:30,2:20,1]/trueN[1:30,1:19,1],1,mean),
          apply(trueN[1:30,2:20,2]/trueN[1:30,1:19,2],1,mean),
          apply(trueN[1:30,2:20,3]/trueN[1:30,1:19,3],1,mean))
true_psi <- abind(loadObject("lmb90/True_psi_3PntsPerTrns")[,,"occtru"],
                  loadObject("lmb95/True_psi_3PntsPerTrns")[,,"occtru"],
                  loadObject("lmb98/True_psi_3PntsPerTrns")[,,"occtru"],along=3)
lmbPsi <- c(apply(true_psi[1:30,2:20,1]/true_psi[1:30,1:19,1],1,mean,na.rm=T),
            apply(true_psi[1:30,2:20,2]/true_psi[1:30,1:19,2],1,mean,na.rm=T),
            apply(true_psi[1:30,2:20,3]/true_psi[1:30,1:19,3],1,mean,na.rm=T))

dat_3pt <- data.frame(cbind(lmbN,lmbPsi))
for(i in 1:30) dat <- rbind(dat, c(1, 1)) # Add representation of no-trend scenario.

dseg <- data.frame(x = 0.9, y = 0.9, xe = 1, ye = 1)
pLmb_3pt <- ggplot(data = dat_3pt,aes(x=lmbN,y=lmbPsi)) +
  geom_segment(aes(x = x, y = y, xend = xe, yend = ye), data = dseg, color = "#D55E00", size = 2) +
  geom_point(size=3, alpha=0.2, color="black") +
  ylab(expression(lambda[psi["3-point transects"]])) + xlab(expression(lambda[N])) +
  scale_x_continuous(breaks=c(0.9, 0.95, 0.98, 1)) +
  theme(axis.text.x=element_text(size=25)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.title.y=element_text(size=35)) +
  theme(axis.text.y=element_text(size=25)) +
  geom_text(aes(x=0.9,y=1),label="D",size=15)

p <- ggdraw() + 
  draw_plot(pPsiN_10pt, x = 0, y = 0.525, width = 0.49, height = 0.475) +
  draw_plot(pLmb_10pt, x = 0.49, y = 0.525, width = 0.49, height = 0.475) +
  draw_plot(pPsiN_3pt, x = 0, y = 0, width = 0.49, height = 0.525) +
  draw_plot(pLmb_3pt, x = 0.49, y = 0, width = 0.49, height = 0.525)

setwd("F:/research stuff/FS_PostDoc/Occupancy_analysis_simulations/WHWO_R6_monitoring/Power_analysis_via_simulation/R6/")
save_plot("manuscript/Figure_truth.tiff", p, ncol = 3, nrow = 3, dpi = 600)

# Mean true occupancy trends #
tapply(dat_10pt$lmbPsi, round(dat_10pt$lmbN, digits = 2), mean)

# 0.90      0.95      0.98
# 0.9525030 0.9854334 0.9963604

tapply(dat_3pt$lmbPsi, round(dat_3pt$lmbN, digits = 2), mean)

# 0.90      0.95      0.98
# 0.9251990 0.9694622 0.9902002