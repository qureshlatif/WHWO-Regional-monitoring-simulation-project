library(rstan)
library(MASS)

# Set values for simulating data
p <- 0.7   #Values for p at intercept (mean of covariate values)
psi <- 0.4 #Values for psi at intercept (mean of covariate value(s))
n.site <- 500 #no. sites
n.visit <- 2 #no. visits (will try n.visit<-c(5,10,20,40) for n.site = 30 and p = 0.1 later)
n.yrs <- 20

# Function for simulating data
expit <-
  function(x){
    exp(x)/(1+exp(x))
  }

logit <-
  function(x){
    log(x/(1-x))
  }

Dat.sim <- function(p0,psi0,n.sts,n.vsts,n.yrs) {
  Y<-array(NA,dim=c(n.sts,n.yrs,n.vsts))
  psi <- psi0
  p <- p0
  for(ii in 1:dim(Y)[1]) {
    z <- rbinom(n.yrs,1,psi)
    for(jj in 1:dim(Y)[2]) Y[ii,jj,] <- rbinom(n.vsts,1,z[jj]*p)
  }
  return(Y)
}

# Simulate data
Y <- Dat.sim(p,psi,n.site,n.visit,n.yrs)
#Y[1:250,seq(2,20,by=2),] <- NA #To test for robustness to missing values
#Y[251:500,seq(1,19,by=2),] <- NA
R <- dim(Y)[1]
K <- dim(Y)[2]
J <- dim(Y)[3]
Data <- list(Y=Y,R=R,K=K,J=J)

# Set initial values for HMC sampling
inits <- function()
  list(psi=runif(n.yrs),p=runif(1))

#Define which parameters to save.
parameters <- c("psi","p")

mod_occupancy <- '
data {
int<lower=0> R; // number of sites
int<lower=0> J; // number of visits
int<lower=0> K; // number of years
int<lower=0,upper=J> Y[R,K,J]; // detection-nondetection data
}

parameters {
real<lower=0,upper=1> psi[K];
real<lower=0,upper=1> p;
}

model {
// local variables to avoid recomputing log(psi) and log(1-psi)
real log_psi[K];
real log1m_psi[K];
for(k in 1:K) {
  log_psi[k] <- log(psi[k]);
  log1m_psi[k] <- log1m(psi[k]);
}

// priors
for(k in 1:K) {
  psi[k] ~ uniform(0,1);
}
p ~ uniform(0,1);

// likelihood
   for (r in 1:R) { 
    for (k in 1:K) { 
     if (sum(Y[r,k]) > 0) 
       increment_log_prob(log_psi[k] + bernoulli_log(Y[r,k],p)); 
     else 
       increment_log_prob(log_sum_exp(log_psi[k] + bernoulli_log(Y[r,k],p),log1m_psi[k])); 
    } 
   } 
} 
'

# Settings for HMC chains
nc <- 3
nb <- 1000
ni <- 2000 
nt <- 10 

#Run model
out <- stan(model_code = mod_occupancy, data = Data, pars = parameters, init = inits,
            iter = ni, chains = nc, warmup = nb, thin = nt)

#Run model using previously compiled C code
out1 <- stan(fit = out, data = Data, pars=parameters, init=inits, iter = ni, chains = nc, warmup=nb, thin=nt)

library(R.utils)
saveObject





