model { 
  # prior distributions
  p~dunif(0,1) # Detection probability (point scale)
  
  for (t in 1:n.yrs) {
    B0[t] ~ dnorm(0,0.1)T(-10,10) # Transect scale logit occupancy for points
    logit(PSI[t]) <- B0[t]
  }
  
  for (i in 1:n.trns) {
    for (t in 1:n.yrs) {
      Z[i,t] ~ dbin(PSI[t],1) # Occupancy state at transects following year 1
      prob.y[i,t] <- Z[i,t]*p     
      Y.mat[i,t] ~ dbin(prob.y[i,t],n.vsts[i,t]) # OBSERVATION MODEL
    }
  }
}
