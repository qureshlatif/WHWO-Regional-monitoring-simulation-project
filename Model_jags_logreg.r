model { 
  # prior distributions
  for (t in 1:n.yrs) {
    B0[t] ~ dnorm(0,0.1)T(-10,10) # Transect scale logit occupancy for points
    logit(PSI[t]) <- B0[t]
  }
  
  for (i in 1:n.trns) {
    for (t in 1:n.yrs) {
      Y.lgreg[i,t] ~ dbin(PSI[t],1) # OBSERVATION MODEL
    }
  }
}
