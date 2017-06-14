# Utility functions used in r6 WHWO monitoring simulations #

Maxtrns <- function(parent.GIS,forest) { ## Finds maximum number of transects allowable in each forest ##
  N <- numeric(length(forest))
  for(f in 1:length(forest)) {
    sets <- read.table(paste(parent.GIS,"/NF_",forest[f],"/","Transect_sets.txt",sep=""),header=T,sep=",")
    N[f] <- max(apply(sets[,2:ncol(sets)],2,sum))
  }
  return(N)
}

logit <- function(x) log(x/(1-x))

expit <- function(x) exp(x)/(1+exp(x))

trend <- function(obs) {
  m <- lm(obs~seq(1,length(obs)))
  return(as.numeric(coefficients(m)[2]))
}
