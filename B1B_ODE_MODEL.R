B1B_model <- function (t, y, interactionMat, growthRate_vector, numTaxa) {
  
  growthRate_vector = parms.growthRateVector
  interactionMat    = parms.interactionMatrix
  numTaxa           = parms.numTaxa
  
  abundancevectemp=c()
  for(k in 1:numTaxa){
    abundancevectemp[k] = y[k];
  }
  
  dydt = c()
  
  for(tx in 1:numTaxa){
    interactionvectemp = interactionMat[tx,1:numTaxa];
    interactions_tx    = interactionvectemp%*%abundancevectemp; 
    dydt[tx]           = y[tx]*(growthRate_vector[tx] + interactions_tx);
  }
  
  ## return result as a list!
  list(dydt)
}

parms.growthRateVector=growthRateVector
parms.interactionMatrix=as.matrix(interactionMatrix)
parms.numTaxa=numTaxa

times <- seq(from=0,to=t_end,by=1)
y0    <- y0_meanSubjects

library(tidyverse)

ode(
  func=B1B_model,
  y=y0,
  times=times,
  parms=parms,
  method='bdf'
) %>% as.data.frame() -> out

colnames(out)=c('day',taxa_array)

out_longer = out %>% pivot_longer(!day,names_to = 'taxa', values_to = 'abundance')
out_longer$type='simulation'
out_longer_day0 = out_longer %>% filter(day==0)


