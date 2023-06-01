rm(list=ls())
library("bayesplot")
library("ggplot2")
library("rstanarm") 
library('rstan')

pathModelOutput='/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/31052023/RDATA';
fileList = list.files(path =pathModelOutput)

##### Pick the file
indexPick     = 3
fileNamePick  = paste0(pathModelOutput,"/",fileList[indexPick])
#################################################

### load the model object - clean up afterwards
temp.space <- new.env()
bar <- load(fileNamePick, temp.space)
loadedModel <- get(bar, temp.space)
rm(temp.space)
#################################################

fit           = loadedModel
fit_summary   = summary(fit)
posterior     = as.array(fit)
list_of_draws = rstan::extract(loadedModel)
paramRoot     = as.name('interactionMat_vector_diag')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]

paramNames=c()
for(i in 1:numOfParams){
  paramNames = c(paramNames,paste0(paramRoot,'[',i,']'))
}

color_scheme_set("brightblue")
mcmc_dens(posterior, pars = paramNames, 
          facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("brightblue")
mcmc_intervals(posterior, pars = paramNames)

color_scheme_set("viridis")
mcmc_trace(posterior, pars = paramNames, 
           facet_args = list(ncol = 1, strip.position = "left"))

