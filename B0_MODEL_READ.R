rm(list=ls())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
library("bayesplot")
library("ggplot2")
library("rstanarm") 
library('rstan')
library('brms')
library('dplyr')
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html

pathModelOutput='/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/04062023/RDATA';
fileList = list.files(path =pathModelOutput, pattern='.RData')

##### Pick the file
indexPick     = 2
fileNamePick  = paste0(pathModelOutput,"/",fileList[indexPick])
#################################################

### load the model object - clean up afterwards
temp.space <- new.env()
bar <- load(fileNamePick, temp.space)
loadedModel <- get(bar, temp.space)
rm(temp.space)
#################################################

fit           = loadedModel
fitSummary    = summary(fit)
posterior     = as.array(fit)
list_of_draws = rstan::extract(loadedModel)
paramRoot     = as.name('interactionMat_vector_diag')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
n_chains <- fit@sim$chains
n_warmup <- fit@sim$warmup2[1]
n_iter   <-fit@sim$iter[[1]]
print(c(n_iter,n_warmup,n_chains))

paramNames=c()
for(i in 1:numOfParams){
  paramNames = c(paramNames,paste0(paramRoot,'[',i,']'))
}

# color_scheme_set("brightblue")
# mcmc_dens(posterior, pars = paramNames, 
#           facet_args = list(ncol = 1, strip.position = "left"))
# 
# color_scheme_set("brightblue")
# mcmc_intervals(posterior, pars = paramNames)
# 
# color_scheme_set("viridis")
# mcmc_trace(posterior, pars = paramNames, 
#            facet_args = list(ncol = 1, strip.position = "left"))

#################################################
source("SETUP.R")
source('RESHAPE_DATA_KNIGHT.R')
taxa_array = colnames(abundanceArray_meanSubjects)


# compartment_names = 'output_pred'
# summaryTable = as.data.frame(summary(loadedModel,compartment_names)[[1]])
# summaryTable$populationNames = rownames(summaryTable)
# summaryTable$t    = sub(",.*","",sub(".*\\[", "", summaryTable$populationNames))  # Extract characters after pattern
# summaryTable$taxa = sub("\\].*","",sub(".*,", "", summaryTable$populationNames))  # Extract characters after pattern
# summaryTable_use = summaryTable[c('t','taxa','mean','2.5%','97.5%','50%')]
# summaryTable_use$t = as.numeric(summaryTable_use$t)
# summaryTable_use$mean = as.numeric(summaryTable_use$mean)
# summaryTable_use$`97.5%` = as.numeric(summaryTable_use$`97.5%`)
# summaryTable_use$`2.5%` = as.numeric(summaryTable_use$`2.5%`)
# summaryTable_use$`50%` = as.numeric(summaryTable_use$`50%`)
# summaryTable_use$taxa = paste0('y_',summaryTable_use$taxa)
# summaryTable_use = summaryTable_use %>% rowwise() %>% mutate(taxa = taxa_array[as.numeric(sub(",.*","",sub(".*\\_", "", taxa)))] )
# 
# graphics.off()
# ggplot() +
#   geom_point(data=abundanceArray_meanSubjects_longer,aes(x=day,y=abundance,fill=taxa),shape=21,size=1,colour = "black", fill = "white") +
#   # geom_ribbon(data=summaryTable_use,aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=taxa),alpha=.5) +
#   geom_line(data=summaryTable_use,aes(x=t,y=`50%`),colour="black") +
#   facet_wrap(~ taxa ,scales="free",nrow=2) +
#   scale_colour_manual(values=c("grey20","grey"),guide=FALSE) +
#   scale_alpha_manual(values=c(1,0),guide=FALSE) +
#   # scale_y_continuous(trans = 'log10')+
#   labs(x="sample index",y="abundance")


# graphics.off()
# ggplot(summaryTable_use, aes(x = t, y = mean, color = taxa)) + geom_point()

# Save the interaction matricies to compare

paramRoot     = as.name('interactionMat_vector_diag')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
paramNames_diag=c()
for(i in 1:numOfParams){
  paramNames_diag = c(paramNames_diag,paste0(paramRoot,'[',i,']'))
  
}
paramRoot     = as.name('interactionMat_vector_nondiag')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
paramNames_nondiag=c()
for(i in 1:numOfParams){
  paramNames_nondiag = c(paramNames_nondiag,paste0(paramRoot,'[',i,']'))
}
paramRoot     = as.name('growthRate_vector')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
paramNames_growth=c()
for(i in 1:numOfParams){
  paramNames_growth = c(paramNames_growth,paste0(paramRoot,'[',i,']'))
}

summaryTable_nondiag = as.data.frame(summary(loadedModel,paramNames_nondiag)[[1]])
summaryTable_diag    = as.data.frame(summary(loadedModel,paramNames_diag)[[1]])
summaryTable_growth  = as.data.frame(summary(loadedModel,paramNames_growth)[[1]])

# https://mc-stan.org/rstan/reference/Rhat.html
# my Rhats are terrible...

############ INTERACTION MATRIX
interactionMatrix      = matrix(0,nrow = length(taxa_array), ncol = length(taxa_array))
interactionMatrix_rhat = matrix(0,nrow = length(taxa_array), ncol = length(taxa_array))
names_diag    = rownames(summaryTable_diag)
names_nondiag = rownames(summaryTable_nondiag)
idx_diag      = as.numeric(sub("\\].*", "",sub(".*\\[", "", names_diag)))
idx_nondiag   = as.numeric(sub("\\].*", "",sub(".*\\[", "", names_nondiag)))
numCols       = length(taxa_array)
correctIndex  = 0

for (r in 1:length(taxa_array)){
  for (c in 1:length(taxa_array)){
    if(r==c){
      interactionMatrix[r,c]      = -1*summaryTable_diag[which(idx_diag==r),]$mean
      interactionMatrix_rhat[r,c] = summaryTable_diag[which(idx_diag==r),]$Rhat
      correctIndex = correctIndex + 1
    }else{
      idx_nondiag_corrected       = idx_nondiag+correctIndex
      interactionMatrix[r,c]      = summaryTable_nondiag[which(idx_nondiag_corrected==((r-1)*numCols) + c),]$mean
      interactionMatrix_rhat[r,c] = summaryTable_nondiag[which(idx_nondiag_corrected==((r-1)*numCols) + c),]$Rhat
    }
  }
}
colnames(interactionMatrix) = taxa_array
rownames(interactionMatrix) = taxa_array
colnames(interactionMatrix_rhat) = taxa_array
rownames(interactionMatrix_rhat) = taxa_array

library(openxlsx)
ts = as.numeric(sub(".*\\_", "",sub("\\.RData.*", "",fileList[indexPick])))
wb <- createWorkbook()
addWorksheet(wb, "mean")
addWorksheet(wb, "Rhat")
writeData(wb, "mean", as.data.frame(interactionMatrix), startRow = 1, startCol = 1)
writeData(wb, "Rhat", as.data.frame(interactionMatrix_rhat), startRow = 1, startCol = 1)
saveWorkbook(wb, file = paste0(pathModelOutput,"/INTERACTIONS_",ts,".xlsx"), overwrite = TRUE)


############ GROWTH
growthMatrix      = matrix(0,nrow = 1, ncol = length(taxa_array))
growthMatrix_rhat = matrix(0,nrow = 1, ncol = length(taxa_array))
names_growth = rownames(summaryTable_growth)
idx_growth   = as.numeric(sub("\\].*", "",sub(".*\\[", "", names_growth)))

for (c in 1:length(taxa_array)){
  growthMatrix[1,c]      = summaryTable_growth[which(idx_growth==c),]$mean
  growthMatrix_rhat[1,c] = summaryTable_growth[which(idx_growth==c),]$Rhat
}
colnames(growthMatrix) = taxa_array
colnames(growthMatrix_rhat) = taxa_array

library(openxlsx)
ts = as.numeric(sub(".*\\_", "",sub("\\.RData.*", "",fileList[indexPick])))
wb <- createWorkbook()
addWorksheet(wb, "mean")
addWorksheet(wb, "Rhat")
writeData(wb, "mean", as.data.frame(growthMatrix), startRow = 1, startCol = 1)
writeData(wb, "Rhat", as.data.frame(growthMatrix_rhat), startRow = 1, startCol = 1)
saveWorkbook(wb, file = paste0(pathModelOutput,"/GROWTH_",ts,".xlsx"), overwrite = TRUE)











