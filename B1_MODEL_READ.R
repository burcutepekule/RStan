rm(list=ls())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
library("bayesplot")
library("ggplot2")
library("rstanarm") 
library('rstan')
library('brms')
library('dplyr')
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html

pathModelOutput='/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/23062023/RDATA';
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
paramRoot     = as.name('interactionMat_vector_Agnostic')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
n_chains <- fit@sim$chains
n_warmup <- fit@sim$warmup2[1]
n_iter   <-fit@sim$iter[[1]]
print(c(n_iter,n_warmup,n_chains))

# paramNames=c()
# for(i in 1:numOfParams){
#   paramNames = c(paramNames,paste0(paramRoot,'[',i,']'))
# }
# 
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
source('RESHAPE_DATA_KNIGHT_L6.R')
taxa_array = colnames(abundanceArray_meanSubjects)

abundanceArray_meanSubjects_use = abundanceArray_meanSubjects
abundanceArray_meanSubjects_use$day = rownames(abundanceArray_meanSubjects_use)
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects_use %>% pivot_longer(!day,names_to = 'taxa', values_to = 'abundance')

compartment_names = 'output_pred'
summaryTable = as.data.frame(summary(loadedModel,compartment_names)[[1]])
summaryTable$populationNames = rownames(summaryTable)
summaryTable$t    = sub(",.*","",sub(".*\\[", "", summaryTable$populationNames))  # Extract characters after pattern
summaryTable$taxa = sub("\\].*","",sub(".*,", "", summaryTable$populationNames))  # Extract characters after pattern
summaryTable_use = summaryTable[c('t','taxa','mean','2.5%','97.5%','50%')]
summaryTable_use$t = as.numeric(summaryTable_use$t)
summaryTable_use$mean = as.numeric(summaryTable_use$mean)
summaryTable_use$`97.5%` = as.numeric(summaryTable_use$`97.5%`)
summaryTable_use$`2.5%` = as.numeric(summaryTable_use$`2.5%`)
summaryTable_use$`50%` = as.numeric(summaryTable_use$`50%`)
summaryTable_use$taxa = paste0('y_',summaryTable_use$taxa)
summaryTable_use = summaryTable_use %>% rowwise() %>% mutate(taxa = taxa_array[as.numeric(sub(",.*","",sub(".*\\_", "", taxa)))] )

# Save the interaction matricies to compare

paramRoot     = as.name('interactionMat_vector_Agnostic')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
paramNames_agnostic=c()
for(i in 1:numOfParams){
  paramNames_agnostic = c(paramNames_agnostic,paste0(paramRoot,'[',i,']'))
  
}
paramRoot     = as.name('interactionMat_vector_Gnostic')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
paramNames_gnostic=c()
for(i in 1:numOfParams){
  paramNames_gnostic = c(paramNames_gnostic,paste0(paramRoot,'[',i,']'))
}
paramRoot     = as.name('growthRate_vector')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
paramNames_growth=c()
for(i in 1:numOfParams){
  paramNames_growth = c(paramNames_growth,paste0(paramRoot,'[',i,']'))
}
paramRoot     = as.name('phi')
numOfParams   = dim(list_of_draws[[paramRoot]])[2]
paramNames_phi=c()
for(i in 1:numOfParams){
  paramNames_phi = c(paramNames_phi,paste0(paramRoot,'[',i,']'))
}

summaryTable_agnostic = as.data.frame(summary(loadedModel,paramNames_agnostic)[[1]])
summaryTable_gnostic  = as.data.frame(summary(loadedModel,paramNames_gnostic)[[1]])
summaryTable_growth   = as.data.frame(summary(loadedModel,paramNames_growth)[[1]])
summaryTable_phi      = as.data.frame(summary(loadedModel,paramNames_phi)[[1]])


graphics.off()
ggplot() +
  geom_point(data=abundanceArray_meanSubjects_longer,aes(x=day,y=abundance,fill=taxa),shape=21,size=1,colour = "black", fill = "white") +
  # geom_ribbon(data=summaryTable_use,aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=taxa),alpha=.5) +
  geom_line(data=summaryTable_use,aes(x=t,y=`50%`),colour="black") +
  facet_wrap(~ taxa ,scales="free",nrow=2) +
  scale_colour_manual(values=c("grey20","grey")) +
  scale_alpha_manual(values=c(1,0)) +
  # scale_y_continuous(trans = 'log10')+
  labs(x="sample index",y="abundance")


# graphics.off()
# ggplot(summaryTable_use, aes(x = t, y = mean, color = taxa)) + geom_point()


# https://mc-stan.org/rstan/reference/Rhat.html
# my Rhats are terrible...

############ INTERACTION MATRIX
interactionMatrix      = matrix(0,nrow = length(taxa_array), ncol = length(taxa_array))
interactionMatrix_rhat = matrix(0,nrow = length(taxa_array), ncol = length(taxa_array))
names_gnostic   = rownames(summaryTable_gnostic)
names_agnostic  = rownames(summaryTable_agnostic)
idx_gnostic     = as.numeric(sub("\\].*", "",sub(".*\\[", "", names_gnostic)))
idx_agnostic    = as.numeric(sub("\\].*", "",sub(".*\\[", "", names_agnostic)))
numCols         = length(taxa_array)
correctIndex    = 0

counter_mask    = 1;
counter_agnostic= 1;
counter_gnostic = 1;

numTaxa = length(taxa_array)
scale   = 1

for(r in 1:numTaxa){
  for(c in 1:numTaxa){
    mask = interactionMask_vector[counter_mask];
    if(mask==0){
      interactionMatrix[r,c]      = scale*summaryTable_agnostic[counter_agnostic,]$mean
      interactionMatrix_rhat[r,c] = summaryTable_agnostic[counter_agnostic,]$Rhat
      counter_agnostic    = counter_agnostic + 1;
    }else{
      interactionMatrix[r,c]      = mask*scale*summaryTable_gnostic[counter_gnostic,]$mean
      interactionMatrix_rhat[r,c] = summaryTable_gnostic[counter_gnostic,]$Rhat
      counter_gnostic     = counter_gnostic + 1;
    }
    counter_mask = counter_mask + 1;
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

############ PHI
phiMatrix      = matrix(0,nrow = 1, ncol = length(taxa_array))
phiMatrix_rhat = matrix(0,nrow = 1, ncol = length(taxa_array))
names_phi = rownames(summaryTable_phi)
idx_phi   = as.numeric(sub("\\].*", "",sub(".*\\[", "", names_phi)))

for (c in 1:length(taxa_array)){
  phiMatrix[1,c]      = summaryTable_phi[which(idx_phi==c),]$mean
  phiMatrix_rhat[1,c] = summaryTable_phi[which(idx_phi==c),]$Rhat
}
colnames(phiMatrix) = taxa_array
colnames(phiMatrix_rhat) = taxa_array

library(openxlsx)
ts = as.numeric(sub(".*\\_", "",sub("\\.RData.*", "",fileList[indexPick])))
wb <- createWorkbook()
addWorksheet(wb, "mean")
addWorksheet(wb, "Rhat")
writeData(wb, "mean", as.data.frame(phiMatrix), startRow = 1, startCol = 1)
writeData(wb, "Rhat", as.data.frame(phiMatrix_rhat), startRow = 1, startCol = 1)
saveWorkbook(wb, file = paste0(pathModelOutput,"/PHI_",ts,".xlsx"), overwrite = TRUE)












