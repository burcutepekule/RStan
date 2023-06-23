rm(list=ls())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
library("bayesplot")
library("ggplot2")
library("rstanarm") 
library('rstan')
library('brms')
library('dplyr')
library(readxl)
library(writexl)
library(deSolve)

# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html

source("SETUP.R")
source('RESHAPE_DATA_KNIGHT_L6.R')
taxa_array = colnames(abundanceArray_meanSubjects)
abundanceArray_meanSubjects_use = abundanceArray_meanSubjects
abundanceArray_meanSubjects_use$day = rownames(abundanceArray_meanSubjects_use)
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects_use %>% pivot_longer(!day,names_to = 'taxa', values_to = 'abundance')

pathModelOutput='/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/23062023/RDATA';
fileList   = list.files(path =pathModelOutput, pattern='.xlsx')
timestamps = unique(sub("\\.xlsx.*","",sub(".*\\_", "",fileList)))  # Extract characters after pattern

##### Pick the file
indexPick          = 2
timestampPick      = timestamps[indexPick]
growthRateVector   = read_excel(paste0(pathModelOutput,"/GROWTH_",timestampPick,".xlsx"))
interactionMatrix  = read_excel(paste0(pathModelOutput,"/INTERACTIONS_",timestampPick,".xlsx"))
numTaxa            = length(taxa_array)
t_end              = max(as.numeric(unique(abundanceArray_meanSubjects_longer$day)))


#################################################

### load the model object - clean up afterwards
fileNamePick = paste0(pathModelOutput,'/MODEL_SMOOTH_0_ODESOLVER_0_B1B_',timestampPick,'.RData')
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

#################################################
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
summaryTable_use$t = days_array[as.numeric(summaryTable_use$t)]
summaryTable_use$mean = as.numeric(summaryTable_use$mean)
summaryTable_use$`97.5%` = as.numeric(summaryTable_use$`97.5%`)
summaryTable_use$`2.5%` = as.numeric(summaryTable_use$`2.5%`)
summaryTable_use$`50%` = as.numeric(summaryTable_use$`50%`)
summaryTable_use$taxa = paste0('y_',summaryTable_use$taxa)
summaryTable_use = summaryTable_use %>% rowwise() %>% mutate(taxa = taxa_array[as.numeric(sub(",.*","",sub(".*\\_", "", taxa)))] )

abundanceArray_meanSubjects_longer$type ='stan'

source("B1B_ODE_MODEL.R")


graphics.off()
ggplot() +
  geom_point(data=out_longer,aes(x=day,y=abundance,fill=taxa),shape=21,size=1,colour = "red", fill = "white") +
  geom_ribbon(data=summaryTable_use,aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=taxa),alpha=.5) +
  geom_line(data=summaryTable_use,aes(x=t,y=`50%`),colour="black") +
  facet_wrap(~ taxa ,scales="free",nrow=2) +
  scale_colour_manual(values=c("grey20","grey")) +
  scale_alpha_manual(values=c(1,0)) +
  # scale_y_continuous(trans = 'log10')+
  labs(x="sample index",y="abundance")

