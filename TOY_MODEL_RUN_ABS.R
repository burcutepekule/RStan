# Setup
rm(list=ls())
# setwd(getwd())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
source("SETUP.R")

t_end=30
source('TOY_MODEL_ODE_ABS.R')

# HUMAN MICROBIOME PROJECT DATA

# y0_meanSubjects
# abundanceArray_meanSubjects

abundanceArray_meanSubjects = round(out_use[2:dim(out_use)[1],]) #make it integer to use neg_binomial_2_lpmf likelihood
# abundanceArray_meanSubjects = out_use[2:dim(out_use)[1],] #can use double for normal_lpdf likelihood
y0_meanSubjects             = out_init
phi_use                     = 1E6
days_array                  = out[2:dim(out)[1],]$time
taxa_array                  = colnames(out_use)


days_array                  = out[2:dim(out)[1],]$time
taxa_array                  = colnames(out_use)

graphics.off()
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects
abundanceArray_meanSubjects_longer$day = as.numeric(rownames(abundanceArray_meanSubjects_longer))-1
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects_longer %>% pivot_longer( !day, names_to = "taxa", values_to = "abundance")
abundanceArray_meanSubjects_longer$taxa=as.factor(abundanceArray_meanSubjects_longer$taxa)
abundanceArray_meanSubjects_longer$day=as.factor(abundanceArray_meanSubjects_longer$day)
abundanceArray_meanSubjects_longer$abundance=as.numeric(abundanceArray_meanSubjects_longer$abundance)
library("ggplot2")
ggplot(abundanceArray_meanSubjects_longer, aes(x = day, y = abundance, color = taxa)) + geom_point() 
dim(abundanceArray_meanSubjects_longer)

#cut short for debugging
# abundanceArray_meanSubjects = abundanceArray_meanSubjects[1:5,]
# days_array = as.numeric(rownames(abundanceArray_meanSubjects))+1

library(hablar)
abundanceArray_meanSubjects = tibble(abundanceArray_meanSubjects)


data_list = list(
  
  numTaxa      = length(taxa_array),
  numTimeSteps = length(days_array),
  y0 = y0_meanSubjects,
  observations = abundanceArray_meanSubjects,

  p_mu           = c(0.15,0.01),
  p_a_intra      = c(1,0.01),
  p_a_inter      = c(0.01,0.1),
  p_phi          = 1/phi_use,
  
  t0        = 0, #starting time
  t_data    = days_array, #time bins of data
  ts_pred   = days_array #time bins of prediction (not doing prediction currently)
)

# RECOMPILE EACH TIME
if(file.exists("MODELS/MODEL_TOY_ABS.rds")){
  file.remove("MODELS/MODEL_TOY_ABS.rds")
}

sinking = 0
wu = 100
ch = 500

M_model  = stan_model("MODELS/MODEL_TOY_ABS.stan")
T_model  = sampling(M_model,data = data_list,warmup=wu,iter=ch,seed=96,chains=4,init="random")

# if (sinking==0){
#   M_model  = stan_model("MODELS/MODEL_TOY_2.stan")
#   T_model  = sampling(M_model,data = data_list,warmup=wu,iter=ch,chains=1,init="random")
# }else{
#   sink("sink-examp.txt")
#   M_model  = stan_model("MODELS/MODEL_TOY_2.stan")
#   T_model  = sampling(M_model,data = data_list,warmup=wu,iter=ch,chains=1,init="random")
#   sink()
#   closeAllConnections()
# }

compartment_names = 'output_pred'
summaryTable = as.data.frame(summary(T_model,compartment_names)[[1]])
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


graphics.off()
ggplot(summaryTable_use, aes(x = t, y = mean, color = taxa)) + geom_point()

graphics.off()
ggplot() +
  geom_point(data=abundanceArray_meanSubjects_longer,aes(x=day,y=abundance,fill=taxa),shape=21,size=1,colour = "black", fill = "white") +
  geom_ribbon(data=summaryTable_use,aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=taxa),alpha=.5) +
  geom_line(data=summaryTable_use,aes(x=t,y=`50%`),colour="black") +
  facet_wrap(~ taxa ,scales="free",nrow=2) +
  scale_colour_manual(values=c("grey20","grey"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  # scale_y_continuous(trans = 'log10')+
  labs(x="sample index",y="abundance")


compartment_names = c('a_11','a_22','a_33')
summaryTable = as.data.frame(summary(T_model,compartment_names)[[1]])




