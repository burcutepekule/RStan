# Setup
rm(list=ls())
# setwd(getwd())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
source("SETUP.R")
source('RESHAPE_DATA_KNIGHT.R')

# HUMAN MICROBIOME PROJECT DATA

# y0_meanSubjects
# abundanceArray_meanSubjects

# data_knight_M5_t, day_grid_M5
# data_knight_F5_t, day_grid_F5 

data_knight_M5_abs = data_knight_M5
data_knight_F5_abs = data_knight_F5

numeric_cols  = names(data_knight_M5)[sapply(data_knight_M5, is.numeric)]
data_knight_M5_abs[, numeric_cols] <- lapply(data_knight_M5[, numeric_cols], function(x) ((10^11)*x))#scale to absolute abundance
numeric_cols  = names(data_knight_F5)[sapply(data_knight_F5, is.numeric)]
data_knight_F5_abs[, numeric_cols] <- lapply(data_knight_F5[, numeric_cols], function(x) ((10^11)*x))#scale to absolute abundance

abundanceArray_meanSubjects = data_knight_M5_t
y0_meanSubjects             = as.numeric(unlist(abundanceArray_meanSubjects[1,]))
y0_meanSubjects[which(y0_meanSubjects==0)] = abs(0.01*rnorm(length(which(y0_meanSubjects==0))))
y0_meanSubjects_abs = 10^11*y0_meanSubjects
names(y0_meanSubjects)=colnames(abundanceArray_meanSubjects)
names(y0_meanSubjects_abs)=colnames(y0_meanSubjects_abs)
days_array = day_grid_M5+1
taxa_array = families

graphics.off()
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects
abundanceArray_meanSubjects_longer$day = rownames(abundanceArray_meanSubjects_longer)
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects_longer %>% pivot_longer( !day, names_to = "taxa", values_to = "abundance")
abundanceArray_meanSubjects_longer$taxa=as.factor(abundanceArray_meanSubjects_longer$taxa)
abundanceArray_meanSubjects_longer$day=as.factor(abundanceArray_meanSubjects_longer$day)
abundanceArray_meanSubjects_longer$abundance=as.numeric(abundanceArray_meanSubjects_longer$abundance)
library("ggplot2")
ggplot(abundanceArray_meanSubjects_longer, aes(x = day, y = abundance, color = taxa)) + geom_point() 


#cut short for debugging
abundanceArray_meanSubjects = abundanceArray_meanSubjects[1:5,]
days_array = as.numeric(rownames(abundanceArray_meanSubjects))+1

data_list = list(
  
  numTaxa      = length(taxa_array),
  numTimeSteps = length(days_array),

  y0 = y0_meanSubjects,
  observations = abundanceArray_meanSubjects,
  
  p_growth_rates = c(0,0.01), # ~ beta
  p_interaction_terms_diag = c(0,0.01),# ~ beta
  p_interaction_terms_nondiag = c(0,0.01),# ~ normal
  p_phi       = 1/1e-1, # ~ exponential - THIS VALUE IS SUPER IMPORTANT! THIS CAN CAUSE LOG PROBABILITY TO FUCK UP AND NOT EVAL THE WHOLE THING! -> # SAMPLING FOR MODEL 'MODEL_B0' NOW (CHAIN 1).Chain 1: Unrecoverable error evaluating the log probability at the initial value. Chain 1: Exception: CVode failed with error flag -1.  (in 'model93061b4041c6_MODEL_B0' at line 123)
  t0        = 0, #starting time
  t0_data   = days_array[1], #index of first sample
  t_sim_end = max(days_array), #total simulation time
  t_data    = days_array, #time bins of data
  ts_pred   = days_array #time bins of prediction (not doing prediction currently)
)

# RECOMPILE EACH TIME
if(file.exists("MODELS/MODEL_B0.rds")){
  file.remove("MODELS/MODEL_B0.rds")
}

sinking = 0

if (sinking==0){
  M_model  = stan_model("MODELS/MODEL_B0.stan")
  T_model  = sampling(M_model,data = data_list,warmup=50,iter=150,chains=1,init="random")
}else{
  sink("sink-examp.txt")
  M_model  = stan_model("MODELS/MODEL_B0.stan")
  T_model  = sampling(M_model,data = data_list,warmup=50,iter=150,chains=1,init="random")
  sink()
}




