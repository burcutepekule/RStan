# Setup
rm(list=ls())
# setwd(getwd())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
source("SETUP.R")
source("PREPARE_MILK.R")
source("RESHAPE_DATA_YAGAHI.R")
##### OUTPUTS:
#### saved_data : df of all data
#### saved_data_milkandsolid : df keeping the day of transition to mixed feeding
#### saved_data_solid : df keeping the day of transition to solid feeding
#### abundanceArray_allSubjects : 3D array of abundances per subject per taxa per time point
#### abundanceArray_allSubjects_mum : 3D array of abundances per subject's mum at sampled time point
#### abundanceArray_allSubjects[days,subject,taxa]
#### abundanceArray_allSubjects_mum[days,subject,taxa]
#### taxa_array : array of taxa used
#### days_array : array of time points (days) used
#### subjects_array : array of subjects used
#### abundanceArray_meanSubjects : 2D array of abundances per taxa per time point (avgd over subjects)

abundanceArray_allSubjects
abundanceArray_allSubjects_mum
saved_data
saved_data_milkandsolid
saved_data_solid
abundanceArray_meanSubjects

ss_coating  = 0.36
coated_y0   = ss_coating*abundanceArray_allSubjects_mum
uncoated_y0 = (1-ss_coating)*abundanceArray_allSubjects_mum
abundanceArray_allSubjects_mum_reshaped = array_reshape(abundanceArray_allSubjects_mum, c(length(subjects_array), length(taxa_array)))
colnames(abundanceArray_allSubjects_mum_reshaped) = taxa_array
rownames(abundanceArray_allSubjects_mum_reshaped) = subjects_array
uncoated_y0 = (1-ss_coating)*abundanceArray_allSubjects_mum_reshaped
coated_y0   = ss_coating*abundanceArray_allSubjects_mum_reshaped
y0_allSubjects = cbind(uncoated_y0,coated_y0)
y0_meanSubjects= colMeans(y0_allSubjects)

growthRate_vector_in    = c(2,2,2,2,2,2,2,2,2,2)
interactionMat_in       = matrix(rexp(length(taxa_array)*length(taxa_array), rate=.1), nrow=length(taxa_array), ncol=length(taxa_array))
diag(interactionMat_in) = -1
interactionMat_vector_in= as.vector((interactionMat_in))
#### sanity check
# interactionMat_check = matrix(interactionMat_vector_in, nrow = length(taxa_array), byrow = length(taxa_array)) #convert array to matrix
# print(sum(interactionMat_in-interactionMat_check))
# ##[1] 0

data_list = list(
  
  numTaxa      = length(taxa_array),
  numTimeSteps = length(days_array),
  t_mixed      = round(mean(unlist(saved_data_milkandsolid$day))),
  t_solid      = round(mean(unlist(saved_data_solid$day))),
  binary_breastmilk     = 1, # all breastfed
  O2Dependency_vector   = c(0,1,0,1,0,1,1,0,0,0),
  HMODependency_vector  = c(1,0,0,0,1,0,0,0,0,0),
  solidDependency_vector= c(0,1,1,1,1,1,1,1,1,1),
  growthRate_vector     = growthRate_vector_in,
  interactionMat_vector = interactionMat_vector_in,
  
  y0 = y0_meanSubjects,
  observations = abundanceArray_meanSubjects,
  p_coating_vector = c(2,5),
  
  t0        = 0, #starting time
  t0_data   = days_array[1], #index of first sample
  t_sim_end = max(days_array), #total simulation time
  t_data    = days_array, #time bins of data
  ts_pred   = days_array #time bins of prediction (not doing prediction currently)
)

# RECOMPILE EACH TIME
if(file.exists("MODELS/MODEL_A2.rds")){
  file.remove("MODELS/MODEL_A2.rds")
}
M_model  = stan_model("MODELS/MODEL_A2.stan")
T_model  = sampling(M_model,data = data_list,warmup=50,iter=150,chains=1,init="random")

compartment_names = 'y'
summaryTable = as.data.frame(summary(T_model,compartment_names)[[1]])
summaryTable$populationNames = rownames(summaryTable)
summaryTable$t    = sub(",.*","",sub(".*y\\[", "", summaryTable$populationNames))  # Extract characters after pattern
summaryTable$taxa = sub("\\].*","",sub(".*,", "", summaryTable$populationNames))  # Extract characters after pattern
summaryTable_use = summaryTable[c('t','taxa','mean')]

for(tx in 1:length(taxa_array)){
  summaryTable_use_tx_uncoated = summaryTable_use %>% filter(taxa==tx)
  summaryTable_use_tx_coated   = summaryTable_use %>% filter(taxa==tx+length(taxa_array))
  summaryTable_use_tx_total    = summaryTable_use_tx_uncoated$mean+summaryTable_use_tx_coated$mean
  plot(summaryTable_use_tx_total)
  
  
}






