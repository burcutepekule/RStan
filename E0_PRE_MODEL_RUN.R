# Setup
rm(list=ls())
# setwd(getwd())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
source("SETUP.R")
source("PREPARE_MILK.R")
useTotalAbundance=0 #if zero, relative abundance is returned
source("RESHAPE_DATA_YAGAHI_RAW.R")
##### OUTPUTS:
# abundanceArray_meanSubjects
# time_grid_prediction

# read these from saved tables
# use 1685924896 in /Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/04062023/RDATA as demo
estimations_growth      = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/04062023/RDATA/GROWTH_1685924896.xlsx', sheet='mean')
estimations_interaction = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/04062023/RDATA/INTERACTIONS_1685924896.xlsx', sheet='mean')
growthRate_vector_in    = unlist(estimations_growth)
interactionMat_vector_in= as.vector(unlist(estimations_interaction))

#### sanity check
# interactionMat_check = matrix(interactionMat_vector_in, nrow = length(taxa_array), byrow = length(taxa_array)) #convert array to matrix
# print(sum(interactionMat_in-interactionMat_check))
# ##[1] 0

######### FITTING THIS PART ONLY MAKES SENSE FOR THE FIRST MONTH, WHERE THERE IS NO RESPONSE.
######### JUST TO CHECK WHETHER SUCH FITTING IS POSSIBLE

time_grid_prediction        = time_grid_prediction[time_grid_prediction<=30]
days_array                  = time_grid_prediction[2:length(time_grid_prediction)] 
abundanceArray_meanSubjects = abundanceArray_meanSubjects %>% filter(day %in% days_array)
abundanceArray_meanSubjects = abundanceArray_meanSubjects[families]
days_array_pred             = days_array

ss_coating      = 0.36 #adult coating ratio
coated_y0       = ss_coating*abundanceArray_meanSubjects[1,]
uncoated_y0     = (1-ss_coating)*abundanceArray_meanSubjects[1,]
y0_meanSubjects = unlist(cbind(uncoated_y0,coated_y0))
taxa_array      = families

index_check_1 = which(days_array==10) #coating ratio checkpoint
index_check_2 = which(days_array==30) #coating ratio checkpoint

phi_use = 1e-1 # for relative abundances

data_list = list(
  
  numTaxa      = length(taxa_array),
  numTimeSteps = length(days_array),
  t_mixed      = round(mean(unlist(saved_data_milkandsolid$day))),
  t_solid      = round(mean(unlist(saved_data_solid$day))),
  icheck_1     = index_check_1,
  icheck_2     = index_check_2,
  binary_breastmilk     = 1, # all breastfed
  O2Dependency_vector   = c(1,1,0,0,0,0,0,1,1),
  HMODependency_vector  = c(0,1,1,1,0,0,0,0,0),
  solidDependency_vector= c(1,0,0,1,1,1,1,1,1),
  TLR9_vector           = c(),
  TLR4_vector           = c(),
  growthRate_vector     = growthRate_vector_in,
  interactionMat_vector = interactionMat_vector_in,
  IgA_halflife          = 1/7, # one week?


  y0           = y0_meanSubjects,
  observations = abundanceArray_meanSubjects,
  
  p_coating_vector    = c(2,2), # beta distribution
  p_O2_vector         = c(2,2), # beta distribution
  p_HMO_vector        = c(2,2), # beta distribution
  p_solid_vector      = c(2,2), # beta distribution
  p_bias_inflammation = c(2,2), # beta distribution
  p_phi               = 1/phi_use,
  
  
  t0        = 0, #starting time
  t0_data   = days_array[1], #index of first sample
  t_sim_end = max(days_array), #total simulation time
  t_data    = days_array, #time bins of data
  ts_pred   = days_array_pred  #time bins of prediction (not doing prediction currently)
)

# RECOMPILE EACH TIME
if(file.exists("MODELS/MODEL_EO_PRE.rds")){
  file.remove("MODELS/MODEL_E0_PRE.rds")
}
M_model  = stan_model("MODELS/MODEL_E0_PRE.stan")

# sink('C0_RUN.txt')

T_model     = sampling(M_model,data = data_list,warmup=250,iter=1000,chains=8,init="random",cores=mc.cores,refresh=10)
todaystr    = format(Sys.Date(), "%d%m%Y");
tstamp      = as.numeric(Sys.time());
direc2save  = paste0("OUT/",todaystr,"/RDATA")
mkdir(direc2save)
save(T_model, file = paste0(direc2save,"/MODEL_E0_PRE_",round(tstamp),".RData"))

# sink()
# closeAllConnections()

# compartment_names = 'y'
# summaryTable = as.data.frame(summary(T_model,compartment_names)[[1]])
# summaryTable$populationNames = rownames(summaryTable)
# summaryTable$t    = sub(",.*","",sub(".*y\\[", "", summaryTable$populationNames))  # Extract characters after pattern
# summaryTable$taxa = sub("\\].*","",sub(".*,", "", summaryTable$populationNames))  # Extract characters after pattern
# summaryTable_use = summaryTable[c('t','taxa','mean')]
# 
# for(tx in 1:length(taxa_array)){
#   summaryTable_use_tx_uncoated = summaryTable_use %>% filter(taxa==tx)
#   summaryTable_use_tx_coated   = summaryTable_use %>% filter(taxa==tx+length(taxa_array))
#   summaryTable_use_tx_total    = summaryTable_use_tx_uncoated$mean+summaryTable_use_tx_coated$mean
#   plot(summaryTable_use_tx_total)
# }


