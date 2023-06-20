#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# setwd(getwd())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
source("SETUP.R")
source('RESHAPE_DATA_KNIGHT.R')

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>0) {
  # default output file
  df_config = read.table(args[1], header=TRUE)
}


# df_config     = read.table('config_raw.txt', header = TRUE, sep = "", dec = ".")

  
scale         = df_config[which(df_config[,1]=='scale'),2]
smoothed      = df_config[which(df_config[,1]=='smoothed'),2]
numWarmup     = df_config[which(df_config[,1]=='numWarmup'),2]
numIterations = df_config[which(df_config[,1]=='numIterations'),2]
numChains     = df_config[which(df_config[,1]=='numChains'),2]
valSeed       = df_config[which(df_config[,1]=='valSeed'),2]
odeSolver     = df_config[which(df_config[,1]=='odeSolver'),2]

# HUMAN MICROBIOME PROJECT DATA

# y0_meanSubjects
# abundanceArray_meanSubjects

# data_knight_M5_t, day_grid_M5
# data_knight_F5_t, day_grid_F5 

## scale up
# scale = 1E11 # TOO BIG
# scale = 0 # DOABLE


if(scale>0){
  y0_meanSubjects             = round(scale*y0_meanSubjects)
  names(y0_meanSubjects)      = colnames(abundanceArray_meanSubjects)
  days_array = as.numeric(rownames(abundanceArray_meanSubjects))
  taxa_array = colnames(abundanceArray_meanSubjects)
  abundanceArray_meanSubjects = as.data.frame(lapply(abundanceArray_meanSubjects, function(x) as.numeric(x))) #normalize and scale to absolute abundance
  abundanceArray_meanSubjects = as.data.frame(lapply(abundanceArray_meanSubjects, function(x) as.integer(round(scale*x)))) #normalize and scale to absolute abundance
  abundanceArray_meanSubjects_longer$abundance = as.integer(round(scale*abundanceArray_meanSubjects_longer$abundance)) #normalize and scale to absolute abundance
  phi_use_exp                 = max(round(log10(y0_meanSubjects)));
  phi_use                     = 10^(phi_use_exp+0) # SEEMS LIKE 1E6 IS THE ONE THAT WORKS, SO NEED TO GET THIS FROM THE ORDER OF MAG. OF THE DATA
  
}else{
  names(y0_meanSubjects)      = colnames(abundanceArray_meanSubjects)
  days_array = as.numeric(rownames(abundanceArray_meanSubjects))
  taxa_array = colnames(abundanceArray_meanSubjects)
  abundanceArray_meanSubjects = as.data.frame(lapply(abundanceArray_meanSubjects, function(x) as.numeric(x))) #normalize and scale to absolute abundance
  phi_use_exp                 = max((log10(y0_meanSubjects)));
  phi_use                     = 10^(phi_use_exp+0) # SEEMS LIKE 1E6 IS THE ONE THAT WORKS, SO NEED TO GET THIS FROM THE ORDER OF MAG. OF THE DATA
}

indexes = 1:9 # ALL
abundanceArray_meanSubjects = abundanceArray_meanSubjects[,indexes]
y0_meanSubjects             = y0_meanSubjects[indexes]
taxa_array                  = taxa_array[indexes]
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects_longer %>% filter(taxa %in% taxa_array)

# maybe smooth before use?
if(smoothed==1){
  abundanceArray_meanSubjects_keep = abundanceArray_meanSubjects
  for (c in 1:dim(abundanceArray_meanSubjects)[2]){
    tempcolumn  = abundanceArray_meanSubjects_keep[,c]
    dayvector   = days_array
    taxa_model  = loess(tempcolumn ~ days_array,  family = c("gaussian"), span=7)
    smooth_data = predict(taxa_model, data.frame(day = days_array), se = TRUE)
    # plot(smooth_data$fit)
    abundanceArray_meanSubjects[,c]=unlist(smooth_data$fit)
  }
}


data_list = list(
  
  numTaxa      = length(taxa_array),
  numTimeSteps = length(days_array),
  y0           = y0_meanSubjects,
  observations = tibble(abundanceArray_meanSubjects),
  
  p_mu           = c(0,2),
  p_a_intra      = c(0,2),
  p_a_inter      = c(0,2),
  p_phi          = 1/phi_use,
  
  t0        = 0, #starting time
  t_data    = days_array, #time bins of data
  ts_pred   = days_array #time bins of prediction (not doing prediction currently)
)

if(scale>0){ # NEGATIVE BINOMIAL LIKELIHOOD
  # RECOMPILE EACH TIME
  if(file.exists("MODELS/MODEL_B0_INT.rds")){
    file.remove("MODELS/MODEL_B0_INT.rds")
  }
  M_model  = stan_model("MODELS/MODEL_B0_INT.stan")
  T_model  = sampling(M_model,data = data_list,warmup=numWarmup,iter=numIterations,chains=numChains,seed=valSeed,init="random", refresh = 10)
}else{ # NORMAL LIKELIHOOD
  # RECOMPILE EACH TIME
  if(odeSolver==45){
    if(file.exists("MODELS/MODEL_B0_FLT_45.rds")){
      file.remove("MODELS/MODEL_B0_FLT_45.rds")
    }
    print('Picked ODE_45 model')
    M_model  = stan_model("MODELS/MODEL_B0_FLT_45.stan")
    print('Compiled the model')
  }else{
    odeSolver = 0 #WILL GIVE THIS AS INPUT BUT STILL TO MAKE SURE IT IS A NUMBER FOR THE FILENAME
    if(file.exists("MODELS/MODEL_B0_FLT_BDF.rds")){
      file.remove("MODELS/MODEL_B0_FLT_BDF.rds")
    }
    print('Picked ODE_BDF model')
    M_model  = stan_model("MODELS/MODEL_B0_FLT_BDF.stan")
    print('Compiled the model')
  }
  
  T_model  = sampling(M_model,data = data_list,warmup=numWarmup,iter=numIterations,chains=numChains,seed=valSeed,init="random", refresh = 10)
}

todaystr    = format(Sys.Date(), "%d%m%Y");
tstamp      = as.numeric(Sys.time());
direc2save  = paste0("OUT/",todaystr,"/RDATA")
mkdir(direc2save)
save(T_model, file = paste0(direc2save,"/MODEL_SMOOTH_",smoothed,"_ODESOLVER_",odeSolver,"_",round(tstamp),".RData"))



