# Setup
rm(list=ls())
setwd(getwd())
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


growthRate_vector_in    = c(2,2,2,2,2,2,2,2,2,2)
interactionMat_in       = matrix(rexp(length(taxa_array)*length(taxa_array), rate=.1), nrow=length(taxa_array), ncol=length(taxa_array))
diag(interactionMat_in) = -1
interactionMat_vector_in= as.vector(t(interactionMat_in))
#### sanity check
# interactionMat_check = matrix(interactionMat_vector_in, nrow = length(taxa_array), byrow = length(taxa_array)) #convert array to matrix
# print(sum(interactionMat_in-interactionMat_check))
# ##[1] 0

data_list = list(
  
  numSubjects  = length(subjects_array),
  numTaxa      = length(taxa_array),
  numTimeSteps = length(days_array),
  t_mixed      = mean(unlist(saved_data_milkandsolid$day)),
  t_solid      = mean(unlist(saved_data_solid$day)),
  binary_breastmilk     = 1, # all breastfed
  O2Dependency_vector   = c(0,1,0,1,0,1,1,0,0,0),
  HMODependency_vector  = c(1,0,0,0,1,0,0,0,0,0),
  growthRate_vector     = growthRate_vector_in,
  interactionMat_vector = interactionMat_vector_in,
  
  y0 = y0_allSubjects,
  observations = abundanceArray_allSubjects,
  p_coating_vector = c(2,5),
  
  t0      = 0, #starting time
  t0_data = days_array[1], #index of first sample
  t_total = seq(0,max(days_array)), #total simulation time
  t_data  = days_array, #time bins of data
  ts_pred = days_array #time bins of prediction (not doing prediction currently)
)

M_model  = stan_model("MODELS/MODEL_A0.stan")

T_model  = sampling(M_model,data = data_list,warmup=50,iter=150,chains=2,init="random")
save(T_model, file =paste0("RDATA/T_model_NASAL_all.RData"))

parameterSummary = summary(T_model, c("mu_a","mu_c","mu_d","alpha_aa","alpha_cc","alpha_dd",
                                      "alpha_ac","alpha_ca","alpha_ad","alpha_da","alpha_dc","alpha_cd",
                                      "phi_a"), probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
write.table(parameterSummary$summary, file = paste0("CSVS/",prename_parameters,"NASAL_all.csv"))
xtable(parameterSummary$summary)


pp = c("comp_SA_0","comp_SA_89","comp_CP_89","comp_SA_34","comp_DP_34","comp_SA_M","comp_CP_M","comp_DP_M")
model_output_all = summary(T_model,pp)[[1]] %>%
  tibble::as_tibble() %>%
  mutate(t=rep(1:data_list$S,8),
         date=1+t-1,
         eta="-100%",
         type=rep(pp,each=data_list$S)) %>%
  mutate(populations=factor(type,levels=pp,
                            labels=c("S.aureus (alone)",
                                     "S.aureus (with C.pseu)", "C.pseudodiphtheriticum (with S. aureus)",
                                     "S.aureus (with D.pigrum)","D.pigrum (with S. aureus)",
                                     "S.aureus (in mix)","C.pseudodiphtheriticum (in mix)","D.pigrum (in mix)")))


numSamples       = numSamples_89+samplesAdd;
samples_sim      = model_output_all$date;
p1  = model_output_all$populations[1:numSamples];
p11 = model_output_all$populations[1];
p2  = model_output_all$populations[(numSamples+1):(2*numSamples)];
p22 = model_output_all$populations[(2*numSamples)];
p3  = model_output_all$populations[(2*numSamples+1):(3*numSamples)];
p33 = model_output_all$populations[(3*numSamples)];
p4  = model_output_all$populations[(3*numSamples+1):(4*numSamples)];
p44 = model_output_all$populations[(4*numSamples)];
p5  = model_output_all$populations[(4*numSamples+1):(5*numSamples)];
p55 = model_output_all$populations[(5*numSamples)];
p6  = model_output_all$populations[(5*numSamples+1):(6*numSamples)];
p66 = model_output_all$populations[(6*numSamples)];
p7  = model_output_all$populations[(6*numSamples+1):(7*numSamples)];
p77 = model_output_all$populations[(7*numSamples)];
p8  = model_output_all$populations[(7*numSamples+1):(8*numSamples)];
p88 = model_output_all$populations[(8*numSamples)];

observed_data_0    = tibble(c(samples_sim[1:numSamples_0]),
                          c(SA_data_0),unlist(list(p1[1:numSamples_0])),
                          .name_repair = ~ c("date", "popsize","populations"))

observed_data_89   = tibble(c(samples_sim[1:numSamples_89]),
                            c(SA_data_89),unlist(list(p2[1:numSamples_89])),
                            .name_repair = ~ c("date", "popsize","populations"))

observed_data_34   = tibble(c(samples_sim[1:numSamples_34]),
                          c(SA_data_34),unlist(list(p4[1:numSamples_34])),
                          .name_repair = ~ c("date", "popsize","populations"))

observed_data_M    = tibble(c(samples_sim[1:numSamples_mix]),
                            c(SA_data_mix),unlist(list(p6[1:numSamples_mix])),
                            .name_repair = ~ c("date", "popsize","populations"))

pp_plot = c("comp_SA_0","comp_SA_89","comp_SA_34","comp_SA_M")
model_output_plot = summary(T_model,pp_plot)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$S,4),
         date=1+t-1,
         eta="-100%",
         type=rep(pp_plot,each=data_list$S)) %>%
  mutate(populations=factor(type,levels=pp_plot,
                            labels=c("S.aureus (alone)","S.aureus (with C.pseu)",
                                     "S.aureus (with D.pigrum)","S.aureus (in mix)")))

ggplot() +
  geom_point(data=observed_data_0,aes(x=date,y=popsize,fill=populations),shape=21,size=1,colour = "black", fill = "white") +
  geom_point(data=observed_data_89,aes(x=date,y=popsize,fill=populations),shape=21,size=1,colour = "black", fill = "white") +
  geom_point(data=observed_data_34,aes(x=date,y=popsize,fill=populations),shape=21,size=1,colour = "black", fill = "white") +
  geom_point(data=observed_data_M,aes(x=date,y=popsize,fill=populations),shape=21,size=1,colour = "black", fill = "white") +
  geom_ribbon(data=model_output_plot,aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=populations),alpha=.5) +
  geom_line(data=model_output_plot,aes(x=date,y=`50%`),colour="black") +
  facet_wrap(~ populations ,scales="free",nrow=2) +
  scale_colour_manual(values=c("grey20","grey"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  # scale_y_continuous(trans = 'log10')+
  labs(x="sample index (30 mins)",y="intensity")

ggsave(paste0("FIGS/",prename_figures,"NASAL_all.png"),width = 15, height = 7)

parameterSummary$summary[,1]

