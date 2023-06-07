# Setup
rm(list=ls())
# setwd(getwd())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
source("SETUP.R")
source("PREPARE_MILK.R")
useTotalAbundance=0 #if zero, relative abundance is returned
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
#### totalAbundance_df : increase in total abundance over time - this is needed to be used as the denominator of the interaction parameters

# abundanceArray_allSubjects
# abundanceArray_allSubjects_mum
# saved_data
# saved_data_milkandsolid
# saved_data_solid
# abundanceArray_meanSubjects



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


# Maybe insert y0 at time point 0?
y0_allSubjectsBoth  = cbind(uncoated_y0+coated_y0)
y0_meanSubjectsBoth = colMeans(y0_allSubjectsBoth)
abundanceArray_meanSubjects = rbind(y0_meanSubjectsBoth,abundanceArray_meanSubjects)

taxa_array     = unique(saved_data_infant$taxa)


pathModelOutput='/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/07062023/RDATA';
fileNamePick  = paste0(pathModelOutput,"/MODEL_D0_1686112768.RData")
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
n_chains <- fit@sim$chains
n_warmup <- fit@sim$warmup2[1]
n_iter   <-fit@sim$iter[[1]]
print(c(n_iter,n_warmup,n_chains))

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
# coated uncoated labelling
# taxa_array = c(paste0(taxa_array,"_uncoated"), paste0(taxa_array,"_coated"))
# summaryTable_use = summaryTable_use %>% rowwise() %>% mutate(taxa = taxa_array[as.numeric(sub(",.*","",sub(".*\\_", "", taxa)))])
summaryTable_use = summaryTable_use %>% rowwise() %>% mutate(taxa_idx = as.numeric(sub(",.*","",sub(".*\\_", "", taxa))))
summaryTable_use = summaryTable_use %>% rowwise() %>% mutate(taxa_idx_name = (taxa_idx-1) %% length(taxa_array) + 1)
summaryTable_use = summaryTable_use %>% rowwise() %>% mutate(taxa = taxa_array[taxa_idx_name])
summaryTable_use_agg_50 = aggregate(`50%` ~ taxa + t,summaryTable_use, FUN=sum)

summaryTable_use_keep = summaryTable_use
summaryTable_use = summaryTable_use_agg_50

abundanceArray_meanSubjects$day = as.numeric(rownames(abundanceArray_meanSubjects))
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects %>% pivot_longer(!day,names_to = "taxa", values_to = "abundance")
 
abundanceArray_meanSubjects_longer = abundanceArray_meanSubjects_longer %>% filter(day<max(summaryTable_use$t))
# summaryTable_use = summaryTable_use %>% filter(t<10)
graphics.off()
ggplot() +
  geom_point(data=abundanceArray_meanSubjects_longer,aes(x=day,y=abundance,fill=taxa),shape=21,size=1,colour = "black", fill = "white") +
  # geom_ribbon(data=summaryTable_use,aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=taxa),alpha=.5) +
  geom_line(data=summaryTable_use,aes(x=t,y=`50%`),colour="blue") +
  facet_wrap(~ taxa ,scales="free",nrow=2) +
  scale_colour_manual(values=c("grey20","grey"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  # scale_y_continuous(trans = 'log10')+
  labs(x="sample index",y="abundance")

# graphics.off()
# ggplot(summaryTable_use, aes(x = t, y = mean, color = taxa)) + geom_point()

# Save the interaction matricies to compare
# abundanceArray_meanSubjects



