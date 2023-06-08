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

pathModelOutput='/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/OUT/07062023/RDATA';
fileNamePick  = paste0(pathModelOutput,"/MODEL_D0_1686172623.RData")
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
# taxa_array = c(paste0(taxa_array,"_uncoated"),paste0(taxa_array,"_coated"))
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
abundanceArray_meanSubjects



