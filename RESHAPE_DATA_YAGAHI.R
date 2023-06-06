# run PREPARE_DATA_YAGAHI.R first if save_data_all.rds does not exist.
if(!file.exists("YAGAHI_DATA/save_data_all.rds")){
  source('PREPARE_DATA_YAGAHI.R')
}
saved_data              = readRDS('YAGAHI_DATA/save_data_all.rds')
saved_data_milkandsolid = readRDS('YAGAHI_DATA/save_data_all_milkandsolid.rds')
saved_data_solid        = readRDS('YAGAHI_DATA/save_data_all_solid.rds')

# REMOVE SUBJECT K -> DON'T KNOW WHEN TRANSITION TO SOLID FOOD HAPPENED
saved_data              = saved_data %>% filter(subject!='K')
saved_data_milkandsolid = saved_data_milkandsolid %>% filter(Subject!='K')
saved_data_solid        = saved_data_solid %>% filter(Subject!='K')

saved_data_mum   = saved_data %>% filter(day==-1)
saved_data_infant= saved_data %>% filter(day>=0)
saved_data_infant$abundance[saved_data_infant$abundance<0] = 0 #fix negative values introduced by loess interpolation
saved_data_mum$abundance[saved_data_mum$abundance<0] = 0 #fix negative values introduced by loess interpolation
colnames(saved_data_infant)[4]='relative_abundance'
colnames(saved_data_mum)[4]='relative_abundance'
saved_data_infant$relative_abundance=as.numeric(saved_data_infant$relative_abundance)

saved_data_infant_bacteria_only = saved_data_infant %>% filter(!(taxa  %in% c('Valerate','Acetate','Butyrate','Formate','Isobutyrate','Lactate','Isovalerate','Succinate','Mitochondria','pH')))
saved_data_infant_bacteria_only_agg = aggregate(relative_abundance~subject+day,data=saved_data_infant_bacteria_only,FUN=sum)
# above gives all around 100- so it is the percentage! 

# CONVERT RELATIVE TO ABSOLUTE ABUNDANCE
# check https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050177
# got the rough values from figure 2

library(tidyverse)
library(broom)

data_time     = c(0,1,2,3,20,50,180) #days
data_abundance= c(10^3,10^4,10^6,10^8,10^9,10^10,10^10) #rDNA copies/g
df_abundance  = tibble(t = data_time, y = data_abundance)
fit_abundance = nls(y ~ SSlogis(t, Asym, xmid, scal), data = df_abundance)
coeffs_abundance = coef(fit_abundance)
Asym        = coeffs_abundance[[1]]
xmid        = coeffs_abundance[[2]]
scal        = coeffs_abundance[[3]]
t_abundance = seq(0,720)
totalAbundance    = Asym/(1+exp((xmid-t_abundance)/scal));
totalAbundance_df = as.data.frame(cbind(t_abundance,totalAbundance))
colnames(totalAbundance_df)=c('day','total_abundance')
# plot(t_abundance,totalAbundance)
# plot(data_time,data_abundance)

saved_data_infant = merge(totalAbundance_df,saved_data_infant,by='day')
saved_data_infant$total_abundance = as.numeric(saved_data_infant$total_abundance)
saved_data_infant$relative_abundance = as.numeric(saved_data_infant$relative_abundance)
if(useTotalAbundance==1){
  saved_data_infant = saved_data_infant %>% dplyr::rowwise() %>% dplyr::mutate(abundance=0.01*total_abundance*relative_abundance)#since it's percentage
}else{
  saved_data_infant = saved_data_infant %>% dplyr::rowwise() %>% dplyr::mutate(abundance=0.01*relative_abundance) #since it's percentage
}
saved_data_infant = saved_data_infant[c('subject','day','taxa','abundance')]

families = c('Bifidobacteriaceae','Enterobacteriaceae','Lachnospiraceae','Streptococcaceae',
                  'Bacteroidaceae','Enterococcaceae','Staphylococcaceae','Peptostreptococcaceae',
                  'Clostridiaceae','Ruminococcaceae')
saved_data_infant_keep = saved_data_infant
saved_data_infant = saved_data_infant %>% filter(taxa %in% families)

saved_data_infant_wide = saved_data_infant %>%
  pivot_wider(names_from = subject, values_from = abundance)

subjects_array = unique(saved_data_infant$subject)
days_array     = unique(saved_data_infant$day)
taxa_array     = unique(saved_data_infant$taxa)
arrays_keep    = c()
for(tx in taxa_array){
  saved_data_infant_wide_tx = saved_data_infant_wide %>% filter(taxa==tx)
  saved_data_infant_wide_tx = saved_data_infant_wide_tx[,3:dim(saved_data_infant_wide_tx)[2]]
  saved_data_infant_wide_tx = as.data.frame(sapply(saved_data_infant_wide_tx,as.numeric))
  row.names    = days_array
  column.names = colnames(saved_data_infant_wide_tx)
  matrix.names = tx
  array_temp = array(unlist(saved_data_infant_wide_tx), 
                     dim = c(dim(saved_data_infant_wide_tx)[1],dim(saved_data_infant_wide_tx)[2],1),
                     dimnames=list(row.names,column.names,matrix.names))
  arrays_keep= cbind(arrays_keep, array_temp)
}

row.names    = days_array
column.names = colnames(saved_data_infant_wide_tx)
matrix.names = taxa_array
abundanceArray_allSubjects = array(arrays_keep, dim=c(length(days_array),length(subjects_array),length(matrix.names)),
                dimnames=list(row.names,column.names,matrix.names))

################### TAKE AVG TO USE IN STAN AS STARTERS
saved_data_infant_use = saved_data_infant_keep
saved_data_infant_use = saved_data_infant_use %>% filter(taxa %in% families)
saved_data_infant_agg = aggregate(abundance~taxa+day,saved_data_infant_use,FUN=mean)
saved_data_infant_agg_wider = saved_data_infant_agg %>%
  pivot_wider(names_from = taxa, values_from = abundance)

abundanceArray_meanSubjects = saved_data_infant_agg_wider[,2:dim(saved_data_infant_agg_wider)[2]]

################### USE MUM'S MICROBES FOR INITIAL SEEDING?
saved_data_mum$day=0
saved_data_mum = merge(totalAbundance_df,saved_data_mum,by='day')
saved_data_mum$total_abundance = as.numeric(saved_data_mum$total_abundance)
saved_data_mum$relative_abundance = as.numeric(saved_data_mum$relative_abundance)
if(useTotalAbundance==1){
  saved_data_mum = saved_data_mum %>% dplyr::rowwise() %>% dplyr::mutate(abundance=0.01*total_abundance*relative_abundance)#since it's percentage
}else{
  saved_data_mum = saved_data_mum %>% dplyr::rowwise() %>% dplyr::mutate(abundance=0.01*relative_abundance) #since it's percentage
}
saved_data_mum = saved_data_mum[c('subject','day','taxa','abundance')]

families = c('Bifidobacteriaceae','Enterobacteriaceae','Lachnospiraceae','Streptococcaceae',
             'Bacteroidaceae','Enterococcaceae','Staphylococcaceae','Peptostreptococcaceae',
             'Clostridiaceae','Ruminococcaceae')
saved_data_mum_keep = saved_data_mum
saved_data_mum = saved_data_mum %>% filter(taxa %in% families)

saved_data_mum_wide = saved_data_mum %>%
  pivot_wider(names_from = subject, values_from = abundance)

subjects_array = unique(saved_data_mum$subject)
days_array_mum = unique(saved_data_mum$day)
taxa_array     = unique(saved_data_mum$taxa)
arrays_keep    = c()

for(tx in taxa_array){
  saved_data_mum_wide_tx = saved_data_mum_wide %>% filter(taxa==tx)
  saved_data_mum_wide_tx = saved_data_mum_wide_tx[,3:dim(saved_data_mum_wide_tx)[2]]
  
  # IF there is no information regarding the mum's abundances, fill that in with the mean of all others
  mean_nonna = mean(unlist(saved_data_mum_wide_tx[which(!is.na(saved_data_mum_wide_tx))]))
  saved_data_mum_wide_tx[which(is.na(saved_data_mum_wide_tx))] = mean_nonna
  saved_data_mum_wide_tx = t(as.data.frame(sapply(saved_data_mum_wide_tx,as.numeric)))
  row.names    = days_array_mum
  column.names = colnames(saved_data_mum_wide_tx)
  matrix.names = tx
  array_temp = array(unlist(saved_data_mum_wide_tx), 
                     dim = c(dim(saved_data_mum_wide_tx)[1],dim(saved_data_mum_wide_tx)[2],1),
                     dimnames=list(row.names,column.names,matrix.names))
  arrays_keep= cbind(arrays_keep, array_temp)
}

row.names    = days_array_mum
column.names = colnames(saved_data_mum_wide_tx)
matrix.names = taxa_array
abundanceArray_allSubjects_mum = array(arrays_keep, dim=c(length(days_array_mum),length(subjects_array),length(matrix.names)),
                                   dimnames=list(row.names,column.names,matrix.names))




# ##### sanity check - difference should be zero!
# all_sum=0
# for(tx in taxa_array){
#   saved_data_infant_wide_tx     = saved_data_infant_wide %>% filter(taxa==tx)
#   saved_data_infant_wide_tx_use = sapply(saved_data_infant_wide_tx[,3:dim(saved_data_infant_wide_tx)[2]],as.numeric)
#   all_sum = all_sum + sum(abundanceArray_allSubjects[,,tx] - saved_data_infant_wide_tx_use)
# }
# print(all_sum)
# ####[1] 0 -> AND INDEED IT IS ZERO!

