# Setup
setwd(getwd())
library(zoo)
library(Rcpp)
library(lubridate)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(prodlim)
library(matrixStats)
library(dplyr)
library(plyr)
library(readxl)
library(writexl)
library(text.alignment)
library(vegan)
library(ggplot2)
library(useful)
library(MatchIt)
library(knitr)
library(mixtools)
library(latex2exp)
library(png)
library('colorspace')
library(patchwork)
library(gridExtra)
library(cowplot)
library(ggpubr)

# MAX INDEX OF DAYS FOR EACH SUBJECT
# subject day
# 1        A 643
# 2        B 726
# 3        C 733
# 4        D 729
# 5        E 717
# 6        F 722
# 7        G 726
# 8        H 729
# 9        I 729
# 10       J 730
# 11       K 728
# 12       L 735

# MIN INDEX OF DAYS FOR EACH SUBJECT
#    subject day
# 1        A   2
# 2        B   1
# 3        C   1
# 4        D   1
# 5        E   1
# 6        F   2
# 7        G   5
# 8        H   2
# 9        I   3
# 10       J   2
# 11       K   2
# 12       L   3


data_yagahi_2 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-2')
data_yagahi_3 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-3')
data_yagahi_4 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-4')
data_yagahi_5 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-5')
data_yagahi_6 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-6')
sort(as.character(unlist(data_yagahi_5[3,])))

avg_abundances = data_yagahi_5[c(1,3),2:102]
taxonomy_above_1percent=as.character(unlist(avg_abundances[2,which(avg_abundances[1,]>1)])) #take the ones with mean abundance above 1%
taxonomy_all = as.character(unlist(avg_abundances[2,which(avg_abundances[1,]>=0)]))
taxonomy_all = setdiff(taxonomy_all,c('uncultured','UNCLASSIFIED','uncultured bacterium','gut metagenome'))

data_yagahi_supp = read_excel('YAGAHI_DATA/41396_2021_937_MOESM2_ESM_edited.xlsx')
data_use = data_yagahi_5 # pick the data on the family level
column_names = c('Sample',as.character(data_use[3,2:dim(data_use)[2]]))
data_use = data_use[4:dim(data_use)[1],]
colnames(data_use) = column_names
data_use = select(data_use, -c(UNCLASSIFIED, uncultured, `gut metagenome`, `uncultured bacterium`))
data_use_merged = merge(data_use, data_yagahi_supp, by='Sample')

Subjects_use = unique(na.omit(data_use_merged)$Subject)

data_use_merged_keep = data_use_merged

# initial seeding - check the mum's sample?

save_data_all <- data.frame(matrix(ncol = 4, nrow = 0))
save_data_all_smooth <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(save_data_all) <- c('subject', 'day', 'taxa','abundance') # day = 0 will be the mum
colnames(save_data_all_smooth) <- c('subject', 'day', 'taxa','abundance') # day = 0 will be the mum

save_data_all_milkandsolid <- data.frame(matrix(ncol = 3, nrow = 0))
save_data_all_solid <- data.frame(matrix(ncol = 3, nrow = 0))

# USE A DAY GRID OF 21 MONTHS (630 DAYS) FROM DAY 5 TO 635
day_grid_prediction = seq(0,max(na.omit(as.numeric(unlist(data_use_merged_keep$day)))),1)


for (subj in Subjects_use){
  
  data_use_merged_mum = data_use_merged_keep %>% filter(Subject==paste0('p',subj,'mo'))
  data_use_merged     = data_use_merged_keep %>% filter(Subject==subj)
  cut0 = which(colnames(data_use_merged)=='Sample')+1
  cut1 = which(colnames(data_use_merged)=='Subject')-1
  cut2 = which(colnames(data_use_merged)=='Acetate')
  cut3 = which(colnames(data_use_merged)=='pH')
  
  taxonomy_use              = colnames(data_use_merged[c(cut0:cut1,cut2:cut3)])
  data_use_merged_sub       = data_use_merged[c('Sample','day','month',taxonomy_use)]
  #####################################################################################################################
  
  data_use_merged_sub2     = data_use_merged[c('day',taxonomy_use)]
  taxonomy_use_plot        = setdiff(colnames(data_use_merged_sub2),c('day'))
  data_use_merged_sub2     = na.omit(data_use_merged_sub2)
  data_use_merged_sub2     = data_use_merged_sub2%>%  mutate_all(as.numeric) # make all values numeric
  data_use_merged_sub2$day = round(data_use_merged_sub2$day)
  #ignore other for now
  for(taxonomy_pick in taxonomy_use_plot){
    
    print(c(subj,taxonomy_pick))
    data_use_merged_pick = data_use_merged_sub2[c('day',taxonomy_pick)]
    if(dim(data_use_merged_mum)[1]>0){
      data_use_merged_mum_pick = data_use_merged_mum[c(taxonomy_pick)]
    }else{
      data_use_merged_mum_pick = NA
    }
    
    data_use_merged_pick_with_mum = rbind(c(day=0,data_use_merged_mum_pick),data_use_merged_pick)
    data_use_merged_pick_with_mum[which((data_use_merged_pick_with_mum[,2])==0),2]=NA
    data_use_merged_pick_with_mum_complete     = matrix(NA,ncol=1,nrow=length(day_grid_prediction))
    data_use_merged_pick_with_mum_complete     = as.data.frame(data_use_merged_pick_with_mum_complete)
    data_use_merged_pick_with_mum_complete$day = day_grid_prediction
    data_use_merged_pick_with_mum_nona     = na.omit(data_use_merged_pick_with_mum)
    data_use_merged_pick_with_mum_complete = merge(data_use_merged_pick_with_mum_complete,data_use_merged_pick_with_mum_nona,by=c('day'),all=T)
    data_use_merged_pick_with_mum_complete = data_use_merged_pick_with_mum_complete[c('day',taxonomy_pick)]
    # 
    
    if(dim(data_use_merged_pick_with_mum_nona)[1]>15){ # smoooth for good amount of datapoints
    # ###### USING LOESS ##################################  
    # taxa_model  = loess(data_use_merged_pick_with_mum_nona[[taxonomy_pick]] ~ day, data_use_merged_pick_with_mum_nona, control = loess.control(surface = "direct"))
    taxa_model  = loess(data_use_merged_pick_with_mum_nona[[taxonomy_pick]] ~ day, data_use_merged_pick_with_mum_nona)
    smooth_data = predict(taxa_model, data.frame(day = day_grid_prediction), se = TRUE)
    
    # raw
    save_data_all_temp <- data.frame(matrix(ncol = 4, nrow = length(day_grid_prediction)))
    save_data_all_temp[,1]=subj
    save_data_all_temp[,2]=day_grid_prediction
    save_data_all_temp[,3]=taxonomy_pick
    save_data_all_temp[,4]=as.numeric(unlist(data_use_merged_pick_with_mum_complete[[taxonomy_pick]]))
    colnames(save_data_all_temp)=c('subject','day','taxa','abundance')
    
    # smooth
    save_data_all_temp_smooth <- data.frame(matrix(ncol = 4, nrow = length(day_grid_prediction)))
    save_data_all_temp_smooth[,1]=subj
    save_data_all_temp_smooth[,2]=day_grid_prediction
    save_data_all_temp_smooth[,3]=taxonomy_pick
    save_data_all_temp_smooth[,4]=as.numeric(unlist(smooth_data$fit))
    colnames(save_data_all_temp_smooth)=c('subject','day','taxa','abundance')
    
    ######################################################
    
    save_data_all = rbind(save_data_all,save_data_all_temp)
    save_data_all_smooth = rbind(save_data_all_smooth,save_data_all_temp_smooth)
    }
    
  }
}
saveRDS(save_data_all,'YAGAHI_DATA/save_data_all_RAW.rds')
write_xlsx(save_data_all,"YAGAHI_DATA/save_data_all_RAW.xlsx")

saveRDS(save_data_all_smooth,'YAGAHI_DATA/save_data_all_SMOOTH.rds')
write_xlsx(save_data_all_smooth,"YAGAHI_DATA/save_data_all_SMOOTH.xlsx")

## sanity check - NA is not included in averaging
# save_data_all_1 = save_data_all %>% filter(day==1 & taxa=='Bifidobacteriaceae') 

save_data_all_agg = aggregate(abundance~taxa+day,save_data_all,FUN=mean)
taxonomy_use  = unique(save_data_all_agg$taxa)
keepAllSmooth = matrix(NA,nrow=length(day_grid_prediction),ncol=length(taxonomy_use))
counter = 1;
for(taxonomy_pick in taxonomy_use){
  save_data_all_agg_pick = save_data_all_agg %>% filter(taxa==taxonomy_pick)
  taxa_model  = loess(abundance ~ day, save_data_all_agg_pick)
  smooth_data = predict(taxa_model, data.frame(day = day_grid_prediction), se = TRUE)
  keepAllSmooth[,counter]=smooth_data$fit
  counter = counter + 1;
  
}
colnames(keepAllSmooth)=taxonomy_use
keepAllSmooth=as.data.frame(keepAllSmooth)

saveRDS(keepAllSmooth,'YAGAHI_DATA/save_data_all_MEAN_SMOOTH.rds')
saveRDS(day_grid_prediction,'YAGAHI_DATA/day_grid_prediction.rds')
