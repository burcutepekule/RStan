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

# USE A DAY GRID OF 21 MONTHS (630 DAYS) FROM DAY 5 TO 635
day_grid      = seq(5,634,1)

data_yagahi_2 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-2')
data_yagahi_3 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-3')
data_yagahi_4 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-4')
data_yagahi_5 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-5')
data_yagahi_6 = read_excel('YAGAHI_DATA/Tsukuda_Yagahi(2021)ISME_16S-data.xlsx', sheet='level-6')

data_use = data_yagahi_6 # pick the data on the genus level
sort(as.character(unlist(data_use[3,])))

avg_abundances = data_use[c(1,3),2:dim(data_use)[2]]
taxonomy_above_1percent=as.character(unlist(avg_abundances[2,which(avg_abundances[1,]>1)])) #take the ones with mean abundance above 1%
taxonomy_all = as.character(unlist(avg_abundances[2,which(avg_abundances[1,]>=0)]))
taxonomy_all = setdiff(taxonomy_all,c('uncultured','UNCLASSIFIED','uncultured bacterium','gut metagenome','unidentified',
                                      "uncultured Peptostreptococcus sp.","uncultured Porphyromonadaceae bacterium","Finegoldia magna"))

data_yagahi_supp = read_excel('YAGAHI_DATA/41396_2021_937_MOESM2_ESM_edited.xlsx')
column_names = c('Sample',as.character(data_use[3,2:dim(data_use)[2]]))
data_use = data_use[4:dim(data_use)[1],]
colnames(data_use) = column_names
data_use = select(data_use, -c(`Finegoldia magna`, UNCLASSIFIED, uncultured, `gut metagenome`, unidentified, `uncultured bacterium`,`uncultured Peptostreptococcus sp.`,`uncultured Porphyromonadaceae bacterium`))
data_use_merged = merge(data_use, data_yagahi_supp, by='Sample')

Subjects_use = unique(na.omit(data_use_merged)$Subject)

data_use_merged_keep = data_use_merged

# initial seeding - check the mum's sample?

save_data_all <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(save_data_all) <- c('subject', 'day', 'taxa','abundance') # day = -1 will mean the mum

save_data_all_milkandsolid <- data.frame(matrix(ncol = 3, nrow = 0))
save_data_all_solid <- data.frame(matrix(ncol = 3, nrow = 0))

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
  
  taxonomy_use_plot        = setdiff(colnames(data_use_merged_sub),c('Sample','day','month'))
  data_use_merged_sub$days = 30*round(data_use_merged_sub$month)
  data_use_merged_sub = select(data_use_merged_sub, -c(Sample)) #drop Sample ID column
  data_use_merged_sub = na.omit(data_use_merged_sub)
  data_use_merged_sub = data_use_merged_sub%>%  mutate_all(as.numeric) # make all values numeric
  data_use_merged_sub$month=as.factor(round(data_use_merged_sub$month))
  data_use_merged_sub$days=as.factor(round(data_use_merged_sub$days))
  
  data_yagahi_supp_milkandsolid  = data_yagahi_supp %>% filter(feeding=='milk & solid')
  data_yagahi_supp_solid         = data_yagahi_supp %>% filter(feeding=='solid')
  
  data_yagahi_supp_milkandsolid_day  = aggregate(day ~ Subject, data_yagahi_supp_milkandsolid, FUN=min)
  data_yagahi_supp_solid_day         = aggregate(day ~ Subject, data_yagahi_supp_solid, FUN=min)
  
  data_yagahi_supp_milkandsolid_day = data_yagahi_supp_milkandsolid_day %>% rowwise() %>% mutate(month=day/30)
  data_yagahi_supp_solid_day = data_yagahi_supp_solid_day %>% rowwise() %>% mutate(month=day/30)
  
  
  data_yagahi_supp_milkandsolid_day = data_yagahi_supp_milkandsolid_day %>% filter(Subject==subj)
  data_yagahi_supp_solid_day        = data_yagahi_supp_solid_day %>% filter(Subject==subj)
  
  # data_use_merged_sub_check = data_use_merged_sub[,3:11]
  # unique(rowSums(data_use_merged_sub_check)) #ok so summation is done correctly - some 99.9% due to exclusion of "unclassified"
  
  #ignore other for now
  for(taxonomy_pick in taxonomy_use_plot){
    
    print(c(subj,taxonomy_pick))
    data_use_merged_pick = data_use_merged_sub[c('day','days','month',taxonomy_pick)]
    if(dim(data_use_merged_mum)[1]>0){
      data_use_merged_mum_pick = data_use_merged_mum[c(taxonomy_pick)]
    }else{
      data_use_merged_mum_pick = NA
    }

    # graphics.off()
    # png(file =paste0("YAGAHI_DATA/FIGS/FAMILY_MODEL_LEVEL_GENUS_",taxonomy_pick,"_",subj,".png"),   # The directory you want to save the file in
    #     width     = 8,
    #     height    = 3,
    #     units     = "in",
    #     res       = 300)
    # p=ggplot()+ geom_point(data=data_use_merged_pick,aes(x=day, y=.data[[taxonomy_pick]])) + theme_minimal()+
    #   geom_vline(data=data_yagahi_supp_milkandsolid_day, mapping=aes(xintercept=day), color="blue") +
    #   geom_vline(data=data_yagahi_supp_solid_day, mapping=aes(xintercept=day), color="blue") +
    #   geom_smooth(data=data_use_merged_pick,aes(x=day, y=.data[[taxonomy_pick]]),color = "black",se=TRUE, size = 0.5,  level=0.95, inherit.aes = T)+
    #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    # print(p)
    # dev.off()
    
    
    save_data_all_milkandsolid = rbind(save_data_all_milkandsolid,data_yagahi_supp_milkandsolid_day)
    save_data_all_solid = rbind(save_data_all_solid,data_yagahi_supp_solid_day)
    
    ###### USING GEOM_SMOOTH  ##################################  
    # smooth_data = ggplot_build(p)$data[[4]]
    # save_data_all_temp <- data.frame(matrix(ncol = 4, nrow = dim(smooth_data)[1])+1)
    # # first row -> day = -1 -> abundance of the mum microbes
    # save_data_all_temp[1,2]=-1
    # save_data_all_temp[1,4]=unlist(data_use_merged_mum_pick)
    # save_data_all_temp[2:(dim(smooth_data)[1]+1),2]=smooth_data$x
    # save_data_all_temp[2:(dim(smooth_data)[1]+1),4]=smooth_data$y
    #############################################################  
    
    ###### USING LOESS ##################################  
    taxa_model  = loess(data_use_merged_pick[[taxonomy_pick]] ~ day, data_use_merged_pick, control = loess.control(surface = "direct"))
    smooth_data = predict(taxa_model, data.frame(day = day_grid), se = TRUE)
    save_data_all_temp <- data.frame(matrix(ncol = 4, nrow = length(day_grid)+1))
    # first row -> day = -1 -> abundance of the mum microbes
    save_data_all_temp[1,2]=-1
    save_data_all_temp[1,4]=unlist(data_use_merged_mum_pick)
    save_data_all_temp[2:(length(day_grid)+1),2]=day_grid
    save_data_all_temp[2:(length(day_grid)+1),4]=as.numeric(unlist(smooth_data$fit))
    ######################################################
    
    save_data_all_temp[,1]=subj
    save_data_all_temp[,3]=taxonomy_pick
    colnames(save_data_all_temp)=c('subject','day','taxa','abundance')
    save_data_all = rbind(save_data_all,save_data_all_temp)
    

  }
}

save_data_all_milkandsolid=unique(save_data_all_milkandsolid)
save_data_all_solid=unique(save_data_all_solid)



# saveRDS(save_data_all,'YAGAHI_DATA/save_data_all.rds')
# saveRDS(save_data_all_milkandsolid,'YAGAHI_DATA/save_data_all_milkandsolid.rds')
# saveRDS(save_data_all_solid,'YAGAHI_DATA/save_data_all_solid.rds')
# 
# write_xlsx(save_data_all,"YAGAHI_DATA/save_data_all.xlsx")
# write_xlsx(save_data_all_milkandsolid,"YAGAHI_DATA/save_data_all_milkandsolid.xlsx")
# write_xlsx(save_data_all_solid,"YAGAHI_DATA/save_data_all_solid.xlsx")


