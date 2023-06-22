# Setup
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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


data_yagahi_2 = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/Tsukuda_Yagahi(2021)ISME_16S-data_edited.xlsx', sheet='level-2')
data_yagahi_3 = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/Tsukuda_Yagahi(2021)ISME_16S-data_edited.xlsx', sheet='level-3')
data_yagahi_4 = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/Tsukuda_Yagahi(2021)ISME_16S-data_edited.xlsx', sheet='level-4')
data_yagahi_5 = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/Tsukuda_Yagahi(2021)ISME_16S-data_edited.xlsx', sheet='level-5')
data_yagahi_6 = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/Tsukuda_Yagahi(2021)ISME_16S-data_edited.xlsx', sheet='level-6')


data_yagahi_supp = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/41396_2021_937_MOESM2_ESM_edited.xlsx')
data_use = data_yagahi_6 # pick the data on the GENUS level
column_names = c('Sample',as.character(data_use[3,2:dim(data_use)[2]]))
data_use = data_use[4:dim(data_use)[1],]
colnames(data_use) = column_names
data_use = select(data_use, -c(UNCLASSIFIED, uncultured, `gut metagenome`, `uncultured bacterium`,`uncultured Firmicutes bacterium`, `uncultured Porphyromonadaceae bacterium`,
                               `uncultured Peptostreptococcus sp.`,unidentified, `Finegoldia magna`))
data_use_merged = merge(data_use, data_yagahi_supp, by='Sample')

cut0  = which(colnames(data_use_merged)=='Sample')+1
cut1  = which(colnames(data_use_merged)=='Subject')-1
cut2  = which(colnames(data_use_merged)=='Acetate')
cut3  = which(colnames(data_use_merged)=='pH')


taxonomy_use   = colnames(data_use_merged[c(cut0:cut1,cut2:cut3)])

# # sanity check
# data_use_merged_sub_check = data_use_merged[c(taxonomy_use)]
# data_use_merged_sub_check = data_use_merged_sub_check[,1:243]
# data_use_merged_sub_check = as.data.frame(data_use_merged_sub_check)
# data_use_merged_sub_check = mutate_all(data_use_merged_sub_check, function(x) as.numeric(as.character(x)))
# max(rowSums(data_use_merged_sub_check)) # 100, phew

data_use_merged_sub = data_use_merged[c('day','month',taxonomy_use)]
data_use_merged     = mutate_all(data_use_merged, function(x) as.numeric(as.character(x)))
data_use_merged_sub$days = 30*round(data_use_merged_sub$month)
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

# taxonomy_pick=taxonomy_use[1]

taxanomy_table = c()
taxanomy_table$genus = c('Escherichia-Shigella','Staphylococcus','Streptococcus','Lactobacillus','Bifidobacterium',
                         'Bacteroides','Anaerostipes','Blautia',"[Eubacterium] hallii group",'Clostridium sensu stricto 1',
                         'Faecalibacterium','Ruminococcus')

taxanomy_table$family = c('Enterobacteriaceae','Staphylococcaceae','Streptococcaceae','Lactobacillaceae','Bifidobacteriaceae',
                         'Bacteroidaceae','Lachnospiraceae','Lachnospiraceae',"Lachnospiraceae",'Clostridiaceae',
                         'Ruminococcaceae','Ruminococcaceae') 

taxanomy_table = as.data.frame(taxanomy_table)
all_families = unique(taxanomy_table$family)

for(family_pick in all_families){
  print(family_pick)
  taxonomy_use = taxanomy_table %>% filter(family==family_pick)
  genus_use = taxonomy_use$genus
  data_use_merged_pick = data_use_merged_sub[c('day','days','month',genus_use)]
  data_use_merged_pick$family = rowSums(data_use_merged_pick[genus_use])
  colnames(data_use_merged_pick)[dim(data_use_merged_pick)[2]] = family_pick
  
  data_use_merged_pick_plot = data_use_merged_pick[c('day','days','month',family_pick)]
  
  graphics.off()
  png(file =paste0("/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/YAGAHI_DATA/FIGS/FAM_AGG_2USE_",family_pick,".png"),   # The directory you want to save the file in
      width     = 8,
      height    = 3,
      units     = "in",
      res       = 300)
  p=ggplot()+ geom_boxplot(data=data_use_merged_pick_plot,aes(x=days, y=.data[[family_pick]]),  fill="lightblue") + theme_minimal()+
    geom_rect(data=data_yagahi_supp_milkandsolid_day,aes(xmin=min(month),xmax=max(month),ymin=0,ymax=Inf),color="gray50",fill='pink',alpha=0.05)+
    geom_rect(data=data_yagahi_supp_solid_day,aes(xmin=min(month),xmax=max(month),ymin=0,ymax=Inf),color="gray50",fill='purple',alpha=0.01)+
    geom_smooth(data=data_use_merged_pick_plot,aes(x=day/30, y=.data[[family_pick]]),color = "black",se=TRUE, size = 0.5,  level=0.95)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  print(p)
  dev.off()
  
  
}