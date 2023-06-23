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
library(ggplot2)
library(knitr)
library(latex2exp)
library(png)
library(colorspace)
library(patchwork)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(data.table) 

# sampled daily 

data_knight_M6 = read_excel('KNIGHT_DATA/M3_feces_L6.xlsx',col_names = TRUE, skip=1)
data_knight_F6 = read_excel('KNIGHT_DATA/F4_feces_L6.xlsx',col_names = TRUE, skip=1)

######## FOR MALE SUBJECT
# TAKE THE STRING AFTER f_
data_knight_M6 = data_knight_M6 %>% dplyr::rowwise() %>% dplyr::mutate(genus=sub(".*g\\__", "", `#OTU ID`))
data_knight_M6 = data_knight_M6 %>% dplyr::rowwise() %>% dplyr::mutate(family=sub(";.*", "",sub(".*f\\__", "", `#OTU ID`)))
data_knight_M6 = data_knight_M6 %>% filter(family!="" & genus!="")
# there are two rows for Erysipelotrichaceae, but that's not a concern for me (not in my focal families)
data_knight_M6 = data_knight_M6[,2:dim(data_knight_M6)[2]]
data_knight_M6 = data_knight_M6 %>% relocate(c(family, genus))
# relative abundance?
numeric_cols  = names(data_knight_M6)[sapply(data_knight_M6, is.numeric)]
day_grid_M6   = as.numeric(numeric_cols) # sampling days
data_knight_M6[, numeric_cols] <- lapply(data_knight_M6[, numeric_cols], function(x) x/sum(x)) #normalize and scale to absolute abundance
# #### sanity check -> all 1
# mean(colSums(data_knight_M6[,numeric_cols]))

# multiply by 100 to meet yagahi 
# data_knight_M6[, numeric_cols] <- lapply(data_knight_M6[, numeric_cols], function(x) 100*x) #normalize and scale to absolute abundance

# FILTER FOR GENUS 

taxanomy_table = c()
taxanomy_table$genus = c('Escherichia','Staphylococcus','Streptococcus','Lactobacillus','Bifidobacterium',
                         'Bacteroides','Anaerostipes','Blautia','Eubacterium','Clostridium',
                         'Faecalibacterium','Ruminococcus')

taxanomy_table$family = c('Enterobacteriaceae','Staphylococcaceae','Streptococcaceae','Lactobacillaceae','Bifidobacteriaceae',
                          'Bacteroidaceae','Lachnospiraceae','Lachnospiraceae',"Lachnospiraceae",'Clostridiaceae',
                          'Ruminococcaceae','Ruminococcaceae') 
taxanomy_table = as.data.frame(taxanomy_table)
data_knight_M6_taxonomy = data_knight_M6[c('genus','family')]

# match_df(taxanomy_table, data_knight_M6_taxonomy, on = NULL)

data_knight_M6_merged = merge(taxanomy_table, data_knight_M6, by=c('genus','family'))

taxanomy_table = as.data.frame(taxanomy_table)
all_families = unique(taxanomy_table$family)

data_use_merged_long  = data_knight_M6_merged %>% pivot_longer(!c('genus','family'), names_to = "days", values_to = "abundance")

# DAILY! - > In this study, two healthy subjects, one male (M3) and one female (F4), 
# one of whom (M3) participated in an earlier survey [6], were sampled daily at three
# body sites (gut (feces), mouth, and skin (left and right palms)), for 15 months (M3) 
# and for 6 months (F4) using an institutional review board-approved protocol.

for(family_pick in all_families){
  print(family_pick)
  taxonomy_use = taxanomy_table %>% filter(family==family_pick)
  genus_use = taxonomy_use$genus
  data_use_merged_pick = data_knight_M6_merged %>% filter(family==family_pick)
  data_use_merged_pick_numeric = data_use_merged_pick[,3:dim(data_use_merged_pick)[2]]
  data_use_merged_pick_numeric = colSums(data_use_merged_pick_numeric)
  data_plot = c()
  data_plot$days = names(data_use_merged_pick_numeric)
  data_plot$fam = data_use_merged_pick_numeric
  data_plot     = as.data.frame(data_plot)
  colnames(data_plot)[2] = family_pick

  data_use_merged_pick_plot = data_plot[c('days',family_pick)]
  data_use_merged_pick_plot$days = as.numeric(data_use_merged_pick_plot$days)
  # 
  # graphics.off()
  # png(file =paste0("/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/KNIGHT_DATA/FIGS/FAM_AGG_2USE_",family_pick,".png"),   # The directory you want to save the file in
  #     width     = 8,
  #     height    = 3,
  #     units     = "in",
  #     res       = 300)
  # p=ggplot()+ geom_point(data=data_use_merged_pick_plot,aes(x=days, y=.data[[family_pick]]),  fill="lightblue") + theme_minimal()+
  #   geom_smooth(data=data_use_merged_pick_plot,aes(x=days, y=.data[[family_pick]]),color = "black",se=TRUE, size = 0.5,  level=0.95)+
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # 
  # print(p)
  # dev.off()

}

data_use_merged_long$days  = as.numeric(data_use_merged_long$days)
data_knight_M6_agg         = aggregate(abundance ~family + days,data_use_merged_long,FUN=sum)
data_knight_M6_wide        = data_knight_M6_agg%>%pivot_wider(names_from = family, values_from = abundance)
days_array = data_knight_M6_wide$days
abundanceArray_meanSubjects= data_knight_M6_wide[,2:dim(data_knight_M6_wide)[2]]
rownames(abundanceArray_meanSubjects) = days_array


y0_meanSubjects             = as.numeric(unlist(abundanceArray_meanSubjects[1,]))
y0_meanSubjects[which(y0_meanSubjects==0)] = abs(0.001*rnorm(length(which(y0_meanSubjects==0))))


abundanceArray_meanSubjects = abundanceArray_meanSubjects[2:dim(abundanceArray_meanSubjects)[1],] #start from day 2, day 0 is the initial condition
rownames(abundanceArray_meanSubjects) = days_array[2:length(days_array)]
taxa_array = colnames(abundanceArray_meanSubjects)


# ## also read the interaction masking matrix
interactionMask = read_excel('KNIGHT_DATA/INTERACTION_MASK.xlsx',col_names = TRUE, skip=0)

interactionMask_df              = as.data.frame(interactionMask)
interactionMask_df[is.na(interactionMask_df)] <- 0

colnames(interactionMask_df)[1] = 'effected'
interactionMask_long     = interactionMask_df %>% pivot_longer(!effected, names_to = "effector", values_to = "direction") 
interactionMask_long     = interactionMask_long %>% rowwise() %>% mutate(sd=ifelse(is.na(direction),4,2))
interactionMask_long     = interactionMask_long %>% rowwise() %>% mutate(taxa_index=10*which(effected==taxa_array)+which(effector==taxa_array))
interactionMask_long_sd  = unlist(interactionMask_long$sd)
interactionMask_long_idx = unlist(interactionMask_long$taxa_index)


interactionMask_df_reorder = interactionMask_df[,c('effected',taxa_array)]

rownames(interactionMask_df_reorder) = interactionMask_df_reorder$effected
interactionMask_df_reorder = interactionMask_df_reorder[taxa_array,]
interactionMask_df_reorder = interactionMask_df_reorder[,2:dim(interactionMask_df_reorder)[2]]
interactionMask_vector = as.vector(t(interactionMask_df_reorder))

numNegative = length(which(interactionMask_vector==-1))
numPositive = length(which(interactionMask_vector==+1))
numAgnostic = length(which(interactionMask_vector==0))
numGnostic  = numNegative+numPositive

