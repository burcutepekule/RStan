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


marino_data_growth_vec      = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/INTERACTION_DATA/sd01.xlsx', sheet='alphas')
marino_data_interaction_mat = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/INTERACTION_DATA/sd01.xlsx', sheet='beta_matrix')
marino_data_otu_mat         = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/INTERACTION_DATA/sd01.xlsx', sheet='OTU_matrix')


marino_data_interaction_mat_pick = marino_data_interaction_mat[c('OTUs','4','9','19','27')]

marino_data_interaction_mat_pick = marino_data_interaction_mat_pick %>% filter(OTUs %in% c(4,9,19,27))


###########
data_yiming   = read_excel('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/INTERACTION_DATA/genus.interaction.xlsx')


families = c('Enterobacteriaceae','Lactobacillaceae','Bifidobacteriaceae',
             'Bacteroidaceae','Lachnospiraceae','Clostridiaceae','Ruminococcaceae',
             'Streptococcaceae','Staphylococcaceae')

genus = c('Enterobacteriaceae','Lactobacillaceae','Bifidobacteriaceae',
             'Bacteroidaceae','Lachnospiraceae','Clostridiaceae','Ruminococcaceae',
             'Streptococcaceae','Staphylococcaceae')

all_genus = sort(unique(c(sort(unique(data_yiming$genus_main)),sort(unique(data_yiming$genus_adjacent)))))
data_genus_family = c()
data_genus_family$genus=all_genus
data_genus_family=as.data.frame(data_genus_family)
all_family = c('Lachnospiraceae','Lachnospiraceae','Bacteroidaceae','Bifidobacteriaceae','Lachnospiraceae','Lachnospiraceae','Lachnospiraceae','Coriobacteriaceae',
               'Lachnospiraceae','Lachnospiraceae','Coriobacteriaceae','Enterobacteriaceae','Ruminococcaceae','Lachnospiraceae','Eggerthellaceae','Lachnospiraceae',
               'Odoribacteraceae','Tannerellaceae','Comamonadaceae','Bacteroidaceae','Pseudomonadaceae','Lachnospiraceae','Ruminococcaceae','Coriobacteriaceae',
               'Coriobacteriaceae','Staphylococcaceae')
data_genus_family$family = all_family
data_genus_family_main = data_genus_family
data_genus_family_adj = data_genus_family
colnames(data_genus_family_main) = c('genus_main','family_main')
colnames(data_genus_family_adj) = c('genus_adjacent','family_adjacent')

data_yiming_family = merge(data_yiming,data_genus_family_main,by='genus_main')
data_yiming_family = merge(data_yiming_family,data_genus_family_adj,by='genus_adjacent')

data_yiming_family_keep = data_yiming_family
data_yiming_family = data_yiming_family[c('family_main','genus_main','family_adjacent','genus_adjacent',"Pvalue_postive_interaction","Pvalue_negative_interaction","log2FC_of_area_sqrt")]
data_yiming_family_pick = data_yiming_family %>% filter(genus_adjacent=='Anaerostipes' & genus_main=='Bifidobacterium')

####################

genus_pick = c('Escherichia/Shigella','Staphylococcus','Streptococcus','Lactobacillus','Bifidobacterium',
               'Bacteroides','Anaerostipes','Blautia','Clostridium IV','Clostridium XlVa','Roseburia',
               'Clostridium sensu stricto','Faecalibacterium','Ruminococcus')

data_yiming_family_pick = data_yiming_family %>% filter(genus_adjacent %in% genus_pick & genus_main %in% genus_pick)

data_yiming_family_pick = data_yiming_family_pick[c('genus_adjacent','genus_main','log2FC_of_area_sqrt',
                                                    'Pvalue_postive_interaction','Pvalue_negative_interaction',
                                                    'family_adjacent', 'family_main')]

data_yiming_family_pick2 = data_yiming_family %>% filter(family_adjacent=='Lachnospiraceae' & family_main=='Bifidobacteriaceae')
