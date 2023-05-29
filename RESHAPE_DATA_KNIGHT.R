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
library(data.table) 

# sampled daily 

data_knight_M5 = read_excel('KNIGHT_DATA/M3_feces_L5.xlsx',col_names = TRUE, skip=1)
data_knight_F5 = read_excel('KNIGHT_DATA/F4_feces_L5.xlsx',col_names = TRUE, skip=1)

######## FOR MALE SUBJECT
# TAKE THE STRING AFTER f_
data_knight_M5 = data_knight_M5 %>% dplyr::rowwise() %>% dplyr::mutate(family=sub(".*f\\__", "", `#OTU ID`))
data_knight_M5 = data_knight_M5 %>% filter(family!="")
# there are two rows for Erysipelotrichaceae, but that's not a concern for me (not in my focal families)
data_knight_M5 = data_knight_M5[,2:dim(data_knight_M5)[2]]
data_knight_M5 = data_knight_M5 %>% relocate(family)
# relative abundance?
numeric_cols  = names(data_knight_M5)[sapply(data_knight_M5, is.numeric)]
day_grid_M5   = as.numeric(numeric_cols) # sampling days
data_knight_M5[, numeric_cols] <- lapply(data_knight_M5[, numeric_cols], function(x) x/sum(x)) #normalize and scale to absolute abundance
# #### sanity check -> all 1
# mean(colSums(data_knight_M5[,numeric_cols]))

######## FOR FEMALE SUBJECT
# TAKE THE STRING AFTER f_
data_knight_F5 = data_knight_F5 %>% dplyr::rowwise() %>% dplyr::mutate(family=sub(".*f\\__", "", `#OTU ID`))
data_knight_F5 = data_knight_F5 %>% filter(family!="")
# there are two rows for Erysipelotrichaceae, but that's not a concern for me (not in my focal families)
data_knight_F5 = data_knight_F5[,2:dim(data_knight_F5)[2]]
data_knight_F5 = data_knight_F5 %>% relocate(family)
# relative abundance?
numeric_cols  = names(data_knight_F5)[sapply(data_knight_F5, is.numeric)]
day_grid_F5   = as.numeric(numeric_cols) # sampling days
data_knight_F5[, numeric_cols] <- lapply(data_knight_F5[, numeric_cols], function(x) x/sum(x)) #normalize and scale to absolute abundance

# #### sanity check -> all 1
# mean(colSums(data_knight_F5[,numeric_cols]))

# FILTER FOR FAMILIES

families = c('Bifidobacteriaceae','Enterobacteriaceae','Lachnospiraceae','Streptococcaceae',
             'Bacteroidaceae','Enterococcaceae','Staphylococcaceae','Peptostreptococcaceae',
             'Clostridiaceae','Ruminococcaceae')

data_knight_M5 = data_knight_M5 %>% filter(family %in% families)
data_knight_M5_keep = data_knight_M5
data_knight_M5_t    = as.data.frame(t(data_knight_M5))
colnames(data_knight_M5_t) = data_knight_M5_t[1,]
data_knight_M5_t = data_knight_M5_t[2:dim(data_knight_M5_t)[1],]

data_knight_F5 = data_knight_F5 %>% filter(family %in% families)
data_knight_F5_keep = data_knight_F5
data_knight_F5_t    = as.data.frame(t(data_knight_F5))
colnames(data_knight_F5_t) = data_knight_F5_t[1,]
data_knight_F5_t = data_knight_F5_t[2:dim(data_knight_F5_t)[1],]

# data_knight_M5_t, day_grid_M5
# data_knight_F5_t, day_grid_F5 






