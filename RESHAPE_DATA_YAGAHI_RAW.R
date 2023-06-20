# run PREPARE_DATA_YAGAHI.R first if save_data_all.rds does not exist.
if(!file.exists("YAGAHI_DATA/save_data_all_RAW.rds")){
  source('PREPARE_DATA_YAGAHI_RAW.R')
}
saved_data_raw          = readRDS('YAGAHI_DATA/save_data_all_RAW.rds')
saved_data_smooth       = readRDS('YAGAHI_DATA/save_data_all_SMOOTH.rds')
saved_data_mean_smooth  = readRDS('YAGAHI_DATA/save_data_all_MEAN_SMOOTH.rds')
saved_data_milkandsolid = readRDS('YAGAHI_DATA/save_data_all_milkandsolid.rds')
saved_data_solid        = readRDS('YAGAHI_DATA/save_data_all_solid.rds')
time_grid_prediction    = readRDS('YAGAHI_DATA/day_grid_prediction.rds')

# REMOVE SUBJECT K -> DON'T KNOW WHEN TRANSITION TO SOLID FOOD HAPPENED
saved_data_raw          = saved_data_raw %>% filter(subject!='K')
saved_data_smooth       = saved_data_smooth %>% filter(subject!='K')
saved_data_milkandsolid = saved_data_milkandsolid %>% filter(Subject!='K')
saved_data_solid        = saved_data_solid %>% filter(Subject!='K')

# families = c('Bifidobacteriaceae','Enterobacteriaceae','Lachnospiraceae','Streptococcaceae',
#              'Bacteroidaceae','Enterococcaceae','Staphylococcaceae','Peptostreptococcaceae',
#              'Clostridiaceae','Ruminococcaceae')

families = c('Enterobacteriaceae','Lactobacillaceae','Bifidobacteriaceae',
             'Bacteroidaceae','Lachnospiraceae','Clostridiaceae','Ruminococcaceae',
             'Streptococcaceae','Staphylococcaceae')

saved_data_mean_smooth      = saved_data_mean_smooth %>% filter(taxa %in% families)
saved_data_mean_smooth      = saved_data_mean_smooth %>% rowwise() %>% mutate(abundance=ifelse(abundance<0,0,0.01*abundance))
saved_data_mean_smooth_wide = saved_data_mean_smooth %>% pivot_wider(names_from = taxa, values_from = abundance)
saved_data_mean_smooth_wide = na.omit(saved_data_mean_smooth_wide)

# max(saved_data_mean_smooth_wide$day)
# [1] 729
# update 

time_grid_prediction        = seq(0,719,1)
saved_data_mean_smooth_wide = saved_data_mean_smooth_wide %>% filter(day %in% time_grid_prediction)
abundanceArray_meanSubjects = saved_data_mean_smooth_wide
# abundanceArray_meanSubjects = saved_data_mean_smooth_wide[families]

