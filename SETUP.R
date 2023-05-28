# Set-up
library(tidyverse)
library(lubridate)
library(rstan)
library(cowplot)
library(readxl)
library(foreign)
library(xtable)
library(rstan)
library(zoo)
library(Rcpp)
library(reticulate)
library(icesTAF)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

theme_set(theme_bw())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

qsum = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))

logit = function(x) log(x/(1-x))

inv.logit = function(x) exp(x)/(1+exp(x))

prename_projections = "model_predictions_";
prename_figures     = "model_figures_";
prename_parameters  = "model_parameters_";
todaystr            = format(Sys.Date(), "%d%m%Y");
direc2save          = paste0("OUT/",todaystr,"/")
mkdir(direc2save)

