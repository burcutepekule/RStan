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
library(pracma)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


sample

theme_set(theme_bw())
mc.cores = parallel::detectCores()
print(c('cores: ',mc.cores))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

qsum = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))

logit = function(x) log(x/(1-x))

inv.logit = function(x) exp(x)/(1+exp(x))

todaystr            = format(Sys.Date(), "%d%m%Y");
direc2save          = paste0("OUT/",todaystr,"/")
mkdir(direc2save)

