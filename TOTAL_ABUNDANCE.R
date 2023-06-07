
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
totalAbundance    = totalAbundance/Asym
totalAbundance_df = as.data.frame(cbind(t_abundance,totalAbundance))
colnames(totalAbundance_df)=c('day','total_abundance')
