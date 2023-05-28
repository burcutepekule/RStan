### Changes in HMO Concentrations throughout Lactation: Influencing Factors, Health Effects and Opportunities

### Studies show that the HMO concentrations are highest in colostrum (average 9–22 g/L), 
### followed by slightly lower concentrations in transitional milk (assume 1wk here) (average 8–19 g/L), 
### with a gradual decline in mature milk as lactation progresses, from 6–15 g/L in breast milk 
### collected within one month of birth, to 4–6 g/L after 6 months. 
library(tidyverse)
library(broom)

data_milk_time          = c(0,7,28,168) #days
data_milk_concentration = c((22+9)/2, (8+19)/2, (6+15)/2, (4+6)/2) #g/L
df_milk       = tibble(t = data_milk_time, y = data_milk_concentration)
fit_milk      = nls(y ~ SSasymp(t, yf, y0, log_alpha), data = df_milk)
coeffs_milk   = coef(fit_milk)
yf_milk       = coeffs_milk[[1]]
y0_milk       = coeffs_milk[[2]]
logalpha_milk = coeffs_milk[[3]]


