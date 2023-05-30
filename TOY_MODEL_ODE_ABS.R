
library(deSolve)

closed.sir.model <- function (t, x, params) {
  ## first extract the state variables
  y_1 <- x[1]
  y_2 <- x[2]
  y_3 <- x[3]
  ## now extract the parameters
  mu_1 <- params["mu_1"]
  mu_2 <- params["mu_2"]
  mu_3 <- params["mu_3"]
  
  a_11 <- params["a_11"]
  a_12 <- params["a_12"]
  a_13 <- params["a_13"]
  a_21 <- params["a_21"]
  a_22 <- params["a_22"]
  a_23 <- params["a_23"]
  a_31 <- params["a_31"]
  a_32 <- params["a_32"]
  a_33 <- params["a_33"]

  dxdt_1 = y_1*(mu_1 - y_1*a_11 + y_2*a_21 + y_3*a_31);
  dxdt_2 = y_2*(mu_2 + y_1*a_12 - y_2*a_22 + y_3*a_32);
  dxdt_3 = y_3*(mu_3 + y_1*a_13 + y_2*a_23 - y_3*a_33);
  
  dxdt <- c(dxdt_1,dxdt_2,dxdt_3)
  
  ## return result as a list!
  list(dxdt)
}

parms = c(
  mu_1=0.2,
  mu_2=0.15,
  mu_3=0.25,
  a_11=1*1E-6,
  a_22=1*1E-6,
  a_33=1*1E-6,
  a_12=0.01*1E-6,
  a_13=-0.05*1E-6,
  a_21=-0.01*1E-6,
  a_23=0.01*1E-6,
  a_31=-0.02*1E-6,
  a_32=0.01*1E-6
)


# Chain 1: growthRate_vector: [0.214713,0.164672,4.9555]
# interactionMat_vector_diag: [1.68919,0.978767,4.78219]
# interactionMat_vector_nondiag: [-0.83365,-0.84416,1.51799,0.307746,1.19608,-0.785384]
# phi: [221457,872478,248303]

# parms = c(
#   mu_1=0.214713,
#   mu_2=0.164672,
#   mu_3=4.9555,
#   a_11=-1.68919,
#   a_22=-0.978767,
#   a_33=-4.78219,
#   a_12=-0.83365,
#   a_13=-0.84416,
#   a_21=1.51799,
#   a_23=0.307746,
#   a_31=1.19608,
#   a_32=-0.785384
# )


times <- seq(from=0,to=t_end,by=1)
xstart <- 10^6*c(y_1=0.05,y_2=0.5,y_3=0.2)

library(tidyverse)
ode(
  func=closed.sir.model,
  y=xstart,
  times=times,
  parms=parms
) %>% as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(linewidth=2)+
  theme_classic()+
  labs(x='time (yr)',y='number of individuals')

out_use        = out[,2:dim(out)[2]]
out_use_noisy  = as.data.frame(lapply(out_use, function(x) x+rnorm(1,0,0.001)))
out_init       = xstart

out_use_noisy_plot = out_use_noisy
out_use_noisy_plot$time = out$time

graphics.off()
out_use_noisy_plot %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(linewidth=2)+
  theme_classic()+
  labs(x='time (yr)',y='number of individuals')




