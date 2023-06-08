Progress 08/06/2023
---------------------------------------------------------------------------------------------------

Check model E0, a lot of updating went in there. In a nutshell, I am conceptualizing the reactivity and the immune thresholds rather than dealing with DC and T cell compartments, since that would require way more parameters and I don't have data to fit them.


Progress 07/06/2023
---------------------------------------------------------------------------------------------------

Ok so I will NEED a metric of total abundance, since it related to the aggressiveness of the IgA response. Since the algorithm cannot work with numbers like 10ˆ11, why not normalize that to 1 instead - and assume 1 is the number of max bacteria in the lumen?

In terms of sampling -> this is numerically constrained by the number of DCs. Given the DC capacity, if the amount sampled is below, leave it as it is. If not, truncate. Do not give priority for now.
(Priority would be, say you have the capacity of 5, and you have 1 1 3 - sure, you sample them all. But if you have 3 3 9, do you sample 0 0 5, or do you sample 1 1 3?)


Progress 06/06/2023
---------------------------------------------------------------------------------------------------

Two things

1) Cannot estimate coating parameters for more than 30 days - so even if I do this partial fitting, data should end at day 30

2) 02 and HMO dependency will have to have a taxa dependent modifier, othwerwise the abundance becomes dependent only on the coating ratio, which does not make sense at all.

To implement both of these, I will move from model version C0 to D0.

Side note, maybe I should focus on writing down the whole model first to focus on what can be estimated where and have a better overview of the parameter space.

Progress 05/06/2023
---------------------------------------------------------------------------------------------------

Currently using B0_MODEL_RUN.R to estimate growth rates and interaction parameters.
C0_MODEL_RUN.R takes the output of B0_MODEL_RUN.R and estimates the coating ratios of maternal IgA.

Got the following error when I tried to run C0_MODEL_RUN.R using absolute abundances
Error in stanc(file = file, model_code = model_code, model_name = model_name,  :   0Semantic error in 'string', line 26, column 16 to column 27:Integer literal cannot be larger than 2_147_483_647.
  
Point is : need to consider another way to include absolute abundances, seems like stan cannot compute it like this - what is the point of absolute abundances anyway? Maybe I can try with relative?

Second problem :  I modify the growth rate with the following argument : (1+O2Dependency+HMODependency-solidDependency) - but seems like this can be zero time to time (and also 15 which is huge)

