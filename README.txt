Progress 05/06/2023
---------------------------------------------------------------------------------------------------

Currently using B0_MODEL_RUN.R to estimate growth rates and interaction parameters.
C0_MODEL_RUN.R takes the output of B0_MODEL_RUN.R and estimates the coating ratios of maternal IgA.

Got the following error when I tried to run C0_MODEL_RUN.R using absolute abundances
Error in stanc(file = file, model_code = model_code, model_name = model_name,  :   0Semantic error in 'string', line 26, column 16 to column 27:Integer literal cannot be larger than 2_147_483_647.
  
Point is : need to consider another way to include absolute abundances, seems like stan cannot compute it like this - what is the point of absolute abundances anyway? Maybe I can try with relative?

Second problem :  I modify the growth rate with the following argument : (1+O2Dependency+HMODependency-solidDependency) - but seems like this can be zero time to time (and also 15 which is huge)

