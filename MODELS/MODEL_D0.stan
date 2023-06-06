functions {
  
  matrix to_matrix_rowwise(real[] v, int m, int n) {
    matrix[m, n] res;
    for (i in 1:n) {
      for (j in 1:m) {
        res[i, j] = v[(i - 1) * n + n];
      }
    }
    return res;
  }
  
  // Total abundance function
  // real function_total_abundance(real t) {
    //   
    //   // ////// I am hardcoding the output here for now, but this can an be input as well. 
    //   // data_time     = c(0,1,2,3,20,50,180) #days
    //   // data_abundance= c(10^3,10^4,10^6,10^8,10^9,10^10,10^10) #rDNA copies/g
    //   // df_abundance  = tibble(t = data_time, y = data_abundance)
    //   // fit_abundance = nls(y ~ SSlogis(t, Asym, xmid, scal), data = df_abundance)
    //   // coeffs_abundance = coef(fit_abundance)
    //   // Asym        = coeffs_abundance[[1]]
    //   // xmid        = coeffs_abundance[[2]]
    //   // scal        = coeffs_abundance[[3]]
    //   // t_abundance = seq(0,720)
    //   // totalAbundance    = Asym/(1+exp((xmid-t_abundance)/scal));
    //   // totalAbundance_df = as.data.frame(cbind(t_abundance,totalAbundance))
    //   // colnames(totalAbundance_df)=c('day','total_abundance')
    //   
    //   real Asym = 10015904335;
    //   real xmid = 28.29924;
    //   real scal = 3.77713;
    //   real totalAbundance;
    //   totalAbundance    = Asym/(1+exp((xmid-t)/scal));
    //   return(totalAbundance);
    // }
    
    // Switch to mixed food - introduction of certain bacterial taxa
    real switch_mixed(real t, int t_mixed) {
      return(0);
    }
    // Switch to solid food - introduction of certain bacterial taxa, end of milk
    real switch_solid(real t, int t_solid) {
      return(0);
    }
    
    // HMO concentration, dependent on breastfeeding
    // real function_HMO(real t, int binary_breastmilk, real yf_milk, real y0_milk, real logalpha_milk) 
    real function_HMO(real t, int binary_breastmilk) {
      // hardcode yf_milk y0_milk logalpha_milk (for now)
      
      real yf_milk = 4.775213;
      real y0_milk = 15.321442;
      real logalpha_milk = -3.791853;
      real HMO_level=0;
      
      // formula yf+(y0−yf)e−exp(logalpha)t
      if(binary_breastmilk==1){
        HMO_level = yf_milk+(y0_milk-yf_milk)*exp(-exp(logalpha_milk)*t);
      }else{
        // not zero of course, need to find an estimate for formula feeding
      }
      HMO_level = HMO_level/y0_milk; // max 1
      return(HMO_level);
    }
    
    // plant-O concentration, dependent on transition to mixed and solid food
    // real function_solid(real t, int t_mixed, int t_solid, real k_solid)
    real function_solid(real t, int t_mixed, int t_solid){
      // hardcode k_solid (for now)
      
      real k_solid = 0.015;
      real solid_level=0;
      if(t>=t_mixed){
        solid_level = 1/(1+exp(-k_solid*(t-t_solid)));
      }
      return(solid_level);
    }
    
    // Maternal IgA concentration, dependent on breastfeeding
    real function_mIgA(real t, int binary_breastmilk, real HMO_level) {
      real mIgA_mult = function_HMO(0,binary_breastmilk); // normalize to 1
      real mIgA_level=0;
      
      if(binary_breastmilk==1){
        mIgA_level = HMO_level/mIgA_mult;
      }else{
        // well this is zero (I think) - formula doesn't have IgA
      }
      return(mIgA_level);
    }
    // O2 concentration, dependent on bacterial consumption
    // real function_O2(real t, real O2_0, real k_O2) 
    real function_O2(real t) {
      // hardcode O2_0, k_O2 (for now)
      
      real O2_0 = 1;
      real k_O2 = 0.1;
      real O2_level=0;
      
      O2_level = O2_0*exp(-k_O2*t);
      return(O2_level);
    }
    
    // EGF
    // GAPs
    // M cells
    // TLR4 
    
    real[] ODE_MODEL(real t,
    real[] y,
    real[] theta,
    real[] x_r,
    int[]  x_i 
    ) {
      
      int numTaxa      = x_i[1]; // number of taxa
      int numTimeSteps = x_i[2]; // number of time steps with observed data
      int t_mixed      = x_i[3]; // time point of transition to mixed feeding
      int t_solid      = x_i[4]; // time point of transition to solid feeding
      int icheck_1     = x_i[5]; // coating ratio checkpoint 1 (10 days)
      int icheck_2     = x_i[6]; // coating ratio checkpoint 2 (30 days)
      int binary_breastmilk = x_i[7];
      
      real O2Dependency_vector[numTaxa]   = x_r[1:numTaxa];
      real HMODependency_vector[numTaxa]  = x_r[(1*numTaxa+1):(2*numTaxa)];
      real solidDependency_vector[numTaxa]= x_r[(2*numTaxa+1):(3*numTaxa)];
      real growthRate_vector[numTaxa]     = x_r[(3*numTaxa+1):(4*numTaxa)];
      real interactionMat_vector[numTaxa*numTaxa] = x_r[(4*numTaxa+1):(4*numTaxa+numTaxa*numTaxa)];
      
      real dydt[2*numTaxa]; // 2*numTaxa because one compartment is for coated, the other is uncoated
      
      
      // taxa_array
      // [1] "Bifidobacteriaceae"    "Enterobacteriaceae"    "Lachnospiraceae"       "Streptococcaceae"      "Bacteroidaceae"       
      // [6] "Enterococcaceae"       "Staphylococcaceae"     "Peptostreptococcaceae" "Clostridiaceae"        "Ruminococcaceae" 
      
      // The following part has to have explicit equations for each taxa since the growth rate of some depend on HMO 
      // and some depend on O2, etc. 
      
      // Growth and interaction parameters will be estimated prior to the infant model, by using the
      // adult microbiome data - here the argument is that during homeostasis things are at a steady state 
      // equilibrium  where one can integrate out the impact of immune reaction (since it will be constant over time)
      // This means that the growth parameters are estimated in an environment with low HMO and low O2 but high plant-O. 
      // So the presence of O2 and HMO will enhance the growth rate, whereas the lack of plant-O will reduce it.
      
      // Interaction matrix -> 10 x 10 (numTaxa x numTaxa)
      // One can define matricies in Stan, but given that I will have constraints on certain parameters (directionality)
      // I better define each parameter explicitly 
      
      
      // Differential Equation Indexing
      // 1:numTaxa -> bacterial taxa, uncoated
      // numTaxa+1:2*numTaxa -> bacterial taxa, coated
      
      real level_O2;
      real level_HMO;
      real level_solid;
      real level_mIgA;
      real total_abundance;
      
      real O2Dependency;
      real HMODependency;
      real solidDependency;
      
      real growthRate;
      real interactions;
      real growthRate_net;
      matrix[numTaxa,numTaxa] interactionMat;
      
      
      // as a toy model, estimate the coating vector?
      // Free parameters -> theta is the coating vector
      row_vector[numTaxa] interactionvectemp;
      row_vector[numTaxa] abundancevectemp;
      real y_uncoated[numTaxa];
      real y_coated[numTaxa];
      
      level_O2       = function_O2(t);
      interactionMat = to_matrix(interactionMat_vector, numTaxa, numTaxa); //convert array to matrix -> this works, checked.
      
      level_HMO       = function_HMO(t,binary_breastmilk);
      level_solid     = function_solid(t,t_mixed,t_solid);
      level_mIgA      = function_mIgA(t,binary_breastmilk,level_HMO); 
      // total_abundance = function_total_abundance(t);
      total_abundance = 1; // stan cannot work with 1E1, will try relative
      
      y_uncoated  = y[1:numTaxa];
      y_coated    = y[(numTaxa+1):(2*numTaxa)];
      
      for(k in 1:numTaxa){
        abundancevectemp[k] = y_uncoated[k]+y_coated[k];
      }
      
      for(tx in 1:numTaxa){
        
        O2Dependency   = O2Dependency_vector[tx]*level_O2;
        HMODependency  = HMODependency_vector[tx]*level_HMO;
        solidDependency= (1-solidDependency_vector[tx]*level_solid);
        growthRate     = growthRate_vector[tx];
        
        interactionvectemp = interactionMat[tx,1:numTaxa];
        
        interactions   = dot_product(interactionvectemp,abundancevectemp);
        // print(t," ",tx," ",O2Dependency," ",HMODependency," ",solidDependency);
        // print(t," ",tx," ",(1+O2Dependency+HMODependency-solidDependency));
        growthRate_net = (1+O2Dependency+HMODependency-solidDependency)*growthRate;
        
        if(round(t)==(10+t_solid)) {
          // print("t: ",t, "tx: ",tx," solidDependency: ", solidDependency," HMODependency: ", HMODependency, " O2Dependency: ", O2Dependency," growthRate: ", growthRate, " growthRate_net: ", growthRate_net, " theta: ", theta[tx], " level_mIgA: ",level_mIgA);
          // print("t: ",t, "tx: ",tx," level_solid: ", level_solid);
        }
        
        // UNCOATED
        dydt[tx]          = (1-theta[tx]*level_mIgA)*y[tx]*(growthRate_net + interactions/total_abundance);
        
        // COATED
        dydt[tx+numTaxa]  = (theta[tx]*level_mIgA)*y[tx+numTaxa]*(growthRate_net + interactions/total_abundance);
        
      }
      // print("dydt: ",dydt)
      // print("t: ",t)
      // print("y: ",y)
      return(dydt);
    }
}

data {
  
  // INPUTS
  int numTaxa; // number of taxa
  int numTimeSteps; // number of time steps with data
  int t_mixed; // time point of transition to mixed feeding
  int t_solid; // time point of transition to solid feeding
  int icheck_1; // coating ratio checkpoint 1 (10 days)
  int icheck_2; // coating ratio checkpoint 2 (30 days)
  int binary_breastmilk;
  
  real O2Dependency_vector[numTaxa];
  real HMODependency_vector[numTaxa];
  real solidDependency_vector[numTaxa];
  real growthRate_vector[numTaxa];
  real interactionMat_vector[numTaxa*numTaxa];
  real y0[2*numTaxa];
  real observations[numTimeSteps,numTaxa];
  
  // priors
  real p_coating_vector[2];
  real p_phi;
  
  // Simulation
  int  t0; //starting time
  int  t0_data; //index of first sample
  int  t_sim_end; //total simulation time
  real t_data[numTimeSteps]; // time bins of data
  real ts_pred[numTimeSteps]; // time bins of prediction (not doing prediction currently)
}

transformed data {
  
  real x_r[(4*numTaxa+numTaxa*numTaxa)]; 
  int  x_i[7]; 
  
  x_i[1]  = numTaxa;
  x_i[2]  = numTimeSteps;
  x_i[3]  = t_mixed;
  x_i[4]  = t_solid;
  x_i[5]  = icheck_1;
  x_i[6]  = icheck_2;
  x_i[7]  = binary_breastmilk;
  
  x_r[1:numTaxa]                 = O2Dependency_vector;
  x_r[(1*numTaxa+1):(2*numTaxa)] = HMODependency_vector;
  x_r[(2*numTaxa+1):(3*numTaxa)] = solidDependency_vector;
  x_r[(3*numTaxa+1):(4*numTaxa)] = growthRate_vector;
  x_r[(4*numTaxa+1):(4*numTaxa+numTaxa*numTaxa)] = interactionMat_vector;
  
}

parameters{
  real<lower=0,upper=1> coating_vector[numTaxa]; // coating parameters
  real<lower=0> phi[numTaxa]; // dispersion parameters
}

transformed parameters {
  real theta[numTaxa]; // vector of parameterss
  real y[numTimeSteps,2*numTaxa]; // raw ODE output 
  real output_mat[numTimeSteps,2*numTaxa];
  
  theta = coating_vector;
  // print("theta: ",theta)
  
  // run ODE solver
  y = integrate_ode_adams(
    ODE_MODEL, // ODE model
    y0, // initial states
    t0, // t0
    t_data, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances and mam_lockmum steps
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    
    for(ti in 1:numTimeSteps){
      for(tx in 1:(2*numTaxa)){
        output_mat[ti,tx]  = y[ti,tx];
      }
    }
    
    // print("y: ",y)
}

model {
  real uncoated;
  real coated;
  int  ti_1=icheck_1;
  int  ti_2=icheck_2;
  real coated_subject;
  real all_subject;
  real ratio_coated;
  
  // priors
  coating_vector ~ beta(p_coating_vector[1],p_coating_vector[2]);
  
  // likelihood - first for total abundances
  for(ti in 1:numTimeSteps){
    for(tx in 1:numTaxa){
      uncoated = output_mat[ti,tx];
      coated   = output_mat[ti,tx+numTaxa];
      // print(uncoated+coated)
      target += normal_lpdf( uncoated+coated | observations[ti,tx], phi[tx]);
    }
  }
  
  // likelihood - second for coating ratio
  // 10 days -> 80%
  // 30 days -> 50%
  // This is the total coated - we don't know how this is distributed among taxa, 
  // and finding this out is the point of this estimation.

  all_subject = 0;
  coated_subject = 0;
  for(tx in 1:numTaxa){
    coated_subject   = coated_subject + output_mat[ti_1,tx+numTaxa];
    all_subject      = all_subject +  output_mat[ti_1,tx] + output_mat[ti_1,tx+numTaxa];
  }
  ratio_coated = coated_subject/all_subject;
  target += normal_lpdf( ratio_coated | 0.8,0.01);
  
  
  all_subject = 0;
  coated_subject = 0;
  for(tx in 1:numTaxa){
    coated_subject   = coated_subject + output_mat[ti_2,tx+numTaxa];
    all_subject      = all_subject +  output_mat[ti_2,tx] + output_mat[ti_2,tx+numTaxa];
  }
  ratio_coated = coated_subject/all_subject;
  target += normal_lpdf( ratio_coated | 0.5,0.01);
  
}

generated quantities{
  
  real y_pred[numTimeSteps,2*numTaxa];// raw ODE output
  matrix[numTimeSteps,2*numTaxa] output_pred;
  
  y_pred = integrate_ode_adams(
    ODE_MODEL, // ODE model
    y0, // initial states
    t0, // t0
    ts_pred, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances and mam_lockmum steps
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(ti in 1:numTimeSteps){
      for(tx in 1:(2*numTaxa)){
        output_pred[ti,tx]  = y_pred[ti,tx];
      }
    }
}


