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
    return(HMO_level);
  }
  
  // plant-O concentration, dependent on transition to mixed and solid food
  // real function_solid(real t, int t_mixed, int t_solid, real k_solid)
  real function_solid(real t, int t_mixed, int t_solid){
    // hardcode k_solid (for now)
    
    real k_solid = 0.5;
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
      mIgA_level = HMO_level*mIgA_mult;
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
  
  real ODE_MODEL(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[]  x_i 
  ) {
    
    int numSubjects  = x_i[1]; // number of subjects (assume same parameters apply for every subject)
    int numTaxa      = x_i[2]; // number of taxa
    int numTimeSteps = x_i[3]; // number of time steps
    int t_mixed[numSubjects]      = x_i[4:(4+numSubjects)]; // time point of transition to mixed feeding
    int t_solid[numSubjects]      = x_i[(4+numSubjects+1):(4+2*numSubjects+1)]; // time point of transition to solid feeding
    int binary_breastmilk[numSubjects] = x_i[(4+2*numSubjects+2):(4+3*numSubjects+2)];
    
    real O2Dependency_vector[numTaxa]   = x_r[1:numTaxa];
    real HMODependency_vector[numTaxa]  = x_r[(numTaxa+1):(2*numTaxa)];
    real growthRate_vector[numTaxa]     = x_r[(2*numTaxa+1):(3*numTaxa)];
    real interactionMat_vector[numTaxa*numTaxa] = x_r[(3*numTaxa+1):(3*numTaxa+numTaxa*numTaxa)];
    
    matrix[numSubjects,2*numTaxa] dydt; // 2*numTaxa because one compartment is for coated, the other is uncoated
    
    
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
    
    real O2Dependency;
    real HMODependency;
    real solidDependency;
    
    real growthRate;
    real interactions;
    real growthRate_net;
    matrix[numTaxa,numSubjects] interactionMat;
 

    // as a toy model, estimate the coating vector?
    // Free parameters -> theta is the coating vector
    row_vector[numTaxa] interactionvectemp;
    row_vector[numTaxa] abundancevectemp;
    real ytemp;


    
    level_O2       = function_O2(t);
    interactionMat = to_matrix_rowwise(interactionMat_vector, numTaxa, numSubjects); //convert array to matrix
    for (s in 1:numSubjects){
      
      level_HMO   = function_HMO(t,binary_breastmilk[s]);
      level_solid = function_solid(t,t_mixed[s],t_solid[s]);
      level_mIgA  = function_mIgA(t,binary_breastmilk[s],level_HMO); 
      
      for(tx in 1:numTaxa){
        
        O2Dependency   = O2Dependency_vector[tx]*level_O2;
        HMODependency  = HMODependency_vector[tx]*level_HMO;
        growthRate     = growthRate_vector[tx];
        ytemp          = y[s];
        
        interactionvectemp = interactionMat[tx,1:tx];
        abundancevectemp   = to_row_vector(ytemp);

        interactions   = dot_product(interactionvectemp,abundancevectemp);
        growthRate_net = (1+O2Dependency*level_O2+HMODependency*level_HMO-solidDependency*level_solid)*growthRate;
        
        // UNCOATED
        dydt[s,tx]          = theta[tx]*(1-level_mIgA)*(growthRate_net * y[s,tx] + interactions);
        
        // COATED
        dydt[s,tx+numTaxa]  = theta[tx]*(level_mIgA)*(growthRate_net * y[s,tx] + interactions);
      }
    }
    return(dydt);
  }
}

data {
  
  // INPUTS
  int numSubjects;// number of subjects (assume same parameters apply for every subject)
  int numTaxa; // number of taxa
  int numTimeSteps; // number of time steps
  int t_mixed; // time point of transition to mixed feeding
  int t_solid; // time point of transition to solid feeding
  int binary_breastmilk;
  
  real O2Dependency_vector[numTaxa];
  real HMODependency_vector[numTaxa];
  real growthRate_vector[numTaxa];
  real interactionMat_vector[numTaxa*numTaxa];
  real y0[numSubjects,2*numTaxa];
  real observations[numTimeSteps,numSubjects,numTaxa];
  
  // priors
  real p_coating_vector[2];
  
  // Simulation
  int  t0; //starting time
  int  t0_data; //index of first sample
  int  t_total; //total simulation time
  real t_data[numTimeSteps]; // time bins of data
  real ts_pred[numTimeSteps]; // time bins of prediction (not doing prediction currently)
}

transformed data {
  
  real x_r[(3*numTaxa+numTaxa*numTaxa)]; 
  int  x_i[6]; 
  
  x_i[1]  = numSubjects;
  x_i[2]  = numTaxa;
  x_i[3]  = numTimeSteps;
  x_i[4]  = t_mixed;
  x_i[5]  = t_solid;
  x_i[6]  = binary_breastmilk;
  
  x_r[1:numTaxa]                 = O2Dependency_vector;
  x_r[(numTaxa+1):(2*numTaxa)]   = HMODependency_vector;
  x_r[(2*numTaxa+1):(3*numTaxa)] = growthRate_vector;
  x_r[(3*numTaxa+1):(3*numTaxa+numTaxa*numTaxa)] = interactionMat_vector;
  
}

parameters{
  real<lower=0,upper=1> coating_vector[numTaxa]; // coating parameters
}

transformed parameters {
  real theta[(3*numTaxa+numTaxa*numTaxa)]; // vector of parameterss
  real y[numTimeSteps,numSubjects,2*numTaxa]; // raw ODE output 
  real output_mat[numTimeSteps,numSubjects,2*numTaxa];
  
  theta = coating_vector;
  
  // run ODE solver
  y = integrate_ode_bdf(
    ODE_MODEL, // ODE model
    y0, // initial states
    t0, // t0
    t_total, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances and mam_lockmum steps
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(ti in 1:t_total){
      for(s in 1:numSubjects){
        for(tx in 1:(2*numTaxa)){
          output_mat[ti,s,tx]  = y[ti,s,tx];
        }
      }
    }
}

model {
  // priors
  coating_vector ~ beta(p_coating_vector[1],p_coating_vector[2]);
  
  // likelihood - first for total abundances
  for(ti in 1:t_total){
    for(s in 1:numSubjects){
      for(tx in 1:numTaxa){
        uncoated = output_mat[ti,s,tx] 
        coated   = output_mat[ti,s,tx+numTaxa] 
        target += normal_lpdf( uncoated+coated | observations[ti,s,tx],0.01);
      }
    }
  }
  
  // likelihood - second for coating ratio
  // 10 days -> 80%
  // 30 days -> 50%
  // This is the total coated - we don't know how this is distributed among taxa, 
  // and finding this out is the point of this estimation.
  
  ti=10
  for(s in 1:numSubjects){
    coated_subject = 0;
    all_subject    = 0;
    for(tx in 1:numTaxa){
      coated_subject   = coated_subject + output_mat[ti,s,tx+numTaxa] 
      all_subject      = all_subject +  output_mat[ti,s,tx] + output_mat[ti,s,tx+numTaxa] 
    }
    ratio_coated = coated_subject/all_subject;
    target += normal_lpdf( ratio_coated | 0.8,0.01);
  }
  
  ti=30
  for(s in 1:numSubjects){
    coated_subject = 0;
    all_subject    = 0;
    for(tx in 1:numTaxa){
      coated_subject   = coated_subject + output_mat[ti,s,tx+numTaxa] 
      all_subject      = all_subject +  output_mat[ti,s,tx] + output_mat[ti,s,tx+numTaxa] 
    }
    ratio_coated = coated_subject/all_subject;
    target += normal_lpdf( ratio_coated | 0.5,0.01);
  }
  
}

generated quantities{
  
  real y_pred[numTimeSteps,numSubjects,2*numTaxa]// raw ODE output
  matrix[numTimeSteps,numSubjects,2*numTaxa] output_pred;
  
  y_pred = integrate_ode_bdf(
    ODE_MODEL, // ODE model
    y0, // initial states
    t0, // t0
    t_total, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances and mam_lockmum steps
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(ti in 1:t_total){
      for(s in 1:numSubjects){
        for(tx in 1:(2*numTaxa)){
          output_pred[ti,s,tx]  = y_pred[ti,s,tx];
        }
      }
    }
}


