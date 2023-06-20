functions {

  // Total abundance function
  real function_total_abundance(real t) {
    // ////// I am hardcoding the output here for now, but this can an be input as well.
    // data_time     = c(0,1,2,3,20,50,180) #days
    // data_abundance= c(10^3,10^4,10^6,10^8,10^9,10^10,10^10) #rDNA copies/g
    // df_abundance  = tibble(t = data_time, y = data_abundance)
    // fit_abundance = nls(y ~ SSlogis(t, Asym, xmid, scal), data = df_abundance)
    // coeffs_abundance = coef(fit_abundance)
    // Asym        = coeffs_abundance[[1]]
    // xmid        = coeffs_abundance[[2]]
    // scal        = coeffs_abundance[[3]]
    // t_abundance = seq(0,720)
    // totalAbundance    = Asym/(1+exp((xmid-t_abundance)/scal));
    // totalAbundance_df = as.data.frame(cbind(t_abundance,totalAbundance))
    // colnames(totalAbundance_df)=c('day','total_abundance')
    // 
    // real Asym = 10015904335; // CANNOT HANDLE VALUES LIKE Asym = 10015904335
    real Asym = 1;
    real xmid = 28.29924;
    real scal = 3.77713;
    real totalAbundance;
    totalAbundance    = Asym/(1+exp((xmid-t)/scal));
    return(totalAbundance);
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

  
  real function_TLR4(real t, int binary_premature, real level_mEGF) {
    // MODIFY THIS
    
    real TLR4_0 = 1;
    real k_TLR4 = 0.1;
    real TLR4_level=0;
    
    TLR4_level = TLR4_0*exp(-k_TLR4*t);
    return(TLR4_level);
  }
  
  real function_TLR9(real t, int binary_premature) {
    // MODIFY THIS
    
    real TLR4_0 = 1;
    real k_TLR4 = 0.1;
    real TLR4_level=0;
    
    TLR4_level = TLR4_0*exp(-k_TLR4*t);
    return(TLR4_level);
  }
  
  
  
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
    int binary_premature  = x_i[8];
    int binary_csection   = x_i[9];

    real O2Dependency_vector[numTaxa]   = x_r[1:numTaxa];
    real HMODependency_vector[numTaxa]  = x_r[(1*numTaxa+1):(2*numTaxa)];
    real solidDependency_vector[numTaxa]= x_r[(2*numTaxa+1):(3*numTaxa)];
    real TLR9_vector                    = x_r[(3*numTaxa+1):(4*numTaxa)];
    real TLR4_vector                    = x_r[(4*numTaxa+1):(5*numTaxa)];
    real growthRate_vector[numTaxa]     = x_r[(5*numTaxa+1):(6*numTaxa)];
    real interactionMat_vector[numTaxa*numTaxa] = x_r[(6*numTaxa+1):(6*numTaxa+numTaxa*numTaxa)];
    real IgA_halflife                           = x_r[(6*numTaxa+numTaxa*numTaxa)+1];

    real dydt[2*numTaxa+1]; // 2*numTaxa because one compartment is for coated, 
    // the other is uncoated, last component tracks the general inflammatory status
    
    // levels, IgA halflife
    real level_O2;
    real level_HMO;
    real level_solid;
    real level_mIgA;
    real level_TLR4;
    real level_TLR9;
    row_vector[numTaxa] TLR9_vectorized;
    row_vector[numTaxa] TLR4_vectorized;
    real TLR9_dot_product;
    real TLR4_dot_product;

    // dependencies
    real O2Dependency;
    real HMODependency;
    real solidDependency;
    
    // inflammation
    real mIgA_coating;
    real biasInflammation;
    real TLR9stimulation_apical;
    real TLR4stimulation_apical;
    real TLR4stimulation_basolateral;

    // growth and interactions
    real growthRate;
    real interactions;
    real growthRate_net;
    matrix[numTaxa,numTaxa] interactionMat;
    row_vector[numTaxa] interactionvectemp;
    row_vector[numTaxa] abundancevectemp_both;
    row_vector[numTaxa] abundancevectemp_uncoated;
    real total_abundance;
    real totalInflammation;
    real y_total;
    real dy_total;
    int index_uncoated;
    int index_coated;
    
    level_O2          = function_O2(t);
    level_HMO         = function_HMO(t,binary_breastmilk);
    level_solid       = function_solid(t,t_mixed,t_solid);
    level_mIgA        = function_mIgA(t,binary_breastmilk,level_HMO); 
    level_TLR4        = function_TLR4(t);
    level_TLR9        = function_TLR9(t);

    total_abundance   = function_total_abundance(t);
    interactionMat    = to_matrix(interactionMat_vector, numTaxa, numTaxa); //convert array to matrix -> this works, checked.
    // total_abundance   = 1; // stan cannot work with 1E1, will try relative
    
    totalInflammation = y[2*numTaxa+1];
    
    for(k in 1:numTaxa){
      abundancevectemp_both[k]     = y[k]+y[k+numTaxa];
      abundancevectemp_uncoated[k] = y[k];
      TLR9_vectorized              = TLR9_vector[k];
      TLR4_vectorized              = TLR4_vector[k];
    }
    
    mIgA_coating     = theta[1:numTaxa];
    // in between are O2, HMO, and solid dependencies.
    biasInflammation = theta[4*numTaxa+1];
    
    for(tx in 1:numTaxa){
      
      mIgA_coating   = theta[0*numTaxa+tx];
      O2Dependency   = O2Dependency_vector[tx]*level_O2*theta[1*numTaxa+tx];
      HMODependency  = HMODependency_vector[tx]*level_HMO*theta[2*numTaxa+tx];
      solidDependency= (1-solidDependency_vector[tx]*level_solid)*theta[3*numTaxa+tx];
      growthRate     = (1/(1+totalInflammation))*growthRate_vector[tx];
      growthRate_net = (1+O2Dependency+HMODependency-solidDependency)*growthRate;
      
      interactionvectemp = interactionMat[tx,1:numTaxa];
      interactions       = dot_product(interactionvectemp,abundancevectemp_both);
      y_total            = y[tx]+y[tx+numTaxa];
      dy_total           = y_total*(growthRate_net + interactions/total_abundance) ;
      
      index_uncoated = 0*numTaxa+tx;
      index_coated   = 1*numTaxa+tx;
      
      // UNCOATED, LUMEN
      dydt[index_uncoated] = + dy_total
                             + y[index_coated]*IgA_halflife // degredation of IgA
                             - y[index_uncoated]*level_mIgA*mIgA_coating[tx]; //mIgA coating
                             
       // UNCOATED, LUMEN
      dydt[index_uncoated] = + 0
                             - y[index_coated]*IgA_halflife // degredation of IgA
                             + y[index_uncoated]*level_mIgA*mIgA_coating[tx]; //mIgA coating                        

      // Assume that the sampled number of bacteria is negligible compared to the amount in the lumen so not subtracted from that.

    }
    
    inflammationCoeff= 1-1/(1+biasInflammation*totalInflammation);
    TLR9_dot_product = dot_product(TLR9_vectorized,abundancevectemp_both);
    TLR4_dot_product = dot_product(TLR4_vectorized,abundancevectemp_uncoated);

    TLR9stimulation_apical      = (1-inflammationCoeff)*level_TLR9*TLR9_dot_product; // uncoated+coated stimulating TLR9. Since this is not about surface but CpG motif load.
    TLR4stimulation_apical      = (1-inflammationCoeff)*level_TLR4*TLR4_dot_product*(1/TLR9stimulation_apical); // uncoated stimulating TLR4, supressed by TLR9
    TLR4stimulation_basolateral = inflammationCoeff*level_TLR4*TLR4_dot_product; // uncoated stimulating TLR4
    // total inflammation
    dydt[2*numTaxa+1]           = +TLR4stimulation_apical+TLR4stimulation_basolateral-TLR9stimulation_apical;
    
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
  int binary_premature;
  int binary_csection;

  real O2Dependency_vector[numTaxa];
  real HMODependency_vector[numTaxa];
  real solidDependency_vector[numTaxa];
  real TLR9_vector[numTaxa];
  real TLR4_vector[numTaxa];
  real growthRate_vector[numTaxa];
  real interactionMat_vector[numTaxa*numTaxa];
  real IgA_halflife;
  
  real y0[2*numTaxa+1];
  real observations[numTimeSteps,numTaxa];
  
  // priors
  real p_coating_vector[2];
  real p_O2_vector[2];
  real p_HMO_vector[2];
  real p_solid_vector[2];
  real p_bias_inflammation[2];
  real p_phi;
  
  // Simulation
  int  t0; //starting time
  int  t0_data; //index of first sample
  int  t_sim_end; //total simulation time
  real t_data[numTimeSteps]; // time bins of data
  real ts_pred[numTimeSteps]; // time bins of prediction (not doing prediction currently)
}

transformed data {
  
  real x_r[(6*numTaxa+numTaxa*numTaxa)+1]; 
  int  x_i[9]; 
  
  x_i[1]  = numTaxa;
  x_i[2]  = numTimeSteps;
  x_i[3]  = t_mixed;
  x_i[4]  = t_solid;
  x_i[5]  = icheck_1;
  x_i[6]  = icheck_2;
  x_i[7]  = binary_breastmilk;
  x_i[8]  = binary_premature;
  x_i[9]  = binary_csection;

  x_r[1:numTaxa]                 = O2Dependency_vector;
  x_r[(1*numTaxa+1):(2*numTaxa)] = HMODependency_vector;
  x_r[(2*numTaxa+1):(3*numTaxa)] = solidDependency_vector;
  x_r[(3*numTaxa+1):(4*numTaxa)] = TLR9_vector;
  x_r[(4*numTaxa+1):(5*numTaxa)] = TLR4_vector;
  x_r[(5*numTaxa+1):(6*numTaxa)] = growthRate_vector;
  x_r[(6*numTaxa+1):(6*numTaxa+numTaxa*numTaxa)] = interactionMat_vector;
  x_r[(6*numTaxa+numTaxa*numTaxa)+1] = IgA_halflife;
}

parameters{
  real<lower=0,upper=1> coating_vector[numTaxa]; // coating parameters
  real<lower=0,upper=1> O2_vector[numTaxa]; // dependencies
  real<lower=0,upper=1> HMO_vector[numTaxa]; // dependencies
  real<lower=0,upper=1> solid_vector[numTaxa]; // dependencies
  real<lower=0> bias_inflammation; // bias_inflammation
  real<lower=0> phi[numTaxa]; // dispersion parameters
}

transformed parameters {
  real theta[4*numTaxa+1]; // vector of parameterss
  real y[numTimeSteps,2*numTaxa]; // raw ODE output 
  real output_mat[numTimeSteps,2*numTaxa];
  
  theta[(0*numTaxa+1):(1*numTaxa)] = coating_vector;
  theta[(1*numTaxa+1):(2*numTaxa)] = O2_vector;
  theta[(2*numTaxa+1):(3*numTaxa)] = HMO_vector;
  theta[(3*numTaxa+1):(4*numTaxa)] = solid_vector;
  theta[1+(4*numTaxa)] = bias_inflammation;

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
  O2_vector ~ beta(p_O2_vector[1],p_O2_vector[2]);
  HMO_vector ~ beta(p_HMO_vector[1],p_HMO_vector[2]);
  solid_vector ~ beta(p_solid_vector[1],p_solid_vector[2]);
  bias_inflammation ~ beta(p_bias_inflammation[1],p_bias_inflammation[2]);
  phi ~ exponential(p_phi);

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


