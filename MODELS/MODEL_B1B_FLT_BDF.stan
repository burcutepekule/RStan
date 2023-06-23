functions {
  real[] ODE_MODEL(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[]  x_i 
  ) {
    int numTaxa      = x_i[1]; // number of taxa
    int numTimeSteps = x_i[2]; // number of time steps with observed data
    int numAgnostic  = x_i[3];
    int numGnostic   = x_i[4];
    // real p_phi       = x_r[numTaxa+1];
    // real p_phi       = x_r[(numTaxa+1):2*numTaxa];

    real dydt[numTaxa];
    real growthRate_vector[numTaxa];
    real interactionMat_vector_Agnostic[numAgnostic];
    real interactionMat_vector_Gnostic[numGnostic];
    
    int interactionMask_vector[numTaxa*numTaxa];
    
    matrix[numTaxa,numTaxa] interactionMat;
    row_vector[numTaxa] interactionvectemp;
    row_vector[numTaxa] abundancevectemp;
    real interactions_tx;
    int counter_agnostic;
    int counter_gnostic;
    int counter_mask;
    int mask;
    real scale; 
    
    for(i in 1:numTaxa*numTaxa){
      interactionMask_vector[i] = x_i[4+i];
    }
    
    growthRate_vector               = theta[1:numTaxa];
    interactionMat_vector_Agnostic  = theta[(1+numTaxa):(numTaxa+numAgnostic)];
    interactionMat_vector_Gnostic   = theta[(numTaxa+numAgnostic+1):(numTaxa+numAgnostic+numGnostic)];

    counter_agnostic = 1;
    counter_gnostic  = 1;
    counter_mask     = 1;
    
    for(r in 1:numTaxa){
      for(c in 1:numTaxa){
        mask = interactionMask_vector[counter_mask];
        if(mask==0){
          interactionMat[r,c] = interactionMat_vector_Agnostic[counter_agnostic];
          counter_agnostic    = counter_agnostic + 1;
        }else{
          interactionMat[r,c] = mask*interactionMat_vector_Gnostic[counter_gnostic];
          counter_gnostic     = counter_gnostic + 1;
        }
        counter_mask = counter_mask + 1; 
      }
    }
    
    for(k in 1:numTaxa){
      abundancevectemp[k] = y[k];
    }
    
    for(tx in 1:numTaxa){
      interactionvectemp = interactionMat[tx,1:numTaxa];
      interactions_tx    = dot_product(interactionvectemp,abundancevectemp); // checked, this works 
      dydt[tx]           = y[tx]*(growthRate_vector[tx] + interactions_tx);
    }
    return(dydt);
  }
}

data {
  
  // INPUTS
  int numTaxa; // number of taxa
  int numTimeSteps; // number of time steps with data
  int numAgnostic;
  int numGnostic;
  int interactionMask_vector[numTaxa*numTaxa];
  real y0[numTaxa];
  real observations[numTimeSteps,numTaxa];
  // int observations[numTimeSteps,numTaxa];
  
  // priors
  real p_mu[2];
  real p_a[2];
  real p_phi[numTaxa];
  
  // Simulation
  int  t0; //starting time
  real t_data[numTimeSteps]; // time bins of data
  real ts_pred[numTimeSteps]; // time bins of prediction (not doing prediction currently)
}

transformed data {
  // real x_r[numTaxa+1]; // y0 and the dispersion parameter 
  real x_r[2*numTaxa]; // y0 and the dispersion parameter 
  int  x_i[(4+numTaxa*numTaxa)]; 
  
  x_i[1]          = numTaxa;
  x_i[2]          = numTimeSteps;
  x_i[3]          = numAgnostic;
  x_i[4]          = numGnostic;
  x_i[5:(4+numGnostic+numAgnostic)] = interactionMask_vector;

  x_r[1:numTaxa]  = y0[1:numTaxa];
  // x_r[numTaxa+1]  = p_phi;
  x_r[(numTaxa+1):2*numTaxa]  = p_phi;

}

parameters{
  real<lower=0> growthRate_vector[numTaxa]; // growth rates
  real interactionMat_vector_Agnostic[numAgnostic]; // interaction, directionality unknown
  real<lower=0> interactionMat_vector_Gnostic[numGnostic]; // interaction, directionality known
  real<lower=0> phi[numTaxa]; // dispersion parameters - BE CAREFUL WITH THIS WHEN YOU CHANGE THE SCALE OF DATA!
}

model {
  
  real theta[numTaxa+numTaxa*numTaxa]; // vector of parameterss
  real y[numTimeSteps,numTaxa]; // raw ODE output 
  matrix[numTimeSteps,numTaxa] output_mat;

  // priors
  growthRate_vector ~ normal(p_mu[1],p_mu[2]);
  interactionMat_vector_Agnostic ~ normal(p_a[1],p_a[2]);
  interactionMat_vector_Gnostic  ~ normal(p_a[1],p_a[2]);
  // phi ~ exponential(p_phi);
  for (i in 1:numTaxa){
    phi[i] ~ exponential(p_phi[i]);
  }

  theta[1:numTaxa] = growthRate_vector;
  theta[(1+numTaxa):(numTaxa+numAgnostic)]=interactionMat_vector_Agnostic;
  theta[(numTaxa+numAgnostic+1):(numTaxa+numAgnostic+numGnostic)]=interactionMat_vector_Gnostic;
    
  // run ODE solver
  y = integrate_ode_bdf(
    // y_pred = integrate_ode_rk45(
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
        for(tx in 1:numTaxa){
          output_mat[ti,tx]  = y[ti,tx];
        }
      }
      
      // likelihood - first for total abundances
      for(ti in 1:numTimeSteps){
        for(tx in 1:numTaxa){
          target += normal_lpdf(output_mat[ti,tx] | observations[ti,tx], phi[tx]);
          // target += neg_binomial_2_lpmf( observations[ti,tx] | output_mat[ti,tx] , phi[tx]);
        }
      }
}

generated quantities{
  
  real y_pred[numTimeSteps,numTaxa];// raw ODE output
  matrix[numTimeSteps,numTaxa] output_pred;
  real theta[numTaxa+numTaxa*numTaxa]; // vector of parameterss
  
  theta[1:numTaxa] = growthRate_vector;
  theta[(1+numTaxa):(numTaxa+numAgnostic)]=interactionMat_vector_Agnostic;
  theta[(numTaxa+numAgnostic+1):(numTaxa+numAgnostic+numGnostic)]=interactionMat_vector_Gnostic;
  
  y_pred = integrate_ode_bdf(
    // y_pred = integrate_ode_rk45(
      ODE_MODEL, // ODE model
      y0, // initial states
      t0, // t0
      ts_pred, // evaluation dates (ts)
      theta, // parameters
      x_r, // real data
      x_i, // integer data
      1.0E-10, 1.0E-10, 1.0E3); // tolerances and mam_lockmum steps
      // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
      for(tx in 1:numTaxa){
        for(ti in 1:numTimeSteps){
          output_pred[ti,tx]  = y_pred[ti,tx]; 
        }
      }
      
}


