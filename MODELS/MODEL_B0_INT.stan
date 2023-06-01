functions {
  real[] ODE_MODEL(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[]  x_i 
  ) {
    int numTaxa      = x_i[1]; // number of taxa
    int numTimeSteps = x_i[2]; // number of time steps with observed data
    real p_phi       = x_r[numTaxa+1];
    
    real dydt[numTaxa];
    // real scale = 1e-6; // NEED TO PUT THIS SCALING FACTOR HERE, DOESN'T WORK WITH DISTRIBUTIONS
    real scale; // LETS SEE IF IT MAKES SENSE TO CONNECT THIS TO PHI (DISPERSION PARAMETER)
    
    real y_use[numTaxa];
    real growthRate_vector[numTaxa];
    real interactionMat_vector_diag[numTaxa];
    real interactionMat_vector_nondiag[numTaxa*numTaxa-numTaxa];
    real interactionMat_vector[numTaxa*numTaxa];
    real phi_vector[numTaxa];
    matrix[numTaxa,numTaxa] interactionMat;
    real interactions_tx;
    row_vector[numTaxa] interactionvectemp;
    row_vector[numTaxa] abundancevectemp;
    real y_uncoated[numTaxa];
    real y_coated[numTaxa];
    int counter_diag;
    int counter_nondiag;

    // Free parameters : separate counter diag vs non-diag because diagonal terms will be negative for sure!
    growthRate_vector             = theta[1:numTaxa];
    interactionMat_vector_diag    = theta[(1+numTaxa):(numTaxa+numTaxa)];
    interactionMat_vector_nondiag = theta[(1+numTaxa+numTaxa):(numTaxa+numTaxa*numTaxa)];
    
    scale = p_phi;  // LETS SEE IF IT MAKES SENSE TO CONNECT THIS TO PHI (DISPERSION PARAMETER)

    // Build the matrix with non-diagonal and diagonal terms
    counter_diag    = 1;
    counter_nondiag = 1;
    for(r in 1:numTaxa){
      for(c in 1:numTaxa){
        if(r==c){
          interactionMat[r,c]=-1*scale*interactionMat_vector_diag[counter_diag];// will sample this from beta distribution so will be strictly positive
          counter_diag = counter_diag + 1;
        }else{
          interactionMat[r,c]=+1*scale*interactionMat_vector_nondiag[counter_nondiag];
          counter_nondiag = counter_nondiag + 1;
        }
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
  real y0[numTaxa];
  // real observations[numTimeSteps,numTaxa];
  int observations[numTimeSteps,numTaxa];
  
  // priors
  real p_mu[2];
  real p_a_intra[2];
  real p_a_inter[2];
  real p_phi;
  
  // Simulation
  int  t0; //starting time
  real t_data[numTimeSteps]; // time bins of data
  real ts_pred[numTimeSteps]; // time bins of prediction (not doing prediction currently)
}

transformed data {
  real x_r[numTaxa+1]; // y0 and the dispersion parameter 
  int  x_i[2]; 
  x_i[1]          = numTaxa;
  x_i[2]          = numTimeSteps;
  x_r[1:numTaxa]  = y0[1:numTaxa];
  x_r[numTaxa+1]  = p_phi;
}

parameters{
  real<lower=0> phi[numTaxa]; // dispersion parameters - BE CAREFUL WITH THIS WHEN YOU CHANGE THE SCALE OF DATA!
  real<lower=0> growthRate_vector[numTaxa]; // growth rates
  real<lower=0> interactionMat_vector_diag[numTaxa]; // -1* when building matrix
  real interactionMat_vector_nondiag[numTaxa*numTaxa-numTaxa]; // 
}

transformed parameters {
  real theta[numTaxa+numTaxa*numTaxa]; // vector of parameterss
  real y[numTimeSteps,numTaxa]; // raw ODE output 
  matrix[numTimeSteps,numTaxa] output_mat;

  theta[1:numTaxa] = growthRate_vector;
  theta[(1+numTaxa):(numTaxa+numTaxa)] = interactionMat_vector_diag;
  theta[(1+numTaxa+numTaxa):(numTaxa+numTaxa*numTaxa)] = interactionMat_vector_nondiag;
  
  // run ODE solver
  y = integrate_ode_bdf(
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
}

model {
  
  // priors
  growthRate_vector ~ normal(p_mu[1],p_mu[2]);
  interactionMat_vector_diag ~ normal(p_a_intra[1],p_a_intra[2]);
  interactionMat_vector_nondiag ~ normal(p_a_inter[1],p_a_inter[2]);
  phi ~ exponential(p_phi);
  
  // likelihood - first for total abundances
  for(ti in 1:numTimeSteps){
    for(tx in 1:numTaxa){
      // target += normal_lpdf(output_mat[ti,tx] | observations[ti,tx], phi[tx]);
      target += neg_binomial_2_lpmf( observations[ti,tx] | output_mat[ti,tx] , phi[tx]);
    }
  }
}

generated quantities{
  
  real y_pred[numTimeSteps,numTaxa];// raw ODE output
  matrix[numTimeSteps,numTaxa] output_pred;
  
  y_pred = integrate_ode_bdf(
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


