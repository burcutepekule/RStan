functions {
  real[] ODE_MODEL(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[]  x_i 
  ) {
    real y_1_IC      = x_r[1];
    real y_2_IC      = x_r[2];
    real y_3_IC      = x_r[3];
    int numTaxa      = x_i[1]; // number of taxa
    int numTimeSteps = x_i[2]; // number of time steps with observed data
    real dydt[numTaxa]; //
    real IC[numTaxa]; 
    
    real mu_1;
    real mu_2;
    real mu_3;
    real a_11;
    real a_22;
    real a_33;
    real a_12;
    real a_13;
    real a_21;
    real a_23;
    real a_31;
    real a_32;
    
    real y_1;
    real y_2;
    real y_3;
    
    // Free parameters
    mu_1   = theta[1];
    mu_2   = theta[2];
    mu_3   = theta[3];
    a_11   = theta[4];
    a_22   = theta[5];
    a_33   = theta[6];
    a_12   = theta[7];
    a_13   = theta[8];
    a_21   = theta[9];
    a_23   = theta[10];
    a_31   = theta[11];
    a_32   = theta[12];
    
    IC[1]   = y_1_IC;
    IC[2]   = y_2_IC;
    IC[3]   = y_3_IC;
    
    y_1 = y[1]+IC[1];
    y_2 = y[2]+IC[2];
    y_3 = y[3]+IC[3];
    
    dydt[1]  = y[1]*(mu_1 - y[1]*a_11 + y[2]*a_21 + y[3]*a_31);
    dydt[2]  = y[2]*(mu_2 + y[1]*a_12 - y[2]*a_22 + y[3]*a_32);
    dydt[3]  = y[3]*(mu_3 + y[1]*a_13 + y[2]*a_23 - y[3]*a_33);
    
    return(dydt);
  }
}

data {
  
  // INPUTS
  int numTaxa; // number of taxa
  int numTimeSteps; // number of time steps with data
  real y0[numTaxa];
  real observations[numTimeSteps,numTaxa];
  
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
  real x_r[numTaxa]; 
  int  x_i[2]; 
  real init[numTaxa] = rep_array(1e-9,numTaxa); // initial values
  x_i[1]  = numTaxa;
  x_i[2]  = numTimeSteps;
  x_r[1]  = y0[1];
  x_r[2]  = y0[2];
  x_r[3]  = y0[3];
}

parameters{
  real<lower=0, upper=0.01> phi[numTaxa]; // dispersion parameters
  real<lower=0> mu_1; 
  real<lower=0> mu_2; 
  real<lower=0> mu_3; 
  real<lower=0> a_11; 
  real<lower=0> a_22; 
  real<lower=0> a_33; 
  real a_12; 
  real a_13;
  real a_21; 
  real a_23;
  real a_31; 
  real a_32;
}

transformed parameters {
  real theta[numTaxa+numTaxa*numTaxa]; // vector of parameterss
  real y[numTimeSteps,numTaxa]; // raw ODE output 
  matrix[numTimeSteps,numTaxa] output_mat;
  
  theta = {mu_1,mu_2,mu_3,a_11,a_22,a_33,a_12,a_13,a_21,a_23,a_31,a_32};
  
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
        // print(output_mat[ti,tx])
      }
    }
}

model {
  // priors
  mu_1 ~ normal(p_mu[1],p_mu[2]);
  mu_2 ~ normal(p_mu[1],p_mu[2]);
  mu_3 ~ normal(p_mu[1],p_mu[2]);
  a_11 ~ normal(p_a_intra[1],p_a_intra[2]);
  a_22 ~ normal(p_a_intra[1],p_a_intra[2]);
  a_33 ~ normal(p_a_intra[1],p_a_intra[2]);
  a_12 ~ normal(p_a_inter[1],p_a_inter[2]);
  a_13 ~ normal(p_a_inter[1],p_a_inter[2]);
  a_21 ~ normal(p_a_inter[1],p_a_inter[2]);
  a_23 ~ normal(p_a_inter[1],p_a_inter[2]);
  a_31 ~ normal(p_a_inter[1],p_a_inter[2]);
  a_32 ~ normal(p_a_inter[1],p_a_inter[2]);
  phi ~ exponential(p_phi);
  
  // likelihood - first for total abundances
  for(ti in 1:numTimeSteps){
    for(tx in 1:numTaxa){
      // print(output_mat[ti,tx])
      target += normal_lpdf(output_mat[ti,tx] | observations[ti,tx], phi[tx]);
      // target += neg_binomial_2_lpmf( observations[ti,tx] | output_mat[ti,tx] , phi[tx]);
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


