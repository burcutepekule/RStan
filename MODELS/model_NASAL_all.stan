functions {
  real[] SEIR(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[]  x_i 
  ) {
    real IC_A_0     = x_r[1]; 
    real IC_A_89    = x_r[2];
    real IC_C_89    = x_r[3];
    real IC_A_34    = x_r[4];
    real IC_D_34    = x_r[5];
    real IC_A_M     = x_r[6];
    real IC_M_M     = x_r[7];
    int  dum        = x_i[1]; //the only fixed integer is the reduction in inf due to being chronic
    real dydt[8];
    real init[8]; 
    
    real mu_a;
    real mu_c;
    real mu_d;
    real alpha_aa;
    real alpha_cc;
    real alpha_dd;
    real alpha_ac;
    real alpha_ca;
    real alpha_ad;
    real alpha_da;
    real alpha_dc;
    real alpha_cd;

    // Free parameters
    mu_a       = theta[1];
    mu_c       = theta[2];
    mu_d       = theta[3];
    alpha_aa   = theta[4];
    alpha_cc   = theta[5];
    alpha_dd   = theta[6];
    alpha_ac   = theta[7];
    alpha_ca   = theta[8];
    alpha_ad   = theta[9];
    alpha_da   = theta[10];
    alpha_dc   = theta[11];
    alpha_cd   = theta[12];
    
    init[1] = IC_A_0;
    init[2] = IC_A_89;
    init[3] = IC_C_89;
    init[4] = IC_A_34;
    init[5] = IC_D_34;
    init[6] = IC_A_M;
    init[7] = 0.0476*IC_M_M; // (1/21)
    init[8] = 0.9520*IC_M_M; //(20/21)

    // EXPERIMENT 0 : S.aureus alone.
    // A (S.aureus)
    dydt[1] = +mu_a*(y[1]+init[1])-alpha_aa*(y[1]+init[1])*(y[1]+init[1]);

    // EXPERIMENT 89 : S.aureus side by side with C.pseudodiphtheriticum.
    // A (S.aureus) 
    dydt[2] = +mu_a*(y[2]+init[2])-alpha_aa*(y[2]+init[2])*(y[2]+init[2])-(alpha_ca)*(y[2]+init[2])*(y[3]+init[3]);
    // C (C.pseudodiphtheriticum)
    dydt[3] = +mu_c*(y[3]+init[3])-alpha_cc*(y[3]+init[3])*(y[3]+init[3])-(alpha_ac)*(y[2]+init[2])*(y[3]+init[3]);

    // EXPERIMENT 34 : S.aureus side by side with D.pigrum.
    // A (S.aureus)  
    dydt[4] = +mu_a*(y[4]+init[4])-alpha_aa*(y[4]+init[4])*(y[4]+init[4])-(alpha_da)*(y[4]+init[4])*(y[5]+init[5]);
    // D (D.pigrum)
    dydt[5] = +mu_d*(y[5]+init[5])-alpha_dd*(y[5]+init[5])*(y[5]+init[5])-(alpha_ad)*(y[4]+init[4])*(y[5]+init[5]);
    
    // EXPERIMENT MIX : S.aureus side by side with a mix of D.pigrum & C.pseudodiphtheriticum.
    // The ratio of commensals in the mix is 20:1, respectively.
    // A (S.aureus) 
    dydt[6] = +mu_a*(y[6]+init[6])-alpha_aa*(y[6]+init[6])*(y[6]+init[6])
    -(y[6]+init[6])*((alpha_ca)*(y[7]+init[7])+(alpha_da)*(y[8]+init[8]));
    // C (C.pseudodiphtheriticum)
    dydt[7] = +mu_c*(y[7]+init[7])-alpha_cc*(y[7]+init[7])*(y[7]+init[7])
    -(y[7]+init[7])*((alpha_ac)*(y[6]+init[6])+(alpha_dc)*(y[8]+init[8]));
    // D (D.pigrum)
    dydt[8] = +mu_d*(y[8]+init[8])-alpha_dd*(y[8]+init[8])*(y[8]+init[8])
    -(y[8]+init[8])*((alpha_ad)*(y[6]+init[6])+(alpha_cd)*(y[7]+init[7]));

    return(dydt);
  }
}

data {
  // data
  real IC_A_0;
  real IC_A_89;
  real IC_C_89;
  real IC_A_34;
  real IC_D_34;
  real IC_A_M;
  real IC_M_M;
  int dum;
  int D_0; // number of samples
  int D_89; // number of samples
  int D_34;
  int D_M;
  real k_SA_0[D_0];
  real k_SA_89[D_89];
  real k_SA_34[D_34];
  real k_SA_M[D_M];

  // priors
  real p_mu[2];
  real p_alpha_inter[2];
  real p_alpha_intra[2];
  real p_phi;

  // Simulation
  real t0; //starting time
  int t_data; //index of first sample
  int S; // total simulation time
  real ts[D_89]; // time bins -> give the longest one
  real ts_pred[S];
}

transformed data {
  real x_r[7]; 
  int  x_i[1]; 
  real init[8] = rep_array(1e-9,8); 
  x_r[1]  = IC_A_0;
  x_r[2]  = IC_A_89;
  x_r[3]  = IC_C_89;
  x_r[4]  = IC_A_34;
  x_r[5]  = IC_D_34;
  x_r[6]  = IC_A_M;
  x_r[7]  = IC_M_M;

  init[1] = x_r[1];
  init[2] = x_r[2];
  init[3] = x_r[3];
  init[4] = x_r[4];
  init[5] = x_r[5];
  init[6] = x_r[6];
  init[7] = 0.0476*x_r[7]; // (1/21)
  init[8] = 0.9520*x_r[7]; //(20/21)
  
  x_i[1]  = dum;
}

parameters{
  real<lower=0> mu_a; 
  real<lower=0> mu_c; 
  real<lower=0> mu_d; 
  real<lower=0> alpha_aa; 
  real<lower=0> alpha_cc; 
  real<lower=0> alpha_dd; 
  // real<lower=0> alpha_ac; 
  // real<lower=0> alpha_ca;
  real alpha_ac; 
  real alpha_ca;
  real alpha_ad; 
  real alpha_da;
  real alpha_dc; 
  real alpha_cd;
  real<lower=0,upper=0.01> phi_a[4]; // dispersion parameters
  // real<lower=0, upper=0.001> phi_a[4]; // dispersion parameters
  // real<lower=0, upper=0.05> phi_c; // dispersion parameters
}
transformed parameters {
  real theta[12]; // vector of parameterss
  real y[D_89,8]; // raw ODE output -> give the longest D
  vector[D_89] output_SA_0;
  vector[D_89] output_SA_89;
  vector[D_89] output_CP_89;
  vector[D_89] output_SA_34;
  vector[D_89] output_DP_34;
  vector[D_89] output_SA_M;
  vector[D_89] output_CP_M;
  vector[D_89] output_DP_M;

  theta = {mu_a,mu_c,mu_d,
  alpha_aa,alpha_cc,alpha_dd,
  alpha_ac,alpha_ca,alpha_ad,alpha_da,alpha_dc,alpha_cd};
  // run ODE solver
  y = integrate_ode_bdf(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances and mam_lockmum steps
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(i in 1:D_89) {
      output_SA_0[i]  = y[i,1];
      output_SA_89[i] = y[i,2];
      output_CP_89[i] = y[i,3];
      output_SA_34[i] = y[i,4];
      output_DP_34[i] = y[i,5];
      output_SA_M[i]  = y[i,6];
      output_CP_M[i]  = y[i,7];
      output_DP_M[i]  = y[i,8];
    }
}

model {
  // priors
  mu_a ~ normal(p_mu[1],p_mu[2]);
  mu_c ~ normal(p_mu[1],p_mu[2]);
  mu_d ~ normal(p_mu[1],p_mu[2]);
  alpha_aa ~ normal(p_alpha_inter[1],p_alpha_inter[2]);
  alpha_cc ~ normal(p_alpha_inter[1],p_alpha_inter[2]);
  alpha_dd ~ normal(p_alpha_inter[1],p_alpha_inter[2]);
  alpha_ac ~ normal(p_alpha_intra[1],p_alpha_intra[2]);
  alpha_ca ~ normal(p_alpha_intra[1],p_alpha_intra[2]);
  alpha_ad ~ normal(p_alpha_intra[1],p_alpha_intra[2]);
  alpha_da ~ normal(p_alpha_intra[1],p_alpha_intra[2]);
  alpha_dc ~ normal(p_alpha_intra[1],p_alpha_intra[2]);
  alpha_cd ~ normal(p_alpha_intra[1],p_alpha_intra[2]);
  phi_a ~ exponential(p_phi);
  // likelihood
  for(i in 1:D_0) {
    target += normal_lpdf( output_SA_0[i]  |k_SA_0[i],phi_a[1]);
    target += normal_lpdf( output_SA_89[i] |k_SA_89[i],phi_a[2]);
    target += normal_lpdf( output_SA_34[i] |k_SA_34[i],phi_a[3]);
    target += normal_lpdf( output_SA_M[i]  |k_SA_M[i],phi_a[4]);
  }
    for(i in (D_0+1):D_34) {
    target += normal_lpdf( output_SA_89[i] |k_SA_89[i],phi_a[2]);
    target += normal_lpdf( output_SA_34[i] |k_SA_34[i],phi_a[3]);
    target += normal_lpdf( output_SA_M[i]  |k_SA_M[i],phi_a[4]);
  }
  for(i in (D_34+1):D_M) {
    target += normal_lpdf( output_SA_89[i] |k_SA_89[i],phi_a[2]);
    target += normal_lpdf( output_SA_M[i]  |k_SA_M[i],phi_a[4]);
  }
  for(i in (D_M+1):D_89) {
    target += normal_lpdf( output_SA_89[i] |k_SA_89[i],phi_a[2]);
  }
}

generated quantities{
  
  real y_pred[S,8]; // raw ODE output

  vector[S] comp_SA_0;
  vector[S] comp_SA_89;
  vector[S] comp_CP_89;
  vector[S] comp_SA_34;
  vector[S] comp_DP_34;
  vector[S] comp_SA_M;
  vector[S] comp_CP_M;
  vector[S] comp_DP_M;

  y_pred = integrate_ode_bdf(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts_pred, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3);
    
    for(i in 1:S) {
      comp_SA_0[i]  = (y_pred[i,1]+x_r[1]);
      comp_SA_89[i] = (y_pred[i,2]+x_r[2]);
      comp_CP_89[i] = (y_pred[i,3]+x_r[3]);
      comp_SA_34[i] = (y_pred[i,4]+x_r[4]);
      comp_DP_34[i] = (y_pred[i,5]+x_r[5]);
      comp_SA_M[i]  = (y_pred[i,6]+x_r[6]);
      comp_CP_M[i]  = (y_pred[i,7]+0.0476*x_r[7]);
      comp_DP_M[i]  = (y_pred[i,8]+0.9520*x_r[7]);
    }
}


