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
  
  real function_GAPs(real t) {
    // MODIFY THIS
    
    real GAPs_0 = 1;
    real k_GAPs = 0.1;
    real GAPs_level=0;
    
    GAPs_level = GAPs_0*exp(-k_GAPs*t);
    return(GAPs_level);
  }
  
  real function_MCells(real t) {
    // MODIFY THIS
    
    real MCells_0 = 1;
    real k_MCells = 0.1;
    real MCells_level=0;
    
    MCells_level = MCells_0*exp(-k_MCells*t);
    return(MCells_level);
  }
  
  real function_EGF(real t) {
    // MODIFY THIS
    
    real EGF_0 = 1;
    real k_EGF = 0.1;
    real EGF_level=0;
    
    EGF_level = EGF_0*exp(-k_EGF*t);
    return(EGF_level);
  }
  
  real function_TLR4(real t) {
    // MODIFY THIS
    
    real TLR4_0 = 1;
    real k_TLR4 = 0.1;
    real TLR4_level=0;
    
    TLR4_level = TLR4_0*exp(-k_TLR4*t);
    return(TLR4_level);
  }
  
  real function_nBCells(real t) { // Naive BCell influx
  // MODIFY THIS
  
  real nBCells_0 = 1;
  real k_nBCells = 0.1;
  real nBCells_level=0;
  
  nBCells_level = nBCells_0*exp(-k_nBCells*t);
  return(nBCells_level);
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
    
    real O2Dependency_vector[numTaxa]   = x_r[1:numTaxa];
    real HMODependency_vector[numTaxa]  = x_r[(1*numTaxa+1):(2*numTaxa)];
    real solidDependency_vector[numTaxa]= x_r[(2*numTaxa+1):(3*numTaxa)];
    real growthRate_vector[numTaxa]     = x_r[(3*numTaxa+1):(4*numTaxa)];
    real interactionMat_vector[numTaxa*numTaxa] = x_r[(4*numTaxa+1):(4*numTaxa+numTaxa*numTaxa)];
    

    real dydt[2*numTaxa]; // 2*numTaxa because one compartment is for coated, the other is uncoated
    
    // check paper ApplyingMathematicalToolstoAccelerateVaccine Development/ModelingShigellaImmuneDynamics for B cell death rates etc.
    
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
    real level_MCells;
    real level_GAPs;
    real level_nBcells;
    real total_abundance;
    
    real O2Dependency;
    real HMODependency;
    real solidDependency;
    
    // Carrying capacities of plasma, memory, and circulating germinal center B cells.
    real carryingCapacity_PC;
    real carryingCapacity_MB;
    real carryingCapacity_GCB;
    
    real growthRate;
    real interactions;
    real growthRate_net;
    matrix[numTaxa,numTaxa] interactionMat;
    
    
    real total_abundance_sum;
    real total_sampling_threshold;
    real total_sampled;
    vector[numTaxa] y_uncoated_lumen;
    vector[numTaxa] y_coated_lumen;
    vector[numTaxa] y_coated_over_uncoated_lumen;
    vector[numTaxa] y_coated_over_uncoated_SED;
    vector[numTaxa] y_coated_over_uncoated_LP;
    vector[numTaxa] y_is_abundant_lumen;
    vector[2*numTaxa] y_sampled_SED_IECs;
    vector[2*numTaxa] y_sampled_SED_MCells;
    vector[2*numTaxa] y_sampled_LP_IECs;
    vector[2*numTaxa] y_sampled_LP_GAPs;
    
    vector[numTaxa] delta_reactivity_SED;
    vector[numTaxa] delta_reactivity_MLN;
    
    row_vector[numTaxa] interactionvectemp;
    row_vector[numTaxa] abundancevectemp;
    
    y_uncoated_lumen  = y[1:numTaxa];
    y_coated_lumen    = y[(numTaxa+1):(2*numTaxa)];
    
    level_O2          = function_O2(t);
    level_HMO         = function_HMO(t,binary_breastmilk);
    level_solid       = function_solid(t,t_mixed,t_solid);
    level_mIgA        = function_mIgA(t,binary_breastmilk,level_HMO); 
    level_MCells      = function_MCells(t);
    level_GAPs        = function_GAPs(t);
    level_nBcells     = function_nBCells(t); // influx of naive B cells
    
    total_abundance   = function_total_abundance(t);
    interactionMat    = to_matrix(interactionMat_vector, numTaxa, numTaxa); //convert array to matrix -> this works, checked.
    // total_abundance   = 1; // stan cannot work with 1E1, will try relative
    
    // DETERMINE THE TOTAL ABUNDANCE VECTOR & ABUNDANCE AS PRESENT / ABSENT
    total_abundance_sum = 0;
    for(k in 1:numTaxa){
      abundancevectemp[k] = y_uncoated_lumen[k]+y_coated_lumen[k];
      total_abundance_sum = total_abundance_sum + abundancevectemp[k];
      if(abundancevectemp[k]>0){
        y_is_abundant_lumen[k] = 1;
      }else{
        y_is_abundant_lumen[k] = 0;
      }
    }
    
    // DETERMINE THE ***LUMEN*** COATED / UNCOATED RATIOS
    for(k in 1:numTaxa){
      if(y_coated_lumen[k]>0 && y_uncoated_lumen[k]>0){
        y_coated_over_uncoated_lumen[k] = y_coated_lumen[k]/(y_uncoated_lumen[k]+y_coated_lumen[k]);
      }else{ // not == in case machine precision or some shit
      y_coated_over_uncoated_lumen[k] = 0;
      }
    }
    
    // SAMPLING RATES
    real tau_UC_GAP; //UNOCATED, VIA GAPS
    real tau_UC_IEC; //UNOCATED, VIA IEC
    real tau_UC_M; //UNOCATED, VIA M-CELLS
    real tau_C_M; //COATED, VIA M-CELLS
    
    // CONSTRAINTS ON SAMPLING RATES
    // tau_C_M      << tau_UC_M     -> M-cell bias of IgA-Ag uptake 
    // gamma_C_M    << gamma_UC_M   -> Tolerogenic bias of IgA-Ag complex (less reactivity increase)
    // gamma_UC_GAP << gamma_UC_IEC -> Tolerogenic bias of GAP proteins (less reactivity increase)
    
    // DETERMINE THE AMOUNT SAMPLED FROM DIFFERENT PATHWAYS (IECs, MCells, GAPs)
    // Assume that the coated IgA-Ag complexes are only sampled by M cells, since the main purpose of such a complex is to keep 
    // the bacteria further from the epithelium, but M cells have special affinity to IgA for that non-inf. feedback loop.
    
    // Sample proportional to the relative abundances - If sampling exceeds a certain threshold, trim (since this is a physical)
    // capacity for the number of DCs available. 
    total_sampled = 0;
    for(k in 1:numTaxa){
      //k : index for uncoated
      //k+numTaxa : index for coated
      
      // SED (Sub epithelial dome)
      y_sampled_SED_IECs[k]           = tau_UC_IEC*(y_uncoated_lumen[k]/total_abundance_sum);
      y_sampled_SED_IECs[k+numTaxa]   = 0;
      
      y_sampled_SED_MCells[k]         = level_MCells*tau_UC_M*(y_uncoated_lumen[k]/total_abundance_sum);
      y_sampled_SED_MCells[k+numTaxa] = level_MCells*tau_C_M*(y_coated_lumen[k]/total_abundance_sum);
      
      // LP (Lamina Propria)
      y_sampled_LP_IECs[k]            = tau_UC_IEC*(y_uncoated_lumen[k]/total_abundance_sum);
      y_sampled_LP_IECs[k+numTaxa]    = 0;
      
      y_sampled_LP_GAPs[k]            = level_GAPs*tau_UC_GAP*(y_uncoated_lumen[k]/total_abundance_sum);
      y_sampled_LP_GAPs[k+numTaxa]    = 0;
      
      total_sampled = total_sampled +
      y_sampled_SED_IECs[k] + y_sampled_SED_IECs[k+numTaxa] +
      y_sampled_SED_MCells[k] + y_sampled_SED_MCells[k+numTaxa] + 
      y_sampled_LP_IECs[k] + y_sampled_LP_IECs[k+numTaxa]+
      y_sampled_LP_GAPs[k] + y_sampled_LP_GAPs[k+numTaxa];
    }
    
    
    // Sample proportional to the relative abundances - If sampling exceeds a certain threshold, trim (since this is a physical)
    // capacity for the number of DCs available. 
    
    if(total_sampled > total_sampling_threshold){
      y_sampled_SED_IECs   = total_sampling_threshold.*y_sampled_SED_IECs./total_sampled;
      y_sampled_SED_MCells = total_sampling_threshold.*y_sampled_SED_MCells./total_sampled;
      y_sampled_LP_IECs    = total_sampling_threshold.*y_sampled_LP_IECs./total_sampled;
      y_sampled_LP_GAPs    = total_sampling_threshold.*y_sampled_LP_GAPs./total_sampled;
    }
    
    
    // More epitopes sampled at the same time, less chances for each of them to exceed the threshold of "irritation"
    // or "aggresiveness" on their own. IMMUNE ACTIVATION THRESHOLD! 
    
    // They don't exceed this during the critical window, because there is a WALL or a SHIELD due to the maternal IgA. 
    
    // sIgA-Ag complex induces tolerance response whereas Ag alone induces inflammatory cascades. 
    // DCs exposed to S. flexneri alone expressed high levels of the inflammatory cytokines interleukin (IL)-12, 
    // tumor necrosis factor (TNF)-a, and KC, DCs incubated with sIgA-S. flexneri complexes had reduced inflammatory
    // cytokine production and instead produced high levels of the anti-inflammatory mediator transforming 
    // growth factor (TGF)-b. Thus, SIgA binding to pathogens can dampen harmful inflammatory processes while still 
    // neutralizing the offending microbe.
    
    // So you need 3 coefficients, 
    
    // Incubating DCs with, 
    // 1) Ag alone : most reactivity increase
    // 2) IgA-Ag complex : reduced reactivity increase
    // 3) Ag + GAP proteins : least reactivity increase (given they almost always lead to tolerance)
    
    // REACTIVITY INCREASE RATES
    real gamma_Ag_0;//Ag alone
    real gamma_Ag_IgA;//IgA-Ag complex
    real gamma_Ag_GAP;//Ag + GAP proteins
    // gamma_Ag_0 > gamma_Ag_IgA > gamma_Ag_GAP
    
    
    // I don't know how much differentiation will go on in the PPs versus MLNs, but PPs should not react to GAP samples
    // because those samples directly drain to the MLN from LP. So maybe I should disect the reactivity in the SED vs MLN.
    
    for(k in 1:numTaxa){
      // Ag alone via IECs and MCells, IgA-Ag complex via MCells
      delta_reactivity_SED[k] = gamma_Ag_0*(y_sampled_SED_IECs[k] + y_sampled_SED_MCells[k]) + gamma_Ag_IgA*(y_sampled_SED_MCells[k+numTaxa]); 
      
      // Ag alone via IECs, g + GAP proteins via GAPs
      delta_reactivity_MLN[k] = gamma_Ag_0*y_sampled_LP_IECs[k] + gamma_Ag_GAP*y_sampled_LP_GAPs[k];
    }
    
    // HERE - FIGURE OUT HOW TO USE DELTA_REACTIVITY TO TRANSLATE IT INTO THE PARAMETRIZATION OF THE FUZZY IGA CURVE
    // THAT PARAMETRIZATION IS YOUR IMMUNE ACTIVATION THREHSOLD!
    
    // Also need to decide on the memory vs plasma cell diff! Think about age of the cells - each delta_t is a day in your case (since parameters are per day) - so you should 
    // think about how many days it takes to leave the GC -> once you leave, what is the threshold of turning memory vs plasma?
    
    // MIGHT NOT NEED TO EVEN DEFINE TOL. VS EFF. DCs BECAUSE THAT REQUIRES ANOTHER LEVEL OF THRESHOLDING
    // CHECK PLOT_SIGNALS.m LINE 90-96, WHERE THIS CURVE IS PARAMETERIZED - THERE IS ONE PARAMETER, WHICH IS biasCoat - that's what u are looking for!
    
    for(tx in 1:numTaxa){
      
      O2Dependency   = O2Dependency_vector[tx]*level_O2*theta[1*numTaxa+tx];
      HMODependency  = HMODependency_vector[tx]*level_HMO*theta[2*numTaxa+tx];
      solidDependency= (1-solidDependency_vector[tx]*level_solid)*theta[3*numTaxa+tx];
      growthRate     = growthRate_vector[tx];
      
      interactionvectemp = interactionMat[tx,1:numTaxa];
      interactions       = dot_product(interactionvectemp,abundancevectemp);
      growthRate_net     = (1+O2Dependency+HMODependency-solidDependency)*growthRate;
      
      
      // UNCOATED, LUMEN
      dydt[0*numTaxa+tx] = (1-theta[tx]*level_mIgA)*y_uncoated_lumen[tx]*(growthRate_net + interactions/total_abundance);
      
      // COATED, LUMEN
      dydt[1*numTaxa+tx] = (theta[tx]*level_mIgA)*y_coated_lumen[tx]*(growthRate_net + interactions/total_abundance);
      
      // Assume that the sampled number of bacteria is negligible compared to the amount in the lumen so not subtracted from that.
      
      // SAMPLED, 
      
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
  real p_O2_vector[2];
  real p_HMO_vector[2];
  real p_solid_vector[2];
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
  real<lower=0,upper=1> O2_vector[numTaxa]; // coating parameters
  real<lower=0,upper=1> HMO_vector[numTaxa]; // coating parameters
  real<lower=0,upper=1> solid_vector[numTaxa]; // coating parameters
  real<lower=0> phi[numTaxa]; // dispersion parameters
}

transformed parameters {
  real theta[4*numTaxa]; // vector of parameterss
  real y[numTimeSteps,2*numTaxa]; // raw ODE output 
  real output_mat[numTimeSteps,2*numTaxa];
  
  theta[(0*numTaxa+1):(1*numTaxa)] = coating_vector;
  theta[(1*numTaxa+1):(2*numTaxa)] = O2_vector;
  theta[(2*numTaxa+1):(3*numTaxa)] = HMO_vector;
  theta[(3*numTaxa+1):(4*numTaxa)] = solid_vector;
  
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


