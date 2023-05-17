//--- Model of infant maternal DENV antibody decay ---//
// naive model not accounting for discretized measurememt processes or lower assay limits

data {
  
  int N_hi1; // number of infants HI (Bangkok)
  int N_hi2; // number of infants HI (KPP)
  int N_prnt; // number of infants PRNT (Bangkok)
  int Nt; // number of time points
  matrix[N_hi1,Nt] obsHI1[4]; // HI titers (Bangkok)
  matrix[N_hi2,Nt] obsHI2[4]; // HI titers (KPP)
  matrix[N_prnt,Nt] obsNT[4]; // PRNT titers (Bangkok)
  matrix[N_hi1,Nt] HIind1[4]; // HI index (0/1 to indicate any NAs)
  matrix[N_hi2,Nt] HIind2[4]; // HI index (0/1 to indicate any NAs)
  matrix[N_prnt,Nt] NTind[4]; // PRNT index (0/1 to indicate any NAs)
  matrix[N_hi1,Nt] ageHI1; // age of infants HI1
  matrix[N_hi2,Nt] ageHI2; // age of infants HI2
  matrix[N_prnt,Nt] ageNT; // age of infants PRNT
  int maxHI[2]; // maximum titer dilution of HI assays
  int Ntfit; // N time steps for prediction
  vector[Ntfit] tfit; // time points for prediction
  int Ndata; // number of data points

}


parameters {
  
  matrix<lower=0>[4,N_hi1] trueT0HI1; // estimated titers at t0 
  matrix<lower=0>[4,N_hi2] trueT0HI2; // estimated titers at t0 
  matrix<lower=0>[4,N_prnt] trueT0NT; // estimated titers at t0
  real gamma[3]; // rate of decay
  real <lower=0> sigma[3]; // error terms

}


transformed parameters{
  
  real pHI1[4,N_hi1,Nt]; // estimated HI titers on log2 scale
  real pHI2[4,N_hi2,Nt]; // estimated HI titers on log2 scale
  real pNT[4,N_prnt,Nt]; // estimated PRNT titers on log2 scale

  
  // model decay of estimated titers
  for(s in 1:4) for(i in 1:N_hi1) for(t in 1:Nt) pHI1[s,i,t] = 1 + log(exp((trueT0HI1[s,i]-1)*log(2)+log(10))*exp(-gamma[1]*ageHI1[i,t])/10)/log(2);
  for(s in 1:4) for(i in 1:N_hi2) for(t in 1:Nt) pHI2[s,i,t] = 1 + log(exp((trueT0HI2[s,i]-1)*log(2)+log(10))*exp(-gamma[2]*ageHI2[i,t])/10)/log(2);
  for(s in 1:4) for(i in 1:N_prnt) for(t in 1:Nt) pNT[s,i,t] = 1 + log(exp((trueT0NT[s,i]-1)*log(2)+log(10))*exp(-gamma[3]*ageNT[i,t])/10)/log(2);

}


model {
  
  //--- priors ---//
  sigma ~ normal(1,1);
  gamma ~ normal(0,1);
  for(s in 1:4) to_vector(trueT0HI1[s,]) ~ uniform(0,20);
  for(s in 1:4) to_vector(trueT0HI2[s,]) ~ uniform(0,20);
  for(s in 1:4) to_vector(trueT0NT[s,]) ~ uniform(0,20);

  
  //--- likelihood ---//
  for(i in 1:N_hi1) for(t in 1:Nt) for(s in 1:4){
    
    if(HIind1[s,i,t]==1) target += normal_lpdf(pHI1[s,i,t] | obsHI1[s,i,t], sigma[1]); // HI cohort 1
    
  }
  for(i in 1:N_hi2) for(t in 1:Nt) for(s in 1:4){
    
    if(HIind2[s,i,t]==1) target += normal_lpdf(pHI2[s,i,t] | obsHI2[s,i,t], sigma[2]); // HI cohort 2
    
  }
  for(i in 1:N_prnt) for(t in 1:Nt) for(s in 1:4){
    
    if(NTind[s,i,t]==1) target += normal_lpdf(pNT[s,i,t] | obsNT[s,i,t], sigma[3]); // PRNT cohort 1
    
  }

  
}


generated quantities {
  
  real log_lik[Ndata];
  real half_life[3];
  real meanT0[3];
  real meanT0s[3,4];
  real decay_fit[3,Ntfit];
  real decay_fit_sero[3,4,Ntfit];
  {
    int ind = 1;

  
  // likelihood
  for(i in 1:N_hi1) for(t in 1:Nt) for(s in 1:4) if(HIind1[s,i,t]==1){
     log_lik[ind] = normal_lpdf(pHI1[s,i,t] | obsHI1[s,i,t], sigma[1]);
     ind = ind+1;
  }
  for(i in 1:N_hi2) for(t in 1:Nt) for(s in 1:4) if(HIind2[s,i,t]==1){
    log_lik[ind] = normal_lpdf(pHI2[s,i,t] | obsHI2[s,i,t], sigma[2]);
    ind = ind+1;
  }
  for(i in 1:N_prnt) for(t in 1:Nt) for(s in 1:4) if(NTind[s,i,t]==1){
    log_lik[ind] = normal_lpdf(pNT[s,i,t] | obsNT[s,i,t], sigma[3]);
    ind=ind+1;
  }

  // half-lives
  for(i in 1:3) half_life[i] = log(2)/gamma[i];
  
  // mean cord titers
  meanT0[1] = mean(exp((trueT0HI1-1)*log(2)+log(10)));
  meanT0[2] = mean(exp((trueT0HI2-1)*log(2)+log(10)));
  meanT0[3] = mean(exp((trueT0NT-1)*log(2)+log(10)));
  for(s in 1:4) meanT0s[1,s] = mean(exp((trueT0HI1[s,]-1)*log(2)+log(10)));
  for(s in 1:4) meanT0s[2,s] = mean(exp((trueT0HI2[s,]-1)*log(2)+log(10)));
  for(s in 1:4) meanT0s[3,s] = mean(exp((trueT0NT[s,]-1)*log(2)+log(10)));
  
  // decay fit
  for(t in 1:Ntfit) decay_fit[1,t] = meanT0[1]*exp(-gamma[1]*tfit[t]);
  for(t in 1:Ntfit) decay_fit[2,t] = meanT0[2]*exp(-gamma[2]*tfit[t]);
  for(t in 1:Ntfit) decay_fit[3,t] = meanT0[3]*exp(-gamma[3]*tfit[t]);
  
  // serotype decay fit
  for(s in 1:4) for(t in 1:Ntfit) decay_fit_sero[1,s,t] = meanT0s[1,s]*exp(-gamma[1]*tfit[t]);
  for(s in 1:4) for(t in 1:Ntfit) decay_fit_sero[2,s,t] = meanT0s[2,s]*exp(-gamma[2]*tfit[t]);
  for(s in 1:4) for(t in 1:Ntfit) decay_fit_sero[3,s,t] = meanT0s[3,s]*exp(-gamma[3]*tfit[t]);
  
  }
}


