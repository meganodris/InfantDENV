//--- Model of infant maternal DENV antibody decay ---//


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
  vector<lower=0>[3] gamma1; // estimated decay rate before change point
  vector<lower=0>[3] gamma2; // estimated decay rate after change point
  real <lower=0> sigma[3]; // error terms
  real <lower=1,upper=10> agebk; // age of decay change

}


transformed parameters{
  
  real pHI1[4,N_hi1,Nt]; // estimated HI titers on log2 scale
  real pHI2[4,N_hi2,Nt]; // estimated HI titers on log2 scale
  real pNT[4,N_prnt,Nt]; // estimated PRNT titers on log2 scale
  real likHI1[4,N_hi1,Nt]; // likelihood values
  real likHI2[4,N_hi2,Nt]; // likelihood values
  real likNT[4,N_prnt,Nt]; // likelihood values
  vector[3] gamma[2]; // rate of decay
  
  gamma[1] = gamma1;
  gamma[2] = gamma2;

  
  // model decay of estimated titers
  for(s in 1:4) for(i in 1:N_hi1) for(t in 1:Nt){
    if(ageHI1[i,t]<agebk) pHI1[s,i,t] = 1 + log(exp((trueT0HI1[s,i]-1)*log(2)+log(10))*exp(-gamma[1,1]*ageHI1[i,t])/10)/log(2);
    else pHI1[s,i,t] = 1 + log(exp((trueT0HI1[s,i]-1)*log(2)+log(10))*exp(-gamma[1,1]*agebk)*exp(-gamma[2,1]*(ageHI1[i,t]-agebk))/10)/log(2);
  }
  for(s in 1:4) for(i in 1:N_hi2) for(t in 1:Nt){
    if(ageHI2[i,t]<agebk) pHI2[s,i,t] = 1 + log(exp((trueT0HI2[s,i]-1)*log(2)+log(10))*exp(-gamma[1,2]*ageHI2[i,t])/10)/log(2);
    else pHI2[s,i,t] = 1 + log(exp((trueT0HI2[s,i]-1)*log(2)+log(10))*exp(-gamma[1,2]*agebk)*exp(-gamma[2,2]*(ageHI2[i,t]-agebk))/10)/log(2);
  }
  for(s in 1:4) for(i in 1:N_prnt) for(t in 1:Nt){
    if(ageNT[i,t]<agebk) pNT[s,i,t] = 1 + log(exp((trueT0NT[s,i]-1)*log(2)+log(10))*exp(-gamma[1,3]*ageNT[i,t])/10)/log(2);
    else pNT[s,i,t] = 1 + log(exp((trueT0NT[s,i]-1)*log(2)+log(10))*exp(-gamma[1,3]*agebk)*exp(-gamma[2,3]*(ageNT[i,t]-agebk))/10)/log(2);
  }
  
  
  //--- HI likelihood cohort 1 (discretized scale titers) ---//
  for(s in 1:4) for(i in 1:N_hi1) for(t in 1:Nt){
    
   if(HIind1[s,i,t]==1){
     
     if(obsHI1[s,i,t]<1){ // for titers under limit of detection
     
       likHI1[s,i,t] = normal_lcdf(1-pHI1[s,i,t] | 0, sigma[1]);
       
     }else if(obsHI1[s,i,t]>maxHI[1]){ // for titers above highest 
     
       likHI1[s,i,t] = log(1-normal_cdf(maxHI[1]-pHI1[s,i,t], 0, sigma[1]));
       
     }else{
       
       likHI1[s,i,t] = log(normal_cdf(obsHI1[s,i,t]+1-pHI1[s,i,t], 0, sigma[1]) - normal_cdf(obsHI1[s,i,t]-pHI1[s,i,t], 0, sigma[1]));
       
     }
    
     if(likHI1[s,i,t]<log(1e-8)) likHI1[s,i,t] = log(1e-8);
     
   }else likHI1[s,i,t] = -999; // no data: will not contribute to likelihood 
  
  }
  
  
  //--- HI likelihood cohort 2 (discretized scale titers) ---//
  for(s in 1:4) for(i in 1:N_hi2) for(t in 1:Nt){
    
   if(HIind2[s,i,t]==1){
     
     if(obsHI2[s,i,t]<1){ // for titers under limit of detection
     
       likHI2[s,i,t] = normal_lcdf(1-pHI2[s,i,t] | 0, sigma[2]);
       
     }else if(obsHI2[s,i,t]>maxHI[2]){ // for titers above highest 
     
       likHI2[s,i,t] = log(1-normal_cdf(maxHI[2]-pHI2[s,i,t], 0, sigma[2]));
       
     }else{
       
       likHI2[s,i,t] = log(normal_cdf(obsHI2[s,i,t]+1-pHI2[s,i,t], 0, sigma[2]) - normal_cdf(obsHI2[s,i,t]-pHI2[s,i,t], 0, sigma[2]));
       
     }
    
     if(likHI2[s,i,t]<log(1e-8)) likHI2[s,i,t] = log(1e-8);
     
   }else likHI2[s,i,t] = -999; // no data: will not contribute to likelihood 
  
  }
  
  
  //--- PRNT likelihood (continuous-scale titers) ---//
  for(s in 1:4) for(i in 1:N_prnt) for(t in 1:Nt){
    
    if(NTind[s,i,t]==1){
      
      if(obsNT[s,i,t]<1){ // for titers under limit of detection
      
        likNT[s,i,t] = normal_lcdf(1-pNT[s,i,t] | 0, sigma[3]);

      }else{
      
        likNT[s,i,t] = normal_lpdf(pNT[s,i,t] | obsNT[s,i,t], sigma[3]);
       
      }
    
      if(likNT[s,i,t]<log(1e-8)) likNT[s,i,t] = log(1e-8);
      
    }else likNT[s,i,t] = -999; // no data: will not contribute to likelihood 

  }

}


model {
  
  //--- priors ---//
  sigma ~ normal(1,1);
  agebk ~ normal(3,0.5);
  gamma1 ~ normal(0,1);
  gamma2 ~ normal(0,1);
  for(s in 1:4) to_vector(trueT0HI1[s,]) ~ uniform(0,20);
  for(s in 1:4) to_vector(trueT0HI2[s,]) ~ uniform(0,20);
  for(s in 1:4) to_vector(trueT0NT[s,]) ~ uniform(0,20);

  
  //--- likelihood ---//
  for(i in 1:N_hi1) for(t in 1:Nt) for(s in 1:4){
    
    if(HIind1[s,i,t]==1) target += likHI1[s,i,t]; // HI cohort 1
    
  }
  for(i in 1:N_hi2) for(t in 1:Nt) for(s in 1:4){
    
    if(HIind2[s,i,t]==1) target += likHI2[s,i,t]; // HI cohort 2
    
  }
  for(i in 1:N_prnt) for(t in 1:Nt) for(s in 1:4){
    
    if(NTind[s,i,t]==1) target += likNT[s,i,t]; // PRNT cohort 1
    
  }

  
}


generated quantities {
  
  real log_lik[Ndata];
  real half_life[2,3];
  real meanT0[3];
  real meanT0s[3,4];
  real decay_fit[3,Ntfit];
  real decay_fit_sero[3,4,Ntfit];
  {
    int ind = 1;

  
  // likelihood
  for(i in 1:N_hi1) for(t in 1:Nt) for(s in 1:4) if(HIind1[s,i,t]==1){
     log_lik[ind] = likHI1[s,i,t];
     ind = ind+1;
  }
  for(i in 1:N_hi2) for(t in 1:Nt) for(s in 1:4) if(HIind2[s,i,t]==1){
    log_lik[ind] = likHI2[s,i,t];
    ind = ind+1;
  }
  for(i in 1:N_prnt) for(t in 1:Nt) for(s in 1:4) if(NTind[s,i,t]==1){
    log_lik[ind] = likNT[s,i,t];
    ind=ind+1;
  }

  // half-lives
  for(i in 1:3) for(b in 1:2) half_life[b,i] = log(2)/gamma[b,i];
  
  // mean cord titers
  meanT0[1] = mean(exp((trueT0HI1-1)*log(2)+log(10)));
  meanT0[2] = mean(exp((trueT0HI2-1)*log(2)+log(10)));
  meanT0[3] = mean(exp((trueT0NT-1)*log(2)+log(10)));
  for(s in 1:4) meanT0s[1,s] = mean(exp((trueT0HI1[s,]-1)*log(2)+log(10)));
  for(s in 1:4) meanT0s[2,s] = mean(exp((trueT0HI2[s,]-1)*log(2)+log(10)));
  for(s in 1:4) meanT0s[3,s] = mean(exp((trueT0NT[s,]-1)*log(2)+log(10)));
  
  // decay fit
  for(t in 1:Ntfit){
    if(tfit[t]<agebk){
      decay_fit[1,t] = meanT0[1]*exp(-gamma[1,1]*tfit[t]);
      decay_fit[2,t] = meanT0[2]*exp(-gamma[1,2]*tfit[t]);
      decay_fit[3,t] = meanT0[3]*exp(-gamma[1,3]*tfit[t]);
    }else{
      decay_fit[1,t] = meanT0[1]*exp(-gamma[1,1]*agebk)*exp(-gamma[2,1]*(tfit[t]-agebk));
      decay_fit[2,t] = meanT0[2]*exp(-gamma[1,2]*agebk)*exp(-gamma[2,2]*(tfit[t]-agebk));
      decay_fit[3,t] = meanT0[3]*exp(-gamma[1,3]*agebk)*exp(-gamma[2,3]*(tfit[t]-agebk));
    }
  }
  
  // serotype decay fit
  for(s in 1:4) for(t in 1:Ntfit){
    if(tfit[t]<agebk){
      decay_fit_sero[1,s,t] = meanT0s[1,s]*exp(-gamma[1,1]*tfit[t]);
      decay_fit_sero[2,s,t] = meanT0s[2,s]*exp(-gamma[1,2]*tfit[t]);
      decay_fit_sero[3,s,t] = meanT0s[3,s]*exp(-gamma[1,3]*tfit[t]);
    }else{
      decay_fit_sero[1,s,t] = meanT0s[1,s]*exp(-gamma[1,1]*agebk)*exp(-gamma[2,1]*(tfit[t]-agebk));
      decay_fit_sero[2,s,t] = meanT0s[2,s]*exp(-gamma[1,2]*agebk)*exp(-gamma[2,2]*(tfit[t]-agebk));
      decay_fit_sero[3,s,t] = meanT0s[3,s]*exp(-gamma[1,3]*agebk)*exp(-gamma[2,3]*(tfit[t]-agebk));
    }
  }

  }
}


