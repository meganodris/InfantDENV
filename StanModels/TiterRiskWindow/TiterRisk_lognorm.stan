//----- Model of lognormal ADE titer-risk scenario -----//

data {
  
  int N[3]; // number of infants
  int Nt; // N time points
  int Ntfit; // N time steps for prediction
  real estHI[N[1],Nt]; // MAb titers
  real estHI2[N[2],Nt]; // MAb titers
  real estNT[N[3],Nt]; // MAb titers
  int cases[13]; // hospital cases by age (months)
  real tfit[Ntfit]; // time points for fits
  real lambda; // monthly FOI

}


parameters {
  
  real mu[3]; 
  real <lower=0> sigma[3];  
  real <lower=0> rhoP;
}


transformed parameters{
  
  matrix[N[1],Nt] pHI;
  matrix[N[1],Nt] probAge;
  matrix[N[2],Nt] pHI2;
  matrix[N[2],Nt] probAge2;
  matrix[N[3],Nt] pNT;
  matrix[N[3],Nt] probAgeNT;
  vector[Nt] probI[3];
  real B = rhoP*1e5;

  // probability density of titers being within risk window 
  for(i in 1:N[1]) for(t in 1:Nt) pHI[i,t] = exp(lognormal_lpdf(estHI[i,t] | mu[1], sigma[1]))/exp(lognormal_lpdf(exp(mu[1]-sigma[1]^2) | mu[1], sigma[1]));
  for(i in 1:N[2]) for(t in 1:Nt) pHI2[i,t] = exp(lognormal_lpdf(estHI2[i,t] | mu[2], sigma[2]))/exp(lognormal_lpdf(exp(mu[2]-sigma[2]^2) | mu[2], sigma[2]));
  for(i in 1:N[3]) for(t in 1:Nt) pNT[i,t] = exp(lognormal_lpdf(estNT[i,t] | mu[3], sigma[3]))/exp(lognormal_lpdf(exp(mu[3]-sigma[3]^2) | mu[3], sigma[3]));


  // conditional prob infection 
  probAge[,1] = pHI[,1]*(1-exp(-lambda));
  probAge2[,1] = pHI2[,1]*(1-exp(-lambda));
  probAgeNT[,1] = pNT[,1]*(1-exp(-lambda));
  for(t in 2:Nt) probAge[,t] = pHI[,t]*((1-exp(-lambda*t))-(1-exp(-lambda*(t-1))));
  for(t in 2:Nt) probAge2[,t] = pHI2[,t]*((1-exp(-lambda*t))-(1-exp(-lambda*(t-1))));
  for(t in 2:Nt) probAgeNT[,t] = pNT[,t]*((1-exp(-lambda*t))-(1-exp(-lambda*(t-1))));
 
  // sum over infants by age
  for(t in 1:Nt) probI[1,t] = mean(probAge[,t]);
  for(t in 1:Nt) probI[2,t] = mean(probAge2[,t]);
  for(t in 1:Nt) probI[3,t] = mean(probAgeNT[,t]);

}


model {
  
  // priors
  mu ~ normal(10,5);
  sigma ~ normal(0,1);
  rhoP ~ normal(0.4,0.5);

  // likelihood
  for(i in 1:3) cases ~ poisson(probI[i,]*B+0.001);

}

generated quantities {
  
  int predCases[3,Nt];
  real predRisk[3,Ntfit];
  real log_lik[3,Nt];

  for(i in 1:3) for(t in 1:Nt) log_lik[i,t] = poisson_lpmf(cases[t] | probI[i,t]*B+0.001);
  for(i in 1:3) predCases[i,] = poisson_rng(probI[i,]*B+0.001);
  for(i in 1:3) for(t in 1:Ntfit) predRisk[i,t] = exp(lognormal_lpdf(tfit[t]| mu[i], sigma[i]))/exp(lognormal_lpdf(exp(mu[i]-sigma[i]^2) | mu[i], sigma[i]));


}

