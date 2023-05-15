#----- Code to fit titer-related mechanisms of infant dengue disease -----#
library(cmdstanr)
library(rstan)
library(ggplot2)
library(loo)


# read in age-specific hospital case data
cases <- read.csv('~/data/InfantCases.csv')

# read in median estimates of infant maternal antibodies by age
MAbs <- readRDS('~data/MAbEstimates.RDS')

# data for model fitting
data <- list()
data$N <- c(nrow(MAbs$estHI_KPP),nrow(MAbs$estHI_BK),nrow(MAbs$estNT_BK)) # N infants per study
data$Nt <- 13 # N ages 
data$estHI <- MAbs$estHI_KPP # KPP HI titer estimates
data$estHI2 <- MAbs$estHI_BK # Bangkok HI titer estimates
data$estNT <- MAbs$estNT_BK # Bangkok PRNT titer estimates
data$cases <- cases$cases # age-specific case data
data$tfit <- seq(0,150,0.5) # titer values for risk prediction
data$Ntfit <- length(data$tfit) # N time points for prediction
data$lambda <- 0.12/12 # assumed monthly FOI
data$delta <- 0.15 # age-specific vulnerability to disease given infection


# model options
mod1 <- cmdstan_model('TiterRisk_lognorm.stan', pedantic=T)
mod2 <- cmdstan_model('TiterRisk_lognorm_age.stan', pedantic=T)
mod3 <- cmdstan_model('TiterRisk_norm.stan', pedantic=T)
mod4 <- cmdstan_model('TiterRisk_norm_age.stan', pedantic=T)
mod5 <- cmdstan_model('TiterRisk_threshold.stan', pedantic=T)
mod6 <- cmdstan_model('TiterRisk_threshold_age.stan', pedantic=T)


# fit model1 for example
fit <- mod1$sample(data=data, chains=4, parallel_chains=4, iter_sampling=4000, refresh=100, iter_warmup=2000)

# Check convergence
stanfit <- rstan::read_stan_csv(fit$output_files())
chains <- rstan::extract(stanfit)
traceplot(stanfit, pars=c('mu','sigma','B','lp__'))

# WAIC 
loglik1 <- extract_log_lik(stanfit, merge_chains=F)
waic(loglik1)

# predicted cases by age
cohorts <- c('KPP [HI]','Bangkok [HI]','Bangkok [PRNT]')
fit <- data.frame(cohort=NA,mo=seq(0,12),obs=data$cases,fit=NA,ciL=NA,ciU=NA)
fits <- list()
for(c in 1:3){
  fits[[c]] <- fit
  fits[[c]]$cohort <- cohorts[c]
  for(t in 1:13) fits[[c]][t,4:6] <- quantile(chains$predCases[,c,t], c(0.5,0.025,0.975))
} 
fits <- do.call('rbind', fits)
ggplot(fits, aes(mo,obs))+ geom_point()+ theme_minimal()+
  geom_line(aes(mo,fit,col=cohort))+ xlab('Age (months)')+ ylab('Hospitalised cases')+
  geom_ribbon(aes(mo,ymin=ciL,ymax=ciU,fill=cohort), alpha=0.2)+
  theme(text=element_text(size=22))

# risk titer distributions
rTs <- list()
rT <- data.frame(titer=data$tfit,p=NA,ciL=NA,ciU=NA)
for(c in 1:3){
  for(i in 1:data$Ntfit) rT[i,2:4] <- quantile(chains$predRisk[,c,i], c(0.5,0.025,0.975))
  rTs[[c]] <- rT
  rTs[[c]]$cohort <- cohorts[c]
}
rTs <- do.call('rbind',rTs)
ggplot(rTs, aes(titer,p))+ geom_line(aes(col=cohort))+ theme_minimal()+ 
  geom_ribbon(aes(ymin=ciL,ymax=ciU,fill=cohort),alpha=0.3)+ scale_x_continuous(trans='log')



