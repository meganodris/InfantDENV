#----- Run maternal DENV anitbody decay models -----#
rm(list=ls())
library(rstan)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(loo)

# source functions
source('TiterDecayModelFunctions.R')

# read in data
data <- readRDS('Data/DecayModelInputs.RDS')

# please note that the age of infants has been rounded and grouped for anonymization.
# therefore, the model results may differ slightly from the published estimates. 


#--- list of inputs for model
dta <- list()
dta$N_hi1 <- dim(data$HItiters1)[2] # N infants [KPP]
dta$N_hi2 <- dim(data$HItiters2)[2] # N infants [Bangkok]
dta$N_prnt <- dim(data$NTtiters)[2] # N infants [Bangkok]
dta$Nt <- dim(data$HItiters1)[3] # N samples per person
dta$maxHI <- c(max(data$HItiters1, na.rm=T),max(data$HItiters2,na.rm=T)) # max HI dilutions
dta$obsHI1 <- 1 + log(data$HItiters1/10)/log(2) # HI titer data [KPP]
dta$obsHI2 <- 1 + log(data$HItiters2/10)/log(2) # HI titer data [Bangkok]
dta$obsNT <- 1 + log(data$NTtiters/10)/log(2) # PRNT titer data [Bangkok]
dta$ageHI1 <- data$ageHI1 # age at sampling [KPP HI]
dta$ageHI2 <- data$ageHI2 # age at sampling [Bangkok HI]
dta$ageNT <- data$ageNT # age at sampling [Bangkok PRNT]
dta$HIind1 <- data$HIind1 # indexes for if titer measurement was available or not
dta$HIind2 <- data$HIind2 # indexes for if titer measurement was available or not
dta$NTind <- data$NTind # indexes for if titer measurement was available or not
dta$tfit <- seq(0,15,0.1) # ages for model fits
dta$Ntfit <- length(dta$tfit) # N timepoints for fitting
dta$Ndata <- sum(dta$HIind1) + sum(dta$HIind2) + sum(dta$NTind) # total N infants for LOO calculation


# replace NAs (will be ignored in model)
dta$obsHI1[which(is.na(dta$obsHI1))] <- 99999
dta$obsHI2[which(is.na(dta$obsHI2))] <- 99999
dta$obsNT[which(is.na(dta$obsNT))] <- 99999
dta$ageHI1[which(is.na(dta$ageHI1))] <- 24
dta$ageHI2[which(is.na(dta$ageHI2))] <- 24
dta$ageNT[which(is.na(dta$ageNT))] <- 24


#--- fit model
check_cmdstan_toolchain(fix=T)
set_cmdstan_path('C:/Users/Megan/Documents/.cmdstanr/cmdstan-2.26.0')
mod <- cmdstan_model('StanModels/MAbDecay/MAbDecay_cohort.stan', pedantic=T)
fit <- mod$sample(data=dta, chains=4, parallel_chains=4, iter_sampling=5000, refresh=100, iter_warmup=1000)
stanfit <- rstan::read_stan_csv(fit$output_files())


# check convergence
chains <- rstan::extract(stanfit)
traceplot(stanfit, pars=c('gamma','sigma','half_life'))


# extract pointwise log-likelihood
log_lik_1 <- extract_log_lik(stanfit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2) 
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
waic(log_lik_1)

# observed vs predicted titers
OvP <- plot_ObsVsPred(chains, dta)
plot_grid(plotlist=OvP)

