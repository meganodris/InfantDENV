#----- Run maternal DENV anitbody model -----#
rm(list=ls())
library(rstan)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(loo)

# source functions
setwd('C:/Users/Megan/OneDrive - University of Cambridge/MaternalDENV/CleanCodeRevisions')
source('FunctionsForDecayModel.R')

# read in data
setwd('C:/Users/Megan/OneDrive - University of Cambridge/MaternalDENV/CleanCodeRevisions/Data/rmMean')
data <- readRDS('modelInputs_meanrm.RDS')
df <- read.csv('titerData_meanrm.csv')

#--- list of inputs
dta <- list()
dta$N_hi1 <- dim(data$HItiters1)[2]
dta$N_hi2 <- dim(data$HItiters2)[2]
dta$N_prnt <- dim(data$NTtiters)[2]
dta$Nt <- dim(data$HItiters1)[3]
dta$maxHI <- c(max(data$HItiters1, na.rm=T),max(data$HItiters2,na.rm=T))
dta$obsHI1 <- 1 + log(data$HItiters1/10)/log(2)
dta$obsHI2 <- 1 + log(data$HItiters2/10)/log(2)
dta$obsNT <- 1 + log(data$NTtiters/10)/log(2)
dta$ageHI1 <- data$ageHI1
dta$ageHI2 <- data$ageHI2
dta$ageNT <- data$ageNT
dta$HIind1 <- data$HIind1
dta$HIind2 <- data$HIind2
dta$NTind <- data$NTind
dta$tfit <- seq(0,15,0.1)
dta$Ntfit <- length(dta$tfit) 
dta$Ndata <- sum(dta$HIind1) + sum(dta$HIind2) + sum(dta$NTind)


# replace NAs (will be ignored in model)
dta$obsHI1[which(is.na(dta$obsHI1))] <- 99999
dta$obsHI2[which(is.na(dta$obsHI2))] <- 99999
dta$obsNT[which(is.na(dta$obsNT))] <- 99999
dta$ageHI1[which(is.na(dta$ageHI1))] <- 24
dta$ageHI2[which(is.na(dta$ageHI2))] <- 24
dta$ageNT[which(is.na(dta$ageNT))] <- 24

# fit model
check_cmdstan_toolchain(fix=T)
set_cmdstan_path('C:/Users/Megan/Documents/.cmdstanr/cmdstan-2.26.0')
setwd('C:/Users/Megan/Documents/GitHub/InfantDENV/StanModels')
mod <- cmdstan_model('MAbDecay_cohort.stan', pedantic=T)
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


# discretized fits
amin <- c(-1,0.5,3.1,6.1,9.1) # min & max of age groups
amax <- c(0.5,3.1,6.1,9.1,15)
iter <- dim(chains$sigma)[1]
dNT <- get_discrNT(dta, chains, iter)
dHI <- get_discrHI(dta, chains, iter, amin, amax)

# fit plots
fitp <- plot_fits(df, dta, chains, dNT, dHI)
fitp$sero

# relative serotype proportions 
serod <- serotype_diff(dta, chains)
serod
rels <- ggplot(serod$reldiff, aes(sero, med,col=est))+ 
  geom_point(position=position_dodge(width=0.5), size=2.5)+
  geom_linerange(aes(ymin=ciL,ymax=ciU,col=est), position=position_dodge(width=0.5))+
  facet_wrap(~cohort)+ theme_minimal()+ labs(col='')+
  theme(text=element_text(size=22), axis.text.x=element_text(angle=60, hjust=1),
        legend.position=c(0.55,0.7), axis.title.x=element_blank())+
  ylab('Relative mean titer')+ 
  scale_color_manual(values=c('grey30','mediumvioletred','dodgerblue2'))
plot_grid(fitp$bprnt, fitp$bhi, fitp$kpp, rels)

# output fit plot
setwd('C:/Users/Megan/OneDrive - University of Cambridge/MaternalDENV/CleanCodeRevisions/Results/Cohort_meanrm')
png(filename='CohortFits_ALL.png', width=40, height=30, res=400, units='cm')
plot_grid(fitp$bprnt, fitp$bhi, fitp$kpp, rels, labels=c('A','B','C','D'), label_size=22)
dev.off()
png(filename='CohortFits_Sero.png', width=40, height=30, res=400, units='cm')
plot(fitp$sero)
dev.off()
png(filename='ObsVPred.png', width=40, height=30, res=400, units='cm')
plot_grid(plotlist=OvP)
dev.off()

setwd('C:/Users/Megan/OneDrive - University of Cambridge/MaternalDENV/CleanCodeRevisions/Figures')
pdf(file='Fig2.pdf', width=13, height=9)
plot_grid(fitp$bprnt, fitp$bhi, fitp$kpp, rels, labels=c('A','B','C','D'), label_size=22)
dev.off()

# param estimates
estimates <- data.frame(par=c('gamma_KPP','gamma_BK_HI','gamma_BK_PRNT','sigma_KPP','sigma_BKHI','sigma_PRNT','half_life_KPP',
                              'half_life_BK_HI','half_life_BK_PRNT'),med=NA,ciL=NA,ciU=NA)
estimates[1,2:4] <- quantile(chains$gamma[,1], c(0.5,0.025,0.975))
estimates[2,2:4] <- quantile(chains$gamma[,2], c(0.5,0.025,0.975))
estimates[3,2:4] <- quantile(chains$gamma[,3], c(0.5,0.025,0.975))
estimates[4,2:4] <- quantile(chains$sigma[,1], c(0.5,0.025,0.975))
estimates[5,2:4] <- quantile(chains$sigma[,2], c(0.5,0.025,0.975))
estimates[6,2:4] <- quantile(chains$sigma[,3], c(0.5,0.025,0.975))
estimates[7,2:4] <- quantile(chains$half_life[,1], c(0.5,0.025,0.975))
estimates[8,2:4] <- quantile(chains$half_life[,2], c(0.5,0.025,0.975))
estimates[9,2:4] <- quantile(chains$half_life[,3], c(0.5,0.025,0.975))
estimates[10:12,2:4] <- estimates[7:9,2:4]*30.4167
estimates$par[10:12] <- c('half_life_KPP_days','half_life_BK_HI_days','half_life_BK_PRNT_days')
estimates

# output loo and waic
waic <- waic(log_lik_1)
ww <- data.frame(v=waic[3:5])
loo <- data.frame(v=loo_1[5:7])
write.csv(ww, 'WAIC.csv', row.names=F)
write.csv(loo, 'LOO.csv', row.names=F)
write.csv(estimates, 'Pars.csv', row.names=F)
saveRDS(chains, 'Chains.RDS')
saveRDS(dta, 'Data.RDS')

#------ plot model fit with naive model
setwd('C:/Users/Megan/OneDrive - University of Cambridge/MaternalDENV/CleanCodeRevisions/Results/CohortNaive_meanrm')
chainsNv <- readRDS('Chains.RDS')
nvfit <- data.frame(age=NA, cohort=NA, med=NA, ciL=NA, ciU=NA)
ind <- 0
for(st in 1:3) for(t in 1:dta$Ntfit){
  nvfit[ind+1,] <- c(dta$tfit[t],st,quantile(chainsNv$decay_fit[,st,t], c(0.5,0.025,0.975)))
  ind <- ind+1
}
nvfit$cohort <- factor(nvfit$cohort, labels=c('KPP HI','Bangkok HI','Bangkok PRNT'))
head(nvfit)

# population mean observed
popm <- PopMeanObs(obs)
popm$meanObs$cohort <- factor(popm$meanObs$cohort, levels=c('KPP','BKHI','BKNT'), labels=c('KPP HI','Bangkok HI','Bangkok PRNT'))

# extract mean decay fits
fit <- data.frame(age=NA, cohort=NA, med=NA, ciL=NA, ciU=NA)
ind <- 0
for(st in 1:3) for(t in 1:dta$Ntfit){
  fit[ind+1,] <- c(dta$tfit[t],st,quantile(chains$decay_fit[,st,t], c(0.5,0.025,0.975)))
  ind <- ind+1
}
fit$cohort <- factor(fit$cohort, labels=c('KPP HI','Bangkok HI','Bangkok PRNT'))

# mean titer fits
prntfit <- ggplot(popm$meanObs[popm$meanObs$cohort=='Bangkok PRNT',], aes(age,mean))+ 
  geom_rect(aes(ymin=0.01, ymax=10, xmin=0, xmax=15), fill='grey90', alpha=0.1)+
  theme_minimal()+ xlab('Age (months)')+ ylab('PRNT titer')+
  geom_line(data=fit[fit$cohort=='Bangkok PRNT',], aes(age, med),col='springgreen3')+
  geom_ribbon(data=fit[fit$cohort=='Bangkok PRNT',], aes(age,y=med,ymin=ciL,ymax=ciU),fill='springgreen3', alpha=0.3)+
  geom_point()+ geom_linerange(aes(ymin=ciL,ymax=ciU))+ 
  geom_point(data=dNT$dfitmean, aes(age+0.3, fit), col='springgreen3')+
  geom_linerange(data=dNT$dfitmean, aes(age+0.3, y=fit, ymin=ciL, ymax=ciU), col='springgreen3')+
  theme(text=element_text(size=22), legend.position='none')+
  scale_x_continuous(breaks=c(0,3,6,9,12))+ 
  geom_line(data=nvfit[nvfit$cohort=='Bangkok PRNT',], aes(age,med), col='grey50', linetype='dashed')+
  geom_ribbon(data=nvfit[nvfit$cohort=='Bangkok PRNT',], aes(age,y=med,ymin=ciL,ymax=ciU), fill='grey70', alpha=0.3)
prntfit
prntfit_log <- prntfit+ scale_y_continuous(trans='log', breaks=c(0.1,1,10,100,1000), labels=comma)+
  theme_bw()+ theme(text=element_text(size=22))
prntf <- ggdraw() +
  draw_plot(prntfit) + draw_plot(prntfit_log, x=.4, y=.4, width=.5, height=.5)
prntf

# bangkok hi
bhifit <- ggplot(popm$meanObs[popm$meanObs$cohort=='Bangkok HI',], aes(age,mean))+ 
  geom_rect(aes(ymin=0.01, ymax=10, xmin=0, xmax=15), fill='grey90', alpha=0.1)+
  theme_minimal()+ xlab('Age (months)')+ ylab('HI titer')+
  geom_line(data=fit[fit$cohort=='Bangkok HI',], aes(age, med),col='springgreen3')+
  geom_ribbon(data=fit[fit$cohort=='Bangkok HI',], aes(age,y=med,ymin=ciL,ymax=ciU),fill='springgreen3', alpha=0.3)+
  geom_point()+ geom_linerange(aes(ymin=ciL,ymax=ciU))+ 
  geom_point(data=dHI$dfitmean2, aes(age+0.3, fit), col='springgreen3')+
  geom_linerange(data=dHI$dfitmean2, aes(age+0.3, y=fit, ymin=ciL, ymax=ciU), col='springgreen3')+
  theme(text=element_text(size=22), legend.position='none')+
  scale_x_continuous(breaks=c(0,3,6,9,12))+
  geom_line(data=nvfit[nvfit$cohort=='Bangkok HI',], aes(age,med), col='grey50', linetype='dashed')+
  geom_ribbon(data=nvfit[nvfit$cohort=='Bangkok HI',], aes(age,y=med,ymin=ciL,ymax=ciU), fill='grey70', alpha=0.3)
bhifit_log <- bhifit+ scale_y_continuous(trans='log', breaks=c(0.1,1,10,100,1000), labels=comma)+
  theme_bw()+ theme(text=element_text(size=22))
bhif <- ggdraw() +
  draw_plot(bhifit) + draw_plot(bhifit_log, x=.4, y=.4, width=.5, height=.5)
bhif
# kpp hi
khifit <- ggplot(popm$meanObs[popm$meanObs$cohort=='KPP HI',], aes(age,mean))+ 
  geom_rect(aes(ymin=0.01, ymax=10, xmin=0, xmax=15), fill='grey90', alpha=0.1)+
  theme_minimal()+ xlab('Age (months)')+ ylab('HI titer')+
  geom_line(data=fit[fit$cohort=='KPP HI',], aes(age, med),col='springgreen3')+
  geom_ribbon(data=fit[fit$cohort=='KPP HI',], aes(age,y=med,ymin=ciL,ymax=ciU),fill='springgreen3', alpha=0.3)+
  geom_point()+ geom_linerange(aes(ymin=ciL,ymax=ciU))+ 
  geom_point(data=dHI$dfitmean1, aes(age+0.3, fit), col='springgreen3')+
  geom_linerange(data=dHI$dfitmean1, aes(age+0.3, y=fit, ymin=ciL, ymax=ciU), col='springgreen3')+
  theme(text=element_text(size=22), legend.position='none')+
  scale_x_continuous(breaks=c(0,3,6,9,12))+
  geom_line(data=nvfit[nvfit$cohort=='KPP HI',], aes(age,med), col='grey50', linetype='dashed')+
  geom_ribbon(data=nvfit[nvfit$cohort=='KPP HI',], aes(age,y=med,ymin=ciL,ymax=ciU), fill='grey70', alpha=0.3)
khifit_log <- khifit+ scale_y_continuous(trans='log', breaks=c(0.1,1,10,100,1000), labels=comma)+
  theme_bw()+ theme(text=element_text(size=22))
khif <- ggdraw() +
  draw_plot(khifit) + draw_plot(khifit_log, x=.4, y=.4, width=.5, height=.5)
khif

setwd('C:/Users/Megan/OneDrive - University of Cambridge/MaternalDENV/CleanCodeRevisions/Results/Cohort_meanrm')
png(filename='CohortFits_ALLnaive.png', width=40, height=30, res=400, units='cm')
plot_grid(prntf, bhif, khif, rels, labels=c('A','B','C','D'), label_size=22)
dev.off()
