#----- Simulate changes in DENV FOI and infant risk profiles over time -----#
library(ggplot2)
library(truncnorm)
library(cowplot)
library(viridis)

# source simulation functions
source('SimulationFunctions.R')

# assumed force of infection between 1980 and 2020
nYears <- 41
FOI <- seq(0.1,0.01,length.out=nYears)
plot(FOI*4)

# force of infection prior to 1980 (assumed constant)
HFOI <- c(rep(FOI[1],100),FOI)

# read in age-specific birth data (https://population.un.org/wpp/Download/Standard/Fertility/)
setwd('~/data')
fta <- read.csv('UNWPP2022_birthsage_singleyears.csv')
fta <- tidyr::gather(fta, key='age', value='rate', 2:36)
colnames(fta)[1] <- 'year'
fta$age <- as.numeric(substr(fta$age,2,3))
fta <- fta[fta$year>1979,]
ggplot(fta, aes(age,rate*1000,group=year,col=year))+ geom_line()+
  theme_minimal()+ theme(text=element_text(size=22), legend.position=c(0.8,0.8))+ 
  xlab('Age')+ ylab('N Births')+ scale_color_viridis_c()+ labs(col='Year') 


#--- Parameters for ADE lognormal simulation
Yrs <- seq(1980,2020) # years
YrsH <- seq(1880,2020) # 100 years prior
ages <- seq(15,49)
nAgeG <- length(ages)
tt <- seq(0,12,0.1) # infant ages (months) for simulation
HL <- log(2)/(28.1/30.44) # assume PRNT half-life of 28.1 days
logmu <- 3.0 # lognormal mean risk titer
logsd <- 1.4 # lognormal standard deviation 

# run simulation
run <- simTrendsADE(Yrs=Yrs, YrsH=YrsH, FOI=FOI, FOIH=HFOI, FTA=fta, nAgeG=nAgeG,
                    tt=tt, HL=HL, logmu=logmu, logsd=logsd)

# plot mothers immune status
stP <- tidyr::gather(run$statusP, key='status',value='prop',2:4)
stP$status <- factor(stP$status, levels=c('Su','Mo','Mu'), labels=c('Susceptible','Monotypic','Multitypic'))
ggplot(stP, aes(years,prop,fill=status))+ geom_col(alpha=0.7)+
  theme_minimal()+ ylab('Proportion')+ xlab('Year')+ labs(fill='')+ 
  theme(text=element_text(size=22), legend.position=c(0.25,0.5),legend.background=element_rect(), legend.title=element_blank())+
  scale_fill_manual(values=c('purple','orange','seagreen'))

# plot yearly infant disease risk
driski <- list()
drisk <- list()
for(y in 1:nYears) driski[[y]] <- data.frame(age=tt, year=Yrs[y], risk=rowMeans(run$Risk[[y]]))
for(y in 1:nYears) drisk[[y]] <- data.frame(age=tt, year=Yrs[y], risk=rowMeans(run$RiskDI[[y]]))
drisk <- do.call('rbind',drisk)
driski <- do.call('rbind',driski)
ggplot(driski, aes(age,risk,group=year))+ geom_line(aes(col=year))+ xlab('age (months)')+
  scale_color_viridis_c()+ theme_minimal()+ ylab('Relative risk of hospitalised disease given infection')
ggplot(drisk, aes(age,risk,group=year))+ geom_line(aes(col=year))+ xlab('age (months)')+
  scale_color_viridis_c()+ theme_minimal()+ ylab('Risk of hospitalised disease')


#--- Parameters for threshold (non-ADE) simulation
Yrs <- seq(1980,2020) # years
YrsH <- seq(1880,2020) # 100 years prior
ages <- seq(15,49)
nAgeG <- length(ages)
tt <- seq(0,12,0.1) # infant ages (months) for simulation
HL <- log(2)/(28.1/30.44) # assume PRNT half-life of 28.1 days
mu <- 11.8 # mean protective threshold titer
sd <- 0.03 # standard deviation 
delta <- 0.15 # age-specific vulnerability to disease (exponential decay rate)

# run simulation
runT <- simTrendsT(Yrs=Yrs, YrsH=YrsH, FOI=FOI, FOIH=HFOI, FTA=fta, nAgeG=nAgeG,
                  tt=tt, HL=HL, mu=mu, sd=sd, delta=delta)


# plot yearly infant disease risk
driski_T <- list()
drisk_T <- list()
for(y in 1:nYears) driski_T[[y]] <- data.frame(age=tt, year=Yrs[y], risk=rowMeans(runT$Risk[[y]]))
for(y in 1:nYears) drisk_T[[y]] <- data.frame(age=tt, year=Yrs[y], risk=rowMeans(runT$RiskDI[[y]]))
drisk_T <- do.call('rbind',drisk_T)
driski_T <- do.call('rbind',driski_T)
ggplot(driski_T, aes(age,risk,group=year))+ geom_line(aes(col=year))+ xlab('age (months)')+
  scale_color_viridis_c()+ theme_minimal()+ ylab('Relative risk of hospitalised disease given infection')
ggplot(drisk_T, aes(age,risk,group=year))+ geom_line(aes(col=year))+ xlab('age (months)')+
  scale_color_viridis_c()+ theme_minimal()+ ylab('Risk of hospitalised disease')


