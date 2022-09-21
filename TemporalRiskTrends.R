#----- Simulate changes in FOI and infant risk -----#
#------ Thailand 1980-2019 ------#
rm(list=ls())
library(ggplot2)
library(truncnorm)


# assumed force of infection between 1980 and 2019
nYears <- 40
FOIrecent <- seq(0.18,0.06,length.out=nYears)/4
plot(FOIrecent)

# force of infection prior to 1980 (assumed constant)
FOI <- c(rep(FOIrecent[1],100),FOIrecent)
years <- seq(1880,2019)

# read in age-specific fertility rates
setwd('C:/Users/Megan/Documents/GitHub/InfantDENV/data')
fta <- read.csv('fertilityAge.csv')
fta <- tidyr::gather(fta, key='age', value='rate', 2:8)
ggplot(fta, aes(period,rate,group=age,col=age))+ geom_line()+
  theme_minimal()
period <- unique(fta$period)[6:14]
yp <- sort(rep(seq(1,9),5))
yp




nsim <- 100
ages <- c('X15.19','X20.24','X25.29','X30.34','X35.39','X40.44','X45.49')
ageL <- c(15,20,25,30,35,40,45)
ageU <- c(19,24,29,34,39,44,49)
years <- seq(1980,2019)
probSus <- probMonos <- probMultis <- matrix(NA,nsim,nYears)
TiterList <- list()

# simulate
nsim <- 1
for(i in 1:nsim){
  ageB <- list()
  titers <- list()
  for(y in 1:nYears){
    ratesA <- fta$rate[fta$period==period[yp[y]]]
    for(a in 1:7) ageB[[a]] <- runif(floor(ratesA[a]),ageL[a],ageU[a]+1) # draw age of mother for each baby
    agesM <- unlist(ageB)
    pS <- pMo <- pMu <- vector()
    tS <- tMo <- tMu <- vector()
    titer <- vector()
    for(a in 1:length(agesM)){ # estimate mothers prob infection
      min <- years[y]-floor(agesM[a])
      max <- years[y]
      ifoi <- sum(lam[which(lamY==min):which(lamY==max)]) 
      pS[a] <- exp(-4*ifoi)
      pMo[a] <- 4*exp(-3*ifoi)*(1-exp(-ifoi))
      pMu[a] <- 1 - pS[a] - pMo[a]
      tS[a] <- rtruncnorm(1,a=0,b=Inf,5,5)
      tMo[a] <- rtruncnorm(1,a=0,b=Inf,40,10)
      tMu[a] <- rtruncnorm(1,a=0,b=Inf,100,20)
      titer[a] <- mean(tS[a]*pS[a] + tMo[a]*pMo[a] + tMu[a]*pMu[a])
    }
    probSus[i,y] <- mean(pS)
    probMonos[i,y] <- mean(pMo)
    probMultis[i,y] <- mean(pMu)
    titers[[y]] <- titer
  }
  #TiterList[i] <- titers
}

# mean titer at birth over time
mea <- vector()
for(i in 1:40) mea[i] <- mean(titers[[i]])
plot(mea)

# model the decay of maternal abs
tt <- seq(0,12,0.1)
log(2)/0.55 * 30.4
HL <- 0.55
logmu <- 3.06
logsd <- 1.18
RiskL <- list()
for(y in 1:nYears){
  
  n <- length(titers[[y]])
  abs <- matrix(NA, nrow=length(tt), ncol=n)
  risk <- abs
  for(i in 1:n){
    for(t in 1:length(tt)){
      abs[t,i] <- titers[[y]][i]*exp(-HL*tt[t])
      risk[t,i] <- dlnorm(abs[t,i], logmu, logsd)
      
    }
  }
  for(t in 1:nrow(risk)) risk[t,1:n] <- risk[t,1:n]*(1-exp(-FOI[y]*tt[t]/12))
  RiskL[[y]] <- risk
}
matplot(RiskL[[1]])

# calculate mean age
ages <- vector()
AgeY <- data.frame(years=seq(1980,2019), mean=NA, ciL=NA, ciU=NA)
SumRisk <- data.frame(years=seq(1980,2019),sum=NA,sumav=NA)
for(y in 1:40){
  for(i in 1:ncol(RiskL[[y]])){
    ages[i] <- tt[which.max(RiskL[[y]][,i])]
  }
  AgeY$mean[y] <- mean(ages)
  AgeY[y,3:4] <- t.test(ages)$conf.int[1:2]
  SumRisk$sum[y] <- sum(RiskL[[y]])
  SumRisk$sumav[y] <- sum(RiskL[[y]])/ncol(RiskL[[y]])
}
AgeY
SumRisk
p1 <- ggplot(AgeY, aes(years, mean))+ geom_point(col='seagreen4', size=2.5)+ ylim(0,8)+
  geom_linerange(aes(years, ymin=ciL,ymax=ciU),col='seagreen4')+ theme_minimal()+
  theme(text=element_text(size=22))+ ylab('Mean age (months)')+ xlab('Year')
p1
p1zoom <- ggplot(AgeY, aes(years, mean))+ geom_point(col='seagreen4', size=2.5)+ #ylim(0,10)+
  geom_linerange(aes(years, ymin=ciL,ymax=ciU),col='seagreen4')+ theme_bw()+
  theme(text=element_text(size=22))+ ylab('Mean age (months)')+ xlab('Year')
p1zoom

PLOT <- ggdraw() + draw_plot(p1) +
  draw_plot(p1zoom, x=.35, y=.15, width=.6, height=.4)
PLOT

ggplot(SumRisk, aes(years,sum))+ geom_point()
ggplot(SumRisk, aes(years,sumav))+ geom_point()



#----- Repeat with constant maternal titers

RiskL2 <- list()
for(y in 1:nYears){
  
  n <- length(titers[[y]])
  abs <- matrix(NA, nrow=length(tt), ncol=n)
  risk <- abs
  for(i in 1:n){
    for(t in 1:length(tt)){
      abs[t,i] <- titers[[1]][i]*exp(-HL*tt[t])
      risk[t,i] <- dlnorm(abs[t,i], logmu, logsd)
      
    }
  }
  for(t in 1:nrow(risk)) risk[t,1:n] <- risk[t,1:n]*(1-exp(-FOI[y]*tt[t]/12))
  RiskL2[[y]] <- risk
}
matplot(RiskL[[1]])

# calculate mean age
ages <- vector()
AgeY <- data.frame(years=seq(1980,2019), mean=NA, ciL=NA, ciU=NA)
for(y in 1:40){
  for(i in 1:ncol(RiskL[[y]])){
    ages[i] <- tt[which.max(RiskL2[[y]][,i])]
  }
  AgeY$mean[y] <- mean(ages)
  AgeY[y,3:4] <- t.test(ages)$conf.int[1:2]
  
}
AgeY
p2 <- ggplot(AgeY, aes(years, mean))+ geom_point()+ ylim(0,8)+
  geom_linerange(aes(years, ymin=ciL,ymax=ciU))+ theme_minimal()+
  theme(text=element_text(size=22))+ ylab('Mean age (months)')+ xlab('Year')
p2
p2zoom <- ggplot(AgeY, aes(years, mean))+ geom_point()+ #ylim(0,10)+
  geom_linerange(aes(years, ymin=ciL,ymax=ciU))+ theme_minimal()+
  theme(text=element_text(size=22))+ ylab('Mean age (months)')+ xlab('Year')
p2zoom

# output
setwd('C:/Users/Megan/OneDrive - University of Cambridge/MaternalDENV/FinalPlots')
png(filename='SimMeanAge_constantTiters.png', width=25, height=12, res=400, units='cm')
plot(p2)
dev.off()
png(filename='SimMeanAge_constantTiters_zoom.png', width=25, height=12, res=400, units='cm')
plot(p2zoom)
dev.off()

