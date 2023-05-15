#----- Function to simulate temporal DENV infant risk trends -----# 



#--- Simulate from ADE model (lognormal)
simTrendsADE <- function(Yrs, YrsH, YrInd, FOI, FOIH, FTA, nAgeG, ageL, ageU, tt, HL, logmu, logsd){
  
  # lists to store simulations
  titers <- list()
  Risk <- list()
  RiskDI <- list()
  statusP <- data.frame(years=Yrs, Su=NA, Mo=NA, Mu=NA)
  
  # for each year
  for(y in 1:length(Yrs)){
    
    # age-specific number of births in year y
    ratesA <- FTA[FTA$year==Yrs[y],]
    
    # age dist of birthing women
    ageB <- list()
    for(a in 1:nAgeG) ageB[[a]] <- rep(ages[a], floor(ratesA$rate[ratesA$age==ages[a]]))
    ageB <- unlist(ageB)
    
    # vectors to store probabilities & titers
    pS <- pMo <- pMu <- vector()
    tS <- tMo <- tMu <- vector()
    Mtiter <- vector()
    
    # for each mother, simulate a DENV titer
    for(i in 1:length(ageB)){ 
      
      # probability of infection, given age
      minT <- Yrs[y]-floor(ageB[i])
      maxT <- Yrs[y]
      ifoi <- sum(FOIH[which(YrsH==minT):which(YrsH==maxT)]) 
      
      # probability of being susecptible/monotypic/multitypic
      pS[i] <- exp(-4*ifoi)
      pMo[i] <- 4*exp(-3*ifoi)*(1-exp(-ifoi))
      pMu[i] <- 1 - pS[i] - pMo[i]
      
      # draw a weighted mean DENV titer 
      tS[i] <- rtruncnorm(1,a=0,b=Inf,20,20)
      tMo[i] <- rtruncnorm(1,a=0,b=Inf,150,100)
      tMu[i] <- rtruncnorm(1,a=0,b=Inf,500,200)
      Mtiter[i] <- mean(tS[i]*pS[i] + tMo[i]*pMo[i] + tMu[i]*pMu[i])
    }
    
    # store annual probabilities & titer values
    statusP$Su[y] <- mean(pS)
    statusP$Mo[y] <- mean(pMo)
    statusP$Mu[y] <- mean(pMu)
    titers[[y]] <- Mtiter
    
    # model the decay of infant maternally-derived abs
    n <- length(titers[[y]])
    Iabs <- Irisk <- IriskFOI <- matrix(NA, nrow=length(tt), ncol=n)
    for(i in 1:n){
      for(t in 1:length(tt)){
        
        # infant titer at each time point
        Iabs[t,i] <- titers[[y]][i]*exp(-HL*tt[t])
        
        # infant titer-risk density at each time point
        Irisk[t,i] <- dlnorm(Iabs[t,i], logmu, logsd)/dlnorm(exp(logmu-logsd^2), logmu, logsd)
      }
    }
    
    # multiply titer-risk density by infant risk of infection
    for(t in 2:nrow(Irisk)) IriskFOI[t,1:n] <- Irisk[t,1:n]*(exp(-4*FOI[y]*tt[t-1]/12))*(4*FOI[y]/12)
    Risk[[y]] <- Irisk
    RiskDI[[y]] <- IriskFOI
  }
  
  return(list(Risk=Risk, RiskDI=RiskDI, Mtiters=titers, statusP=statusP))
}



#--- Simulate from threshold (non-ADE) model 
simTrendsT <- function(Yrs, YrsH, YrInd, FOI, FOIH, FTA, nAgeG, ageL, ageU, tt, HL, mu, sd, alpha){
  
  # lists to store simulations
  titers <- list()
  Risk <- list()
  RiskDI <- list()
  statusP <- data.frame(years=Yrs, Su=NA, Mo=NA, Mu=NA)
  
  # for each year
  for(y in 1:length(Yrs)){
    
    # age-specific number of births in year y
    ratesA <- FTA[FTA$year==Yrs[y],]
    
    # age dist of birthing women
    ageB <- list()
    for(a in 1:nAgeG) ageB[[a]] <- rep(ages[a], floor(ratesA$rate[ratesA$age==ages[a]]))
    ageB <- unlist(ageB)
    
    # vectors to store probabilities & titers
    pS <- pMo <- pMu <- vector()
    tS <- tMo <- tMu <- vector()
    Mtiter <- vector()
    
    # for each mother, simulate a DENV titer
    for(i in 1:length(ageB)){ 
      
      # probability of infection, given age
      minT <- Yrs[y]-floor(ageB[i])
      maxT <- Yrs[y]
      ifoi <- sum(FOIH[which(YrsH==minT):which(YrsH==maxT)]) 
      
      # probability of being susecptible/monotypic/multitypic
      pS[i] <- exp(-4*ifoi)
      pMo[i] <- 4*exp(-3*ifoi)*(1-exp(-ifoi))
      pMu[i] <- 1 - pS[i] - pMo[i]
      
      # draw a weighted mean DENV titer 
      tS[i] <- rtruncnorm(1,a=0,b=Inf,20,20)
      tMo[i] <- rtruncnorm(1,a=0,b=Inf,150,100)
      tMu[i] <- rtruncnorm(1,a=0,b=Inf,500,200)
      Mtiter[i] <- mean(tS[i]*pS[i] + tMo[i]*pMo[i] + tMu[i]*pMu[i])
    }

        # store annual probabilities & titer values
    statusP$Su[y] <- mean(pS)
    statusP$Mo[y] <- mean(pMo)
    statusP$Mu[y] <- mean(pMu)
    titers[[y]] <- Mtiter
    
    # model the decay of infant maternally-derived abs
    n <- length(titers[[y]])
    Iabs <- Irisk <- IriskFOI <- matrix(NA, nrow=length(tt), ncol=n)
    for(i in 1:n){
      for(t in 1:length(tt)){
        
        # infant titer at each time point
        Iabs[t,i] <- titers[[y]][i]*exp(-HL*tt[t])
        
        # infant titer-risk density at each time point
        Irisk[t,i] <- 1-pnorm(Iabs[t,i], mu, sd)
      }
    }
    
    # multiply titer-risk density by infant risk of infection
    ageR <- exp(-alpha*tt)
    for(t in 2:nrow(Irisk)) IriskFOI[t,1:n] <- Irisk[t,1:n]*ageR[t]*(exp(-4*FOI[y]*tt[t-1]/12))*(4*FOI[y]/12)
    Risk[[y]] <- Irisk
    RiskDI[[y]] <- IriskFOI
  }
  
  return(list(Risk=Risk, RiskDI=RiskDI, Mtiters=titers, statusP=statusP))
}
