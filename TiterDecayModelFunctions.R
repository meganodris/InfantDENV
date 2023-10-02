#----- Functions for Maternal Antibody Decay Model -----#
library(scales)
library(ggplot2)


#--- Function to plot observed vs predicted titers
plot_ObsVsPred <- function(chains, data){
  
  # extract predicted titers
  pred <- pred2 <- pred3 <- data.frame(cohort=NA, sero=NA, obs=NA, pred=NA, ciL=NA, ciU=NA)
  for(i in 1:data$N_hi1) for(s in 1:4) for(t in 1:5) if(data$HIind1[s,i,t]==1){
    pred[nrow(pred)+1,2:6] <- c(s, data$obsHI1[s,i,t], quantile(chains$pHI1[,s,i,t], c(0.5,0.025,0.975)))
  }
  for(i in 1:dta$N_hi2) for(s in 1:4) for(t in 1:5) if(data$HIind2[s,i,t]==1){
    pred2[nrow(pred2)+1,2:6] <- c(s, data$obsHI2[s,i,t], quantile(chains$pHI2[,s,i,t], c(0.5,0.025,0.975)))
  }
  for(i in 1:dta$N_prnt) for(s in 1:4) for(t in 1:5) if(data$NTind[s,i,t]==1){
    pred3[nrow(pred3)+1,2:6] <- c(s, data$obsNT[s,i,t], quantile(chains$pNT[,s,i,t], c(0.5,0.025,0.975)))
  }
  
  # format data
  pred$cohort <- 'KPP'
  pred2$cohort <- 'BKHI'
  pred3$cohort <- 'BKNT'
  pred <- rbind(pred, pred2, pred3)
  pred <- pred[!is.na(pred$sero),]
  
  # transform titers to linear scale
  pred$obs <- exp((pred$obs-1)*log(2)+log(10))
  pred$pred <- exp((pred$pred-1)*log(2)+log(10))
  pred$ciL <- exp((pred$ciL-1)*log(2)+log(10))
  pred$ciU <- exp((pred$ciU-1)*log(2)+log(10))
  
  # plot PRNTs
  prnt <- ggplot(pred[pred$cohort=='BKNT',], aes(obs, pred))+ 
    geom_rect(aes(xmin=3,xmax=10,ymin=1e-3,ymax=10),fill='grey90')+ geom_point(alpha=0.3, col='navy')+ 
    theme_minimal()+ theme(text=element_text(size=22))+ ylab('Predicted titer')+ xlab('Observed titer')+
    scale_y_continuous(trans='log', breaks=c(1,10,100,1000,10000),labels=comma)+
    scale_x_continuous(trans='log', breaks=c(1,10,100,1000,10000),labels=comma)+
    geom_line(aes(obs, obs), linetype='dashed')+ labs(subtitle='Bangkok PRNT')
  
  # plot HIs
  predHI <- pred[pred$cohort %in% c('KPP','BKHI'),]
  predHI$obs <- factor(predHI$obs, labels=c('<10','10','20','40','80','160','320','640','1280','2560','5120'))
  kpp <- ggplot(predHI[predHI$cohort=='KPP',], aes(obs,pred))+
    geom_point(alpha=0.3, col='navy')+ 
    theme_minimal()+ theme(text=element_text(size=22))+ ylab('Predicted titer')+ xlab('Observed titer')+
    scale_y_continuous(trans='log', breaks=c(10,40,160,640,2560,5120),labels=comma)+
    labs(subtitle='KPP HI')
  bkHI <- ggplot(predHI[predHI$cohort=='BKHI',], aes(obs,pred))+
    geom_point(alpha=0.3, col='navy')+ 
    theme_minimal()+ theme(text=element_text(size=22))+ ylab('Predicted titer')+ xlab('Observed titer')+
    scale_y_continuous(trans='log', breaks=c(10,40,160,640,2560,5120),labels=comma)+
    labs(subtitle='Bangkok HI')
  
  return(list(HI1=kpp, HI2=bkHI, NT=prnt))
}



