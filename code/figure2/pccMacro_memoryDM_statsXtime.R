########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
###               Apply to timeseries data              ####
###                        Mem - Choice                 ####
########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)
library(zoo)

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure2'
baseDir='<path2folder>/data/memoryDM'


#Load the choice data
pccDat=read.csv(sprintf('%s/pccMacro_memoryBasedDM_onsetLocked.csv',baseDir),header=T)
head(pccDat)

#### Main effect going by anatomy
regWofLoc=summaryBy(psc~time*elecID*locName,data=pccDat,FUN=c(mean))
names(regWofLoc)=c('time','elecID','locName','psc')
str(regWofLoc)

choicePeriodLoc=regWofLoc[(regWofLoc$time>1700) & (regWofLoc$time<2200),]
head(choicePeriodLoc)

cSum=summaryBy(psc~elecID*locName,FUN=c(mean),data=choicePeriodLoc)
cSum$locName=as.factor(cSum$locName)
cSum$Subject=substr(cSum$elecID,1,6)
cSum$Subject=as.factor(cSum$Subject)
str(cSum)

lmChoice0=lmer(psc.mean~locName+(1|Subject),data=cSum)
summary(lmChoice0)
confint(lmChoice0,method="Wald")

cd0=cSum[cSum$locName=='dPCC',]
cd0lm=lmer(psc.mean~1+(1|Subject),data=cd0)
summary(cd0lm) 
confint(cd0lm,method='Wald')

cv0=cSum[cSum$locName=='vPCC',]
cv0lm=lmer(psc.mean~1+(1|Subject),data=cv0)
summary(cv0lm) 
confint(cv0lm,method='Wald')

#######
## Main effect by cluster
### The base response for this time period (across safe and risky trials)
regWofMain1=summaryBy(psc~time*elecID*blockValue*dendClust,data=pccDat,FUN=c(mean))
str(regWofMain1)
regWofMain2=summaryBy(psc.mean~time*elecID*dendClust,data=regWofMain1,FUN=c(mean))
names(regWofMain2)=c('time','elecID','clustVal','psc')
str(regWofMain2)

choicePeriodMain=regWofMain2[(regWofMain2$time>1700) & (regWofMain2$time<2200),]
head(choicePeriodMain)

regWofMain3=summaryBy(psc~elecID*clustVal,FUN=c(mean),data=choicePeriodMain)
regWofMain3$clustVal=as.factor(regWofMain3$clustVal)
regWofMain3b=regWofMain3[regWofMain3$clustVal=='1' | regWofMain3$clustVal=='2',]
regWofMain3b$Subject=substr(regWofMain3b$elecID,1,6)
regWofMain3b$Subject=as.factor(regWofMain3b$Subject)
str(regWofMain3b)

lmChoice0=lmer(psc.mean~clustVal+(1|Subject),data=regWofMain3b)
summary(lmChoice0)
confint(lmChoice0,method="Wald")

cd0=regWofMain3b[regWofMain3b$clustVal=='1',]
cd0lm=lmer(psc.mean~1+(1|Subject),data=cd0)
summary(cd0lm) #
confint(cd0lm,method='Wald')

cv0=regWofMain3b[regWofMain3b$clustVal=='2',]
cv0lm=lmer(psc.mean~1+(1|Subject),data=cv0)
summary(cv0lm)
confint(cv0lm,method='Wald')


########
### Risk-modulation by cluster
regWof1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName*dendClust,data=pccDat,FUN=c(mean))
str(regWof1)
regWof2=summaryBy(psc.mean~time*elecID*wasRisky*locName*dendClust,data=regWof1,FUN=c(mean))
names(regWof2)=c('time','elecID','wasRisky','locName','dendClust','psc')


#Now calculate safe vs. risky during the choice period
choicePeriod=regWof2[(regWof2$time>1700) & (regWof2$time<2200),]
head(choicePeriod)
choiceSum1=summaryBy(psc~elecID*locName*wasRisky*dendClust,data=choicePeriod,FUN=c(mean))
choiceSum1$timeperiod='onset'
choiceSum1$timeperiod=as.factor(choiceSum1$timeperiod)
choiceSum1$wasRisky=as.factor(choiceSum1$wasRisky)
choiceSum1$locName=as.factor(choiceSum1$locName)
names(choiceSum1)=c('elecID','locName','wasRisky','clustVal','psc','timeperiod')
choiceSum1$clustVal=as.factor(choiceSum1$clustVal)
str(choiceSum1)

choiceSum2=dcast(data=choiceSum1,elecID+clustVal+locName+timeperiod~wasRisky,value.var="psc")
names(choiceSum2)=c('elecID','clustVal','locName','timeperiod','safe','risky')
choiceSum2$diff=choiceSum2$risky-choiceSum2$safe
choiceSum2$locClust=paste(choiceSum2$locName,choiceSum2$clustVal,sep='_')
choiceSum2$locClust=as.factor(choiceSum2$locClust)

#Safe versus risky for the entire cluster (1/2)
choiceSum2$Subject=substr(choiceSum2$elecID,1,6)
choiceSum2$Subject=as.factor(choiceSum2$Subject)
choiceSum2all=choiceSum2[choiceSum2$clustVal=='1' |choiceSum2$clustVal=='2',]

lmChoice0=lmer(diff~clustVal+(1|Subject),data=choiceSum2all)
summary(lmChoice0)
confint(lmChoice0,method="Wald")

cd0=choiceSum2all[choiceSum2all$clustVal=='1',]
cd0lm=lmer(diff~1+(1|Subject),data=cd0)
summary(cd0lm)
confint(cd0lm,method='Wald')

cv0=choiceSum2all[choiceSum2all$clustVal=='2',]
cv0lm=lmer(diff~1+(1|Subject),data=cv0)
summary(cv0lm)
confint(cv0lm,method='Wald')


### Quick method for visualizing results
bothRiskH=summaryBy(psc~time*dendClust*wasRisky,data=regWof2,FUN=c(nanmean,se,length))
names(bothRiskH)=c('time','clustVal','wasRisky','psc','se','len')
bothRiskH=bothRiskH[bothRiskH$clustVal==1 | bothRiskH$clustVal==2,]

smoothRiskH=data.frame()
locListH=unique(bothRiskH$clustVal)
for(locI in locListH){
  locDat=bothRiskH[bothRiskH$clustVal==locI,]
  
  for(riskI in c("0","1")){
    riskDat=locDat[locDat$wasRisky==riskI,]
    
    riskDat$smoothDat=rollapply(riskDat$psc, width=200, FUN=mean, align="center", partial=TRUE)
    riskDat$smoothSE=rollapply(riskDat$se, width=200, FUN=mean, align="center", partial=TRUE)
    smoothRiskH=rbind(smoothRiskH,riskDat)
  }
}

smoothRiskH$wasRisky=as.factor(smoothRiskH$wasRisky)
smoothRiskH$clustVal=as.factor(smoothRiskH$clustVal)
str(smoothRiskH)

smoothRiskHb=smoothRiskH[smoothRiskH$time>=1000,]
smoothRiskHb=smoothRiskHb[smoothRiskHb$time<=3000,]

bothRiskyPlotH1=ggplot() + theme_bw() + geom_line(data=smoothRiskHb,aes(x=time,y=smoothDat,color=clustVal,linetype=wasRisky)) + geom_ribbon(data=smoothRiskHb,aes(x=time,ymin=smoothDat-smoothSE,ymax=smoothDat+smoothSE,fill=clustVal,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='red') + scale_y_continuous(limits=c(-12,22),breaks=c(-10,-5,0,5,10,15,20)) + theme(panel.grid.minor = element_blank())
bothRiskyPlotH1
riskyH1name=sprintf('%s/macroPCC_memory_choice_smoothPSC_xRisky_xClust.pdf',outDir)
ggsave(riskyH1name,height=7,width=10,bothRiskyPlotH1)


# Keep the summary stat for later plotting across all timeperiods
allSum=choiceSum2all
keep='allSum'
rm(list = setdiff(ls(), keep))


########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
###               Apply to timeseries data              ####
###                       Mem - Delay                   ####
########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure2'
baseDir='<path2folder>/data/memoryDM'

#Load the choice data
pccDat=read.csv(sprintf('%s/pccMacro_memoryBasedDM_responseLocked.csv',baseDir),header=T)
head(pccDat)

#### Main effect going by anatomy
regWofLoc=summaryBy(psc~time*elecID*locName,data=pccDat,FUN=c(mean))
names(regWofLoc)=c('time','elecID','locName','psc')

respPeriodLoc=regWofLoc[(regWofLoc$time>1500) & (regWofLoc$time<2000),]
head(respPeriodLoc)

rSum=summaryBy(psc~elecID*locName,FUN=c(mean),data=respPeriodLoc)
rSum$locName=as.factor(rSum$locName)
rSum$Subject=substr(rSum$elecID,1,6)
rSum$Subject=as.factor(rSum$Subject)
str(rSum)

lmResp0=lmer(psc.mean~locName+(1|Subject),data=rSum)
summary(lmResp0)
confint(lmResp0,method="Wald")

rd0=rSum[rSum$locName=='dPCC',]
rd0lm=lmer(psc.mean~1+(1|Subject),data=rd0)
summary(rd0lm) 
confint(rd0lm,method='Wald')

rv0=rSum[rSum$locName=='vPCC',]
rv0lm=lmer(psc.mean~1+(1|Subject),data=rv0)
summary(rv0lm) 
confint(rv0lm,method='Wald')

#######
## Main effect by cluster
### The base response for this time period (across safe and risky trials)
regWofMain1=summaryBy(psc~time*elecID*blockValue*dendClust,data=pccDat,FUN=c(mean))
str(regWofMain1)
regWofMain2=summaryBy(psc.mean~time*elecID*dendClust,data=regWofMain1,FUN=c(mean))
names(regWofMain2)=c('time','elecID','clustVal','psc')
str(regWofMain2)

respPeriodMain=regWofMain2[(regWofMain2$time>1500) & (regWofMain2$time<2000),]
head(respPeriodMain)

regWofMain3=summaryBy(psc~elecID*clustVal,FUN=c(mean),data=respPeriodMain)
regWofMain3$clustVal=as.factor(regWofMain3$clustVal)
regWofMain3b=regWofMain3[regWofMain3$clustVal=='1' | regWofMain3$clustVal=='2',]
regWofMain3b$Subject=substr(regWofMain3b$elecID,1,6)
regWofMain3b$Subject=as.factor(regWofMain3b$Subject)
str(regWofMain3b)

lmResp0b=lmer(psc.mean~clustVal+(1|Subject),data=regWofMain3b)
summary(lmResp0b)
confint(lmResp0b,method="Wald")

rd0m=regWofMain3b[regWofMain3b$clustVal=='1',]
rdLMm=lmer(psc.mean~1+(1|Subject),data=rd0m)
summary(rdLMm)
confint(rdLMm,method='Wald')

rv0m=regWofMain3b[regWofMain3b$clustVal=='2',]
rvLMm=lmer(psc.mean~1+(1|Subject),data=rv0m)
summary(rvLMm)
confint(rvLMm,method='Wald')


########
### Risk-modulation by cluster
regWof1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName*dendClust,data=pccDat,FUN=c(mean))
str(regWof1)
regWof2=summaryBy(psc.mean~time*elecID*wasRisky*locName*dendClust,data=regWof1,FUN=c(mean))
names(regWof2)=c('time','elecID','wasRisky','locName','dendClust','psc')


#Now calculate safe vs. risky during the choice period
respPeriod=regWof2[(regWof2$time>1500) & (regWof2$time<2000),]
head(respPeriod)
respSum1=summaryBy(psc~elecID*locName*wasRisky*dendClust,data=respPeriod,FUN=c(mean))
respSum1$timeperiod='response'
respSum1$timeperiod=as.factor(respSum1$timeperiod)
respSum1$wasRisky=as.factor(respSum1$wasRisky)
respSum1$locName=as.factor(respSum1$locName)
names(respSum1)=c('elecID','locName','wasRisky','clustVal','psc','timeperiod')
respSum1$clustVal=as.factor(respSum1$clustVal)
str(respSum1)

respSum2=dcast(data=respSum1,elecID+clustVal+locName+timeperiod~wasRisky,value.var="psc")
names(respSum2)=c('elecID','clustVal','locName','timeperiod','safe','risky')
respSum2$diff=respSum2$risky-respSum2$safe
respSum2$locClust=paste(respSum2$locName,respSum2$clustVal,sep='_')
respSum2$locClust=as.factor(respSum2$locClust)

#Safe versus risky for the entire cluster (1/2)
respSum2$Subject=substr(respSum2$elecID,1,6)
respSum2$Subject=as.factor(respSum2$Subject)
respSum2all=respSum2[respSum2$clustVal=='1' |respSum2$clustVal=='2',]

lmResp3=lmer(diff~clustVal+(1|Subject),data=respSum2all)
summary(lmResp3)
confint(lmResp3,method="Wald")

rd0r=respSum2all[respSum2all$clustVal=='1',]
rdLMr=lmer(diff~1+(1|Subject),data=rd0r)
summary(rdLMr)
confint(rdLMr,method='Wald')

rv0r=respSum2all[respSum2all$clustVal=='2',]
rvLMr=lmer(diff~1+(1|Subject),data=rv0r)
summary(rvLMr)
confint(rvLMr,method='Wald')


### Quick method for visualizing results
bothRiskH=summaryBy(psc~time*dendClust*wasRisky,data=regWof2,FUN=c(nanmean,se,length))
names(bothRiskH)=c('time','clustVal','wasRisky','psc','se','len')
bothRiskH=bothRiskH[bothRiskH$clustVal==1 | bothRiskH$clustVal==2,]

smoothRiskH=data.frame()
locListH=unique(bothRiskH$clustVal)
for(locI in locListH){
  locDat=bothRiskH[bothRiskH$clustVal==locI,]
  
  for(riskI in c("0","1")){
    riskDat=locDat[locDat$wasRisky==riskI,]
    
    riskDat$smoothDat=rollapply(riskDat$psc, width=200, FUN=mean, align="center", partial=TRUE)
    riskDat$smoothSE=rollapply(riskDat$se, width=200, FUN=mean, align="center", partial=TRUE)
    smoothRiskH=rbind(smoothRiskH,riskDat)
  }
}

smoothRiskH$wasRisky=as.factor(smoothRiskH$wasRisky)
smoothRiskH$clustVal=as.factor(smoothRiskH$clustVal)
str(smoothRiskH)

smoothRiskHb=smoothRiskH[smoothRiskH$time>=500,]

bothRiskyPlotH1=ggplot() + theme_bw() + geom_line(data=smoothRiskHb,aes(x=time,y=smoothDat,color=clustVal,linetype=wasRisky)) + geom_ribbon(data=smoothRiskHb,aes(x=time,ymin=smoothDat-smoothSE,ymax=smoothDat+smoothSE,fill=clustVal,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='red') + scale_y_continuous(limits=c(-12,22),breaks=c(-10,-5,0,5,10,15,20)) + theme(panel.grid.minor = element_blank())
bothRiskyPlotH1
riskyH1name=sprintf('%s/macroPCC_memory_response_smoothPSC_xRisky_xClust.pdf',outDir)
ggsave(riskyH1name,height=7,width=10,bothRiskyPlotH1)

# Keep the summary stat for later plotting across all timeperiods
allSum=rbind(allSum,respSum2all)
keep='allSum'
rm(list = setdiff(ls(), keep))

########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
###               Apply to timeseries data              ####
###                       Mem - Feedback                   ####
########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure2'
baseDir='<path2folder>/data/memoryDM'

#Load the choice data
pccDat=read.csv(sprintf('%s/pccMacro_memoryBasedDM_feedbackLocked.csv',baseDir),header=T)
head(pccDat)

#### Main effect going by anatomy
regMemLoc=summaryBy(psc~time*elecID*locName,data=pccDat,FUN=c(mean))
names(regMemLoc)=c('time','elecID','locName','psc')
str(regMemLoc)

feedPeriodLoc=regMemLoc[(regMemLoc$time>1500) & (regMemLoc$time<2000),]
head(feedPeriodLoc)

fSum=summaryBy(psc~elecID*locName,FUN=c(mean),data=feedPeriodLoc)
fSum$locName=as.factor(fSum$locName)
fSum$Subject=substr(fSum$elecID,1,6)
fSum$Subject=as.factor(fSum$Subject)
str(fSum)

lmResp0=lmer(psc.mean~locName+(1|Subject),data=fSum)
summary(lmResp0)
confint(lmResp0,method="Wald")

fd0=fSum[fSum$locName=='dPCC',]
fd0lm=lmer(psc.mean~1+(1|Subject),data=fd0)
summary(fd0lm) 
confint(fd0lm,method='Wald')

fv0=fSum[fSum$locName=='vPCC',]
fv0lm=lmer(psc.mean~1+(1|Subject),data=fv0)
summary(fv0lm) 
confint(fv0lm,method='Wald')

#######
## Main effect by cluster
### The base response for this time period (across safe and risky trials)
regWofMain1=summaryBy(psc~time*elecID*blockValue*dendClust,data=pccDat,FUN=c(mean))
str(regWofMain1)
regWofMain2=summaryBy(psc.mean~time*elecID*dendClust,data=regWofMain1,FUN=c(mean))
names(regWofMain2)=c('time','elecID','clustVal','psc')
str(regWofMain2)

feedPeriodMain=regWofMain2[(regWofMain2$time>1500) & (regWofMain2$time<2000),]
head(feedPeriodMain)

regWofMain3=summaryBy(psc~elecID*clustVal,FUN=c(mean),data=feedPeriodMain)
regWofMain3$clustVal=as.factor(regWofMain3$clustVal)
regWofMain3b=regWofMain3[regWofMain3$clustVal=='1' | regWofMain3$clustVal=='2',]
regWofMain3b$Subject=substr(regWofMain3b$elecID,1,6)
regWofMain3b$Subject=as.factor(regWofMain3b$Subject)
str(regWofMain3b)

lmFeed0=lmer(psc.mean~clustVal+(1|Subject),data=regWofMain3b)
summary(lmFeed0)
confint(lmFeed0,method="Wald")

fd0=regWofMain3b[regWofMain3b$clustVal=='1',]
fd0lm=lmer(psc.mean~1+(1|Subject),data=fd0)
summary(fd0lm) #
confint(fd0lm,method='Wald')

fv0=regWofMain3b[regWofMain3b$clustVal=='2',]
fv0lm=lmer(psc.mean~1+(1|Subject),data=fv0)
summary(fv0lm)
confint(fv0lm,method='Wald')


########
### Risk-modulation by cluster
regMem1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName*dendClust,data=pccDat,FUN=c(mean))
str(regMem1)
regMem2=summaryBy(psc.mean~time*elecID*wasRisky*locName*dendClust,data=regMem1,FUN=c(mean))
names(regMem2)=c('time','elecID','wasRisky','locName','dendClust','psc')


#Now calculate safe vs. risky during the choice period
feedPeriod=regMem2[(regMem2$time>1500) & (regMem2$time<2000),]
head(feedPeriod)
feedSum1=summaryBy(psc~elecID*locName*wasRisky*dendClust,data=feedPeriod,FUN=c(mean))
feedSum1$timeperiod='feedback'
feedSum1$timeperiod=as.factor(feedSum1$timeperiod)
feedSum1$wasRisky=as.factor(feedSum1$wasRisky)
feedSum1$locName=as.factor(feedSum1$locName)
names(feedSum1)=c('elecID','locName','wasRisky','clustVal','psc','timeperiod')
feedSum1$clustVal=as.factor(feedSum1$clustVal)
str(feedSum1)

feedSum2=dcast(data=feedSum1,elecID+clustVal+locName+timeperiod~wasRisky,value.var="psc")
names(feedSum2)=c('elecID','clustVal','locName','timeperiod','safe','risky')
feedSum2$diff=feedSum2$risky-feedSum2$safe
feedSum2$locClust=paste(feedSum2$locName,feedSum2$clustVal,sep='_')
feedSum2$locClust=as.factor(feedSum2$locClust)

#Safe versus risky for the entire cluster (1/2)
feedSum2$Subject=substr(feedSum2$elecID,1,6)
feedSum2$Subject=as.factor(feedSum2$Subject)
feedSum2all=feedSum2[feedSum2$clustVal=='1' |feedSum2$clustVal=='2',]

lmFeed1=lmer(diff~clustVal+(1|Subject),data=feedSum2all)
summary(lmFeed1)
confint(lmFeed1,method="Wald")

fd0r=feedSum2all[feedSum2all$clustVal=='1',]
fd0rLM=lmer(diff~1+(1|Subject),data=fd0r)
summary(fd0rLM)
confint(fd0rLM,method='Wald')

fv0r=feedSum2all[feedSum2all$clustVal=='2',]
fv0rLM=lmer(diff~1+(1|Subject),data=fv0r)
summary(fv0rLM)
confint(fv0rLM,method='Wald')


### Quick method for visualizing results
bothRiskH=summaryBy(psc~time*dendClust*wasRisky,data=regWof2,FUN=c(nanmean,se,length))
names(bothRiskH)=c('time','clustVal','wasRisky','psc','se','len')
bothRiskH=bothRiskH[bothRiskH$clustVal==1 | bothRiskH$clustVal==2,]

smoothRiskH=data.frame()
locListH=unique(bothRiskH$clustVal)
for(locI in locListH){
  locDat=bothRiskH[bothRiskH$clustVal==locI,]
  
  for(riskI in c("0","1")){
    riskDat=locDat[locDat$wasRisky==riskI,]
    
    riskDat$smoothDat=rollapply(riskDat$psc, width=100, FUN=mean, align="center", partial=TRUE)
    riskDat$smoothSE=rollapply(riskDat$se, width=100, FUN=mean, align="center", partial=TRUE)
    smoothRiskH=rbind(smoothRiskH,riskDat)
  }
}

smoothRiskH$wasRisky=as.factor(smoothRiskH$wasRisky)
smoothRiskH$clustVal=as.factor(smoothRiskH$clustVal)
str(smoothRiskH)

smoothRiskHb=smoothRiskH[smoothRiskH$time>=500,]

bothRiskyPlotH1b=ggplot() + theme_bw() + geom_line(data=smoothRiskHb,aes(x=time,y=smoothDat,color=clustVal,linetype=wasRisky)) + geom_ribbon(data=smoothRiskHb,aes(x=time,ymin=smoothDat-smoothSE,ymax=smoothDat+smoothSE,fill=clustVal,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='solid',color='black') + scale_y_continuous(limits=c(-12,22),breaks=c(-10,-5,0,5,10,15,20)) + theme(panel.grid.minor = element_blank())
bothRiskyPlotH1b
riskyH1name=sprintf('%s/macroPCC_mem_feedback_smoothPSC_xRisky_xClust.pdf',outDir)
ggsave(riskyH1name,height=7,width=10,bothRiskyPlotH1b)

######
# Plot risk-modulation across time periods WoF
allSum=rbind(allSum,feedSum2all)

#Plot risky-safe differences all together
head(allSum)

allSum$plotX=NaN
allSum$plotX=ifelse((allSum$timeperiod=='onset') & (allSum$clustVal=='1'),1,allSum$plotX)
allSum$plotX=ifelse((allSum$timeperiod=='onset') & (allSum$clustVal=='2'),2,allSum$plotX)
allSum$plotX=ifelse((allSum$timeperiod=='response') & (allSum$clustVal=='1'),4,allSum$plotX)
allSum$plotX=ifelse((allSum$timeperiod=='response') & (allSum$clustVal=='2'),5,allSum$plotX)
allSum$plotX=ifelse((allSum$timeperiod=='feedback') & (allSum$clustVal=='1'),7,allSum$plotX)
allSum$plotX=ifelse((allSum$timeperiod=='feedback') & (allSum$clustVal=='2'),8,allSum$plotX)

allSum2=summaryBy(diff~clustVal*plotX,data=allSum,FUN=c(mean,se))

allDiffPlot=ggplot() + theme_bw() + geom_jitter(data=allSum,aes(x=plotX,y=diff,color=clustVal),width=.1,height=0,size=3,alpha=.25) + geom_bar(data=allSum2,stat='identity',aes(fill=clustVal,x=plotX,y=diff.mean),alpha=.3) + geom_errorbar(data=allSum2,stat='identity',aes(x=plotX,ymin=diff.mean-diff.se,ymax=diff.mean+diff.se),width=.1,color='black') + scale_y_continuous(limits=c(-10,15),breaks=c(-10,-5,0,5,10,15))
allDiffPlot
allDiffName=sprintf('%s/macroPCC_memDM_riskModSummary_allTimes_clust12.pdf',outDir)
ggsave(allDiffName,height=7,width=10,allDiffPlot)
