rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

#Load the data
outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/valueDM/hippocampus'

hippVal <- fread(sprintf('%s/hippMacro_valueBasedDM_choiceLocked.csv',baseDir), data.table = FALSE)
head(hippVal)

regWof1=summaryBy(psc~time*elecID*blockValue*locName,data=hippVal,FUN=c(mean))
str(regWof1)

regWof2=summaryBy(psc.mean~time*elecID*locName,data=regWof1,FUN=c(mean))
names(regWof2)=c('time','elecID','locName','psc')
regWof2$elecID=as.factor(regWof2$elecID)
str(regWof2)


#Now also keep the average difference for barplots
choicePeriod=regWof2[(regWof2$time>1700) & (regWof2$time<2201),]
choiceSum2=summaryBy(psc~elecID*locName,data=choicePeriod,FUN=c(mean))
names(choiceSum2)=c('elecID','locName','psc')
head(choiceSum2)
choiceSum2$Subject=substr(choiceSum2$elecID,1,6)
choiceSum2$Subject=as.factor(choiceSum2$Subject)
choiceSum2$locName=as.factor(choiceSum2$locName)
str(choiceSum2)

lmC1=lmer(psc~locName+(1|Subject),data=choiceSum2)
summary(lmC1) 
confint(lmC1,method='Wald')

pDat=choiceSum2[choiceSum2$locName=='postHipp',]
lmCp=lmer(psc~1+(1|Subject),data=pDat)
summary(lmCp) 
confint(lmCp,method='Wald')

aDat=choiceSum2[choiceSum2$locName=='antHipp',]
lmCa=lmer(psc~1+(1|Subject),data=aDat)
summary(lmCa)
confint(lmCa,method='Wald')


############
############
### Risk Modulation Hippocampus Choice
############
############
to_keep <- c("outDir", "baseDir","nanmean","se","hippVal")
rm(list = setdiff(ls(), to_keep))
regWof1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName,data=hippVal,FUN=c(mean))
str(regWof1)

#Smooth the data for plotting
smoothRisk=data.frame()
smoothRisk <- do.call(rbind, lapply(split(regWof1,with(regWof1, paste(elecID, blockValue, wasRisky, sep = "_"))),function(chunk) {
  chunk <- chunk[order(chunk$time), ]
  chunk$smoothDat <- ksmooth(chunk$time, chunk$psc.mean, kernel = "normal", bandwidth = 100)$y
  chunk
}))
#Make sure data is in correct order of time
smoothRisk <- smoothRisk[order(smoothRisk$elecID, 
                               smoothRisk$blockValue, 
                               smoothRisk$wasRisky), ]
rownames(smoothRisk) <- NULL

#Collapse for plotting
regWof2=summaryBy(smoothDat~time*elecID*wasRisky*locName,data=smoothRisk,FUN=c(mean))
names(regWof2)=c('time','elecID','wasRisky','locName','psc')
regWof2$elecID=as.factor(regWof2$elecID)
str(regWof2)

head(regWof2)
bothRisk3=summaryBy(psc~time*locName*wasRisky,data=regWof2,FUN=c(nanmean,se))
names(bothRisk3)=c('time','locName','wasRisky','psc','se')
head(bothRisk3)


bothRisk3$wasRisky=as.factor(bothRisk3$wasRisky)
bothRisk3b=bothRisk3[bothRisk3$time>1000 & bothRisk3$time<=3000,]
riskChoicePlot=ggplot() + theme_bw() + geom_line(data=bothRisk3b,aes(x=time,y=psc,color=locName,linetype=wasRisky),lwd=2) + geom_ribbon(data=bothRisk3b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black') + scale_y_continuous(limits=c(-4,8),breaks=c(-4,0,4,8)) + theme(panel.grid.minor = element_blank())
riskChoicePlot
plotName=sprintf('%s/hippMacro_valueDM_byRisk_choice.pdf',outDir)
ggsave(plotName,height=7,width=10,riskChoicePlot)

#Run stats on risk modulation
choicePeriod=regWof2[(regWof2$time>1700) & (regWof2$time<2200),]
choiceSum1=summaryBy(psc~elecID*locName*wasRisky,data=choicePeriod,FUN=c(mean))
choiceSum1$timeperiod='choice'

choiceSum2=reshape2::dcast(data=choiceSum1,elecID+locName+timeperiod~wasRisky,value.var="psc.mean")
names(choiceSum2)=c('elecID','locName','timeperiod','safe','risky')
choiceSum2$diff=choiceSum2$risky-choiceSum2$safe
choiceSum2$Subject=substr(choiceSum2$elecID,1,6)
choiceSum2$Subject=as.factor(choiceSum2$Subject)
head(choiceSum2)

rlmC1=lmer(diff~locName+(1|Subject),data=choiceSum2)
summary(rlmC1) #
confint(rlmC1,method='Wald')

pDat=choiceSum2[choiceSum2$locName=='postHipp',]
rlmCp=lmer(diff~1+(1|Subject),data=pDat)
summary(rlmCp)
confint(rlmCp,method='Wald')

aDat=choiceSum2[choiceSum2$locName=='antHipp',]
rlmCa=lmer(diff~1+(1|Subject),data=aDat)
summary(rlmCa)
confint(rlmCa,method='Wald')

############
##################
########################
## Delay
#######################
##################
############
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

#Load the data
outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/valueDM/hippocampus'

hippVal <- fread(sprintf('%s/hippMacro_valueBasedDM_respLocked.csv',baseDir), data.table = FALSE)
head(hippVal)

regWof1=summaryBy(psc~time*elecID*blockValue*locName,data=hippVal,FUN=c(mean))
str(regWof1)

regWof2=summaryBy(psc.mean~time*elecID*locName,data=regWof1,FUN=c(mean))
names(regWof2)=c('time','elecID','locName','psc')
regWof2$elecID=as.factor(regWof2$elecID)
str(regWof2)

#Stats for response against baseline during response window
delayPeriod=regWof2[(regWof2$time>1500) & (regWof2$time<=2000),]
delaySum2=summaryBy(psc~elecID*locName,data=delayPeriod,FUN=c(mean))
names(delaySum2)=c('elecID','locName','psc')
delaySum2$Subject=substr(delaySum2$elecID,1,6)
delaySum2$Subject=as.factor(delaySum2$Subject)
delaySum2$locName=as.factor(delaySum2$locName)
str(delaySum2)

lmD1=lmer(psc~locName+(1|Subject),data=delaySum2)
summary(lmD1)
confint(lmD1,method='Wald')

pDat=delaySum2[delaySum2$locName=='postHipp',]
lmDp=lmer(psc~1+(1|Subject),data=pDat)
summary(lmDp) 
confint(lmDp,method='Wald')

aDat=delaySum2[delaySum2$locName=='antHipp',]
lmDa=lmer(psc~1+(1|Subject),data=aDat)
summary(lmDa) 
confint(lmDa,method='Wald')


############
############
### Risk Modulation Hippocampus Response
############
############
to_keep <- c("outDir", "baseDir","nanmean","se","hippVal")
rm(list = setdiff(ls(), to_keep))

regWof1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName,data=hippVal,FUN=c(mean))
str(regWof1)

#Smooth data for plotting
smoothRisk=data.frame()
smoothRisk <- do.call(rbind, lapply(split(regWof1,with(regWof1, paste(elecID, blockValue, wasRisky, sep = "_"))),function(chunk) {
  chunk <- chunk[order(chunk$time), ]
  chunk$smoothDat <- ksmooth(chunk$time, chunk$psc.mean, kernel = "normal", bandwidth = 100)$y
  chunk
}))
#Make sure order is correct
smoothRisk <- smoothRisk[order(smoothRisk$elecID, 
                               smoothRisk$blockValue, 
                               smoothRisk$wasRisky), ]
rownames(smoothRisk) <- NULL

regWof2=summaryBy(smoothDat~time*elecID*wasRisky*locName,data=smoothRisk,FUN=c(mean))
names(regWof2)=c('time','elecID','wasRisky','locName','psc')
regWof2$elecID=as.factor(regWof2$elecID)
str(regWof2)

head(regWof2)
riskDat=summaryBy(psc~time*locName*wasRisky,data=regWof2,FUN=c(nanmean,se))
names(riskDat)=c('time','locName','wasRisky','psc','se')
head(riskDat)

riskDat$wasRisky=as.factor(riskDat$wasRisky)
riskRespPlot=ggplot() + theme_bw() + geom_line(data=riskDat,aes(x=time,y=psc,color=locName,linetype=wasRisky),lwd=2) + geom_ribbon(data=riskDat,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black') + scale_y_continuous(limits=c(-4,8),breaks=c(-4,0,4,8)) + theme(panel.grid.minor = element_blank())
riskRespPlot


riskDatb=riskDat[riskDat$time>500 & riskDat$time<=2500,]
riskRespPlot=ggplot() + theme_bw() + geom_line(data=riskDatb,aes(x=time,y=psc,color=locName,linetype=wasRisky),lwd=2) + geom_ribbon(data=riskDatb,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black') + scale_y_continuous(limits=c(-4,8),breaks=c(-4,0,4,8)) + theme(panel.grid.minor = element_blank())
riskRespPlot
plotName=sprintf('%s/hippMacro_valueDM_byRisk_response.pdf',outDir)
ggsave(plotName,height=7,width=10,riskRespPlot)

#Now also keep the average difference for barplots
respPeriod=regWof2[(regWof2$time>1500) & (regWof2$time<=2000),]
respSum1=summaryBy(psc~elecID*locName*wasRisky,data=respPeriod,FUN=c(mean))
respSum1$timeperiod='delay'

respSum2=reshape2::dcast(data=respSum1,elecID+locName+timeperiod~wasRisky,value.var="psc.mean")
names(respSum2)=c('elecID','locName','timeperiod','safe','risky')
respSum2$diff=respSum2$risky-respSum2$safe
respSum2$Subject=substr(respSum2$elecID,1,6)
respSum2$Subject=as.factor(respSum2$Subject)
head(respSum2)

lmD1r=lmer(diff~locName+(1|Subject),data=respSum2)
summary(lmD1r)
confint(lmD1r,method='Wald')

pDat=respSum2[respSum2$locName=='postHipp',]
lmDpr=lmer(diff~1+(1|Subject),data=pDat)
summary(lmDpr)
confint(lmDpr,method='Wald')

aDat=respSum2[respSum2$locName=='antHipp',]
lmDar=lmer(diff~1+(1|Subject),data=aDat)
summary(lmDar) #  No sig diff (p=.252)
confint(lmDar,method='Wald')


############
##################
########################
## Feedback
#######################
##################
############
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

#Load the data
outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/valueDM/hippocampus'

hippVal <- fread(sprintf('%s/hippMacro_valueBasedDM_feedbackLocked.csv',baseDir), data.table = FALSE)
head(hippVal)

regWof1=summaryBy(psc~time*elecID*blockValue*locName,data=hippVal,FUN=c(mean))
str(regWof1)

regWof2=summaryBy(psc.mean~time*elecID*locName,data=regWof1,FUN=c(mean))
names(regWof2)=c('time','elecID','locName','psc')
regWof2$elecID=as.factor(regWof2$elecID)
str(regWof2)

feedPeriod=regWof2[(regWof2$time>1500) & (regWof2$time<=2000),]
feedSum2=summaryBy(psc~elecID*locName,data=feedPeriod,FUN=c(mean))
feedSum2$Subject=substr(feedSum2$elecID,1,6)
feedSum2$Subject=as.factor(feedSum2$Subject)
feedSum2$locName=as.factor(feedSum2$locName)
str(feedSum2)
names(feedSum2)=c('elecID','locName','psc','Subject')

lmF1=lmer(psc~locName+(1|Subject),data=feedSum2)
summary(lmF1)
confint(lmF1,method='Wald')

pDatF=feedSum2[feedSum2$locName=='postHipp',]
lmFp=lmer(psc~1+(1|Subject),data=pDatF)
summary(lmFp)
confint(lmFp,method='Wald')

aDatF=feedSum2[feedSum2$locName=='antHipp',]
lmFa=lmer(psc~1+(1|Subject),data=aDatF)
summary(lmFa)
confint(lmFa,method='Wald')

############
############
### Risk Modulation Hippocampus Response
############
############
to_keep <- c("outDir", "baseDir","nanmean","se","hippVal")
rm(list = setdiff(ls(), to_keep))

regWof1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName,data=hippVal,FUN=c(mean))
str(regWof1)

#Smooth for plotting
smoothRisk3=data.frame()
smoothRisk3 <- do.call(rbind, lapply(split(regWof1,with(regWof1, paste(elecID, blockValue, wasRisky, sep = "_"))),function(chunk) {
  chunk <- chunk[order(chunk$time), ]
  chunk$smoothDat <- ksmooth(chunk$time, chunk$psc.mean, kernel = "normal", bandwidth = 100)$y
  chunk
}))
#Make sure in correct order
smoothRisk3 <- smoothRisk3[order(smoothRisk3$elecID, 
                                 smoothRisk3$blockValue, 
                                 smoothRisk3$wasRisky), ]
rownames(smoothRisk3) <- NULL

regWof2=summaryBy(smoothDat~time*elecID*wasRisky*locName,data=smoothRisk3,FUN=c(mean))
names(regWof2)=c('time','elecID','wasRisky','locName','psc')
regWof2$elecID=as.factor(regWof2$elecID)
str(regWof2)

head(regWof2)
bothRisk3=summaryBy(psc~time*locName*wasRisky,data=regWof2,FUN=c(nanmean,se))
names(bothRisk3)=c('time','locName','wasRisky','psc','se')
bothRisk3$locName=as.factor(bothRisk3$locName)
bothRisk3$wasRisky=as.factor(bothRisk3$wasRisky)

bothRisk3b=bothRisk3[bothRisk3$time>500 & bothRisk3$time<=2500,]
riskFeedPlot=ggplot() + theme_bw() + geom_line(data=bothRisk3b,aes(x=time,y=psc,color=locName,linetype=wasRisky),lwd=2) + geom_ribbon(data=bothRisk3b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black') + scale_y_continuous(limits=c(-4,8),breaks=c(-4,0,4,8)) + theme(panel.grid.minor = element_blank())
riskFeedPlot
plotName=sprintf('%s/hippMacro_valueDM_byRisk_feedback.pdf',outDir)
ggsave(plotName,height=7,width=10,riskFeedPlot)

#Now stats for the feedback period for risk modulation
feedPeriod=regWof2[(regWof2$time>1500) & (regWof2$time<=2000),]
feedSum1=summaryBy(psc~elecID*locName*wasRisky,data=feedPeriod,FUN=c(mean))
feedSum1$timeperiod='feedback'

feedSum2=reshape2::dcast(data=feedSum1,elecID+locName+timeperiod~wasRisky,value.var="psc.mean")
names(feedSum2)=c('elecID','locName','timeperiod','safe','risky')
feedSum2$diff=feedSum2$risky-feedSum2$safe
feedSum2$Subject=substr(feedSum2$elecID,1,6)
feedSum2$Subject=as.factor(feedSum2$Subject)
head(feedSum2)

lmF1=lmer(diff~locName+(1|Subject),data=feedSum2)
summary(lmF1)
confint(lmF1,method='Wald')

pDatF=feedSum2[feedSum2$locName=='postHipp',]
lmFp=lmer(diff~1+(1|Subject),data=pDatF)
summary(lmFp)
confint(lmFp,method='Wald')

aDatF=feedSum2[feedSum2$locName=='antHipp',]
lmFa=lmer(diff~1+(1|Subject),data=aDatF)
summary(lmFa)
confint(lmFa,method='Wald')